#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Shared DGM / estimator / design module for the ratio mean-difference
meta-analysis simulation study.

Single source of truth for run_simulations.py, compute_targets.py and
analyze_results.py. Implements exactly:

  - the simulation design documented in config/simulation_design.yaml;
  - the seed allocation documented in config/seeds.md
    (SeedSequence(20260925, spawn_key=(family, scen, dcode, M)));
  - estimation procedures P1-P4 and comparators C1-C3 documented in
    docs/METHODS.md.

Draw order per cell (fixed; see config/seeds.md):
  N -> q -> s2p -> (u) -> within-study noise.
  D4 replaces the s2p draw and the arm-mean noise draws by patient-level
  draws (s2p and the arm means are REALIZED from the patient data); its
  consumed order is N -> q -> (u) -> patient noise (control arm block
  first, then treatment arm block, per chunk of studies).
"""
import math
import numpy as np
from scipy import stats

MASTER_SEED = 20260925
Z975 = stats.norm.ppf(0.975)

# ------------------------------------------------------------------ laws ----
# actual-centered moments; provenance: results/targets/exact_targets.json
LAWS = {
    "primary": dict(mu0=3.0, sigma0=0.75, clamp=10.0,
                    mu_c=27.132678202789585, sd_c=22.68795062543078,
                    t2mt3=0.08361854460463167),
    "D2": dict(mu0=3.0, sigma0=0.5, clamp=10.0,
               mu_c=22.91262808847064, sd_c=11.95106165139892,
               t2mt3=0.05215927917676345),
    "clamp5": dict(mu0=3.0, sigma0=0.75, clamp=5.0,
                   mu_c=26.644973989079126, sd_c=23.088268241977325,
                   t2mt3=0.08665149476760768),
}
# scenario -> (tau, gamma) of the actual-centered effect law
SCEN = {1: (0.0, 0.0), 2: (0.5, 0.0), 3: (0.5, 0.1), 4: (1.0, 0.0), 5: (1.0, 0.1)}
DELTA_CODE = {0.0: 1, 0.5: 2}

# comparator probability-limit shifts T_DL - Delta (provenance:
# results/targets/comparator_targets.json; normal branch = corrected DGM,
# mixture = D1). D4 is computed by MC under the general C2 form
# (compute_targets.py -> results/targets/d4_targets.json).
T_DL_SHIFT = {
    "core":    {1: 0.0, 2: 0.0, 3: 0.018306227057702698, 4: 0.0, 5: 0.0073860233323441915},
    "D1":      {1: 0.0, 2: 0.0, 3: 0.01923982029783763, 4: 0.0, 5: 0.007520030420393065},
    "D2":      {1: 0.0, 2: 0.0, 3: 0.01585368268435662, 4: 0.0, 5: 0.006585806871547252},
    "D3":      {1: 0.0, 2: 0.0, 3: 0.018852497881809517, 4: 0.0, 5: 0.0076790497518302425},
    "D4":      None,   # computed by MC under the general C2 form (compute_targets.py)
    "clamp5":  {1: 0.0, 2: 0.0, 3: 0.019444654709078243, 4: 0.0, 5: 0.008529793582389896},
}

FAMILY_ID = {"core": 1, "D1": 2, "D2": 3, "D3": 4, "D4": 5, "clamp5": 6,
             "targets": 7, "verification": 8}


def cell_seed(family, scen, dcode, M):
    """Collision-free seed allocation (config/seeds.md)."""
    return np.random.SeedSequence(MASTER_SEED, spawn_key=(FAMILY_ID[family], scen, dcode, M))


# --------------------------------------------------------------- manifest ----
def build_manifest():
    """The 152-cell simulation manifest (config/simulation_design.yaml)."""
    Ms_core = [3, 5, 10, 30, 100, 300]
    Ms_sens = [5, 10, 30]
    cells = []
    for scen in range(1, 6):
        for d in (0.0, 0.5):
            for M in Ms_core:
                cells.append(dict(family="core", scen=scen, delta=d, M=M, R=100_000,
                                  law="primary", branch="normal", qdist="uniform"))
    for scen in range(1, 6):
        for d in (0.0, 0.5):
            for M in Ms_core:
                if M == 300:
                    continue
                cells.append(dict(family="D1", scen=scen, delta=d, M=M, R=20_000,
                                  law="primary", branch="mixture", qdist="uniform"))
    for scen in (1, 3):
        for d in (0.0, 0.5):
            for M in Ms_sens:
                cells.append(dict(family="D2", scen=scen, delta=d, M=M, R=20_000,
                                  law="D2", branch="normal", qdist="uniform"))
    for scen in (1, 3):
        for d in (0.0, 0.5):
            for M in Ms_sens:
                cells.append(dict(family="D3", scen=scen, delta=d, M=M, R=20_000,
                                  law="primary", branch="normal", qdist="d3"))
    for d in (0.0, 0.5):
        for M in Ms_sens:
            cells.append(dict(family="D4", scen=2, delta=d, M=M, R=20_000,
                              law="primary", branch="patient", qdist="uniform"))
    for scen in (1, 3):
        for d in (0.0, 0.5):
            for M in Ms_sens:
                cells.append(dict(family="clamp5", scen=scen, delta=d, M=M, R=20_000,
                                  law="clamp5", branch="normal", qdist="uniform"))
    for c in cells:
        c["dcode"] = DELTA_CODE[c["delta"]]
        c["key"] = (c["family"], c["scen"], c["dcode"], c["M"])
    return cells


def targets_for(cell, tdl_d4=None):
    """Population targets (T2, T3, T_DL) for a cell; Delta added."""
    fam, scen, d = cell["family"], cell["scen"], cell["delta"]
    law = LAWS[cell["law"]]
    t2mt3 = law["t2mt3"] if scen in (3, 5) else 0.0
    t2, t3 = d + t2mt3, d
    if fam == "D4":
        if tdl_d4 is None:
            raise ValueError("D4 requires the MC comparator target (tdl_d4)")
        tdl = d + tdl_d4
    else:
        tdl = d + T_DL_SHIFT[fam][scen]
    return dict(T2=t2, T3=t3, T_DL=tdl)


# -------------------------------------------------------------- generators ----
def _q_draw(rng, shape, qdist):
    if qdist == "uniform":
        return rng.choice([1/3, 1/2, 2/3], size=shape)
    if qdist == "d3":  # allocation sensitivity: P(q=1/3)=0.6, P(q=2/3)=0.4
        return np.where(rng.random(shape) < 0.6, 1/3, 2/3)
    raise ValueError(qdist)


def draw_cell(rng, R, M, cell, tau=None, gamma=None):
    """Summary-statistic DGMs: normal (corrected primary), mixture (D1)."""
    law = LAWS[cell["law"]]
    if tau is None:
        tau, gamma = SCEN[cell["scen"]]
    N = np.rint(np.exp(law["mu0"] + law["sigma0"] * rng.standard_normal((R, M))))
    N = np.maximum(N, law["clamp"])
    q = _q_draw(rng, (R, M), cell["qdist"])
    nT = np.rint(q * N)
    nC = N - nT
    nu = np.maximum(N - 2.0, 1.0)
    s2p = rng.chisquare(nu) / nu
    if cell["scen"] == 1:
        Theta = np.full((R, M), cell["delta"])
    elif cell["scen"] in (2, 4):
        Theta = cell["delta"] + tau * rng.standard_normal((R, M))
    else:
        Theta = (cell["delta"] + gamma * (N - law["mu_c"]) / law["sd_c"]
                 + tau * rng.standard_normal((R, M)))
    c = 1.0 / nT + 1.0 / nC
    if cell["branch"] == "normal":
        # corrected primary DGM: arm means normal, independent of s2p
        th = Theta + np.sqrt(c) * rng.standard_normal((R, M))
    elif cell["branch"] == "mixture":
        # D1 variance-coupled normal scale mixture (original scheme)
        th = Theta + np.sqrt(s2p * c) * rng.standard_normal((R, M))
    else:
        raise ValueError(cell["branch"])
    sj2 = s2p * c
    return N, th, sj2


def draw_cell_d4(rng, R, M, cell, tau=None):
    """D4 patient-level DGM: realized arm statistics."""
    law = LAWS[cell["law"]]
    if tau is None:
        tau, _ = SCEN[cell["scen"]]
    N = np.rint(np.exp(law["mu0"] + law["sigma0"] * rng.standard_normal((R, M))))
    N = np.maximum(N, law["clamp"])
    q = _q_draw(rng, (R, M), cell["qdist"])
    nT = np.rint(q * N)
    nC = N - nT
    if cell["scen"] in (2, 4):
        Theta = cell["delta"] + tau * rng.standard_normal((R, M))
    elif cell["scen"] == 1:
        Theta = np.full((R, M), cell["delta"])
    else:
        Theta = (cell["delta"] + SCEN[cell["scen"]][1] * (N - law["mu_c"]) / law["sd_c"]
                 + tau * rng.standard_normal((R, M)))
    NT = nT.ravel().astype(np.int64)
    NC = nC.ravel().astype(np.int64)
    TH = Theta.ravel()

    def arm_stats(nflat, theta=None):
        starts = np.concatenate(([0], np.cumsum(nflat)[:-1]))
        tot = int(nflat.sum())
        Z = rng.standard_normal(tot)
        U = (np.exp(Z) - math.exp(0.5)) / math.sqrt(math.e * (math.e - 1))
        if theta is not None:
            U = U + np.repeat(theta, nflat)
        s1 = np.add.reduceat(U, starts)
        s2 = np.add.reduceat(U * U, starts)
        mean = s1 / nflat
        var = (s2 - nflat * mean ** 2) / (nflat - 1)
        return mean, var

    mC, sC2 = arm_stats(NC, None)
    mT, sT2 = arm_stats(NT, TH)
    s2p = ((NT - 1) * sT2 + (NC - 1) * sC2) / (NT + NC - 2)
    th = (mT - mC).reshape(R, M)
    sj2 = (s2p * (1.0 / NT + 1.0 / NC)).reshape(R, M)
    return N, th, sj2


# --------------------------------------------------------------- estimators ----
def corr_rows(A, B):
    Am = A - A.mean(axis=1, keepdims=True)
    Bm = B - B.mean(axis=1, keepdims=True)
    num = (Am * Bm).sum(axis=1)
    den = np.sqrt((Am ** 2).sum(axis=1) * (Bm ** 2).sum(axis=1))
    return np.where(den > 0, num / np.maximum(den, 1e-300), 0.0)


def estimators(N, th, sj2):
    """All six procedures on the same realized data. Returns dict of 1-d arrays.

    P1: published Steps 1-5; P2: P1 with df M-1; P3: P1 with the -2*rho term
    dropped; P4: center R_hat - Khat/M (Cochran plug-in
    Khat = R*sz^2/mz^2 - rho*sy*sz/mz^2) retaining P1's V and t_{M-2};
    C1: IVW-DL (z); C2: conventional Hartung-Knapp (t_{M-1}); C3: equal-weight.
    """
    M = N.shape[1]
    Y = N * th
    ybar = Y.mean(axis=1)
    zbar = N.mean(axis=1)
    sy = Y.std(axis=1, ddof=1)
    sz = N.std(axis=1, ddof=1)
    rho = corr_rows(Y, N)
    Rhat = Y.sum(axis=1) / N.sum(axis=1)
    V2 = (sy ** 2 / zbar ** 2 + (ybar * sz / zbar ** 2) ** 2
          - 2.0 * rho * ybar * sy * sz / zbar ** 3) / M
    V = np.sqrt(np.maximum(V2, 0.0))
    Vn2 = (sy ** 2 / zbar ** 2 + (ybar * sz / zbar ** 2) ** 2) / M
    Vn = np.sqrt(np.maximum(Vn2, 0.0))
    Khat = Rhat * sz ** 2 / zbar ** 2 - rho * sy * sz / zbar ** 2
    corr = Khat / M
    t_m2 = stats.t.ppf(0.975, M - 2)
    t_m1 = stats.t.ppf(0.975, M - 1)
    out = dict(
        Rhat=Rhat, V=V, rho=rho, corr=corr,
        p1_lo=Rhat - t_m2 * V, p1_hi=Rhat + t_m2 * V,          # P1 (df M-2)
        p2_lo=Rhat - t_m1 * V, p2_hi=Rhat + t_m1 * V,          # P2 (df M-1)
        p3_lo=Rhat - t_m2 * Vn, p3_hi=Rhat + t_m2 * Vn,        # P3 (no -2 rho)
        p4_c=Rhat - corr,
        p4_lo=Rhat - corr - t_m2 * V, p4_hi=Rhat - corr + t_m2 * V,
    )
    # C1/C2
    w = 1.0 / sj2
    tw = w.sum(axis=1)
    hat_w = (w * th).sum(axis=1) / tw
    Q = (w * (th - hat_w[:, None]) ** 2).sum(axis=1)
    cden = tw - (w ** 2).sum(axis=1) / tw
    tau2 = np.maximum(0.0, (Q - (M - 1)) / np.maximum(cden, 1e-300))
    wr = 1.0 / (sj2 + tau2[:, None])
    twr = wr.sum(axis=1)
    hat_DL = (wr * th).sum(axis=1) / twr
    se_DL = np.sqrt(1.0 / twr)
    se_HK = np.sqrt((wr * (th - hat_DL[:, None]) ** 2).sum(axis=1) / ((M - 1) * twr))
    hat_eq = th.mean(axis=1)
    se_eq = th.std(axis=1, ddof=1) / math.sqrt(M)
    out.update(dl=hat_DL, se_dl=se_DL, se_hk=se_HK, tau2=tau2, eq=hat_eq, se_eq=se_eq,
               c1_lo=hat_DL - Z975 * se_DL, c1_hi=hat_DL + Z975 * se_DL,
               c2_lo=hat_DL - t_m1 * se_HK, c2_hi=hat_DL + t_m1 * se_HK,
               c3_lo=hat_eq - t_m1 * se_eq, c3_hi=hat_eq + t_m1 * se_eq)
    return out
