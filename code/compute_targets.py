#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""D4 comparator probability-limit targets via the general C2 form (MC).

docs/METHODS.md, comparator targets: D4 (patient-level standardized-lognormal
outcomes) has no closed-form comparator limit because s2_p and theta_hat are
dependent under the realized-data scheme. tau2_* / T_DL are computed by Monte
Carlo over the STUDY-level (urn) distribution with documented MCSE:

  A = 1/V, V = s2_p * c, c = 1/nT + 1/nC   (realized patient-level statistics)
  mu_A  = E[A Theta_hat]/E[A]
  q_*   = E[A (Theta_hat - mu_A)^2] = E[A Th^2] - 2 mu_A E[A Th] + mu_A^2 E[A]
  tau2_* = max(0, (q_* - 1)/E[A])                       (C1)
  T_DL  = E[w* Theta_hat]/E[w*],  w* = 1/(V + tau2_*)   (C2)

Scenario: S2 (Theta = Delta + u, tau = 0.5), Delta = 0 reference (shift only).
T2 = T3 = Delta exactly for D4/S2 (E[Theta_hat | N] = Delta).

Heavy-tail caveat: A involves 1/s2_p with arm sizes >= 3, so high moments of A
are heavy; MCSEs are computed from batch variability over independent batches
and a second, independent run (different seed) is reported for stability.

Outputs: results/targets/d4_targets.json
"""
import json, math, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import methods as rf

BATCH, NBATCH = 400_000, 25          # 1e7 studies per pass
TAU = 0.5                             # S2


def run(seed_offset):
    rng = np.random.default_rng(np.random.SeedSequence(
        20260925, spawn_key=(rf.FAMILY_ID["targets"], 4, 0, seed_offset)))
    law = rf.LAWS["primary"]
    sums = dict(A=0.0, ATh=0.0, ATh2=0.0, w0=0.0, w0N=0.0, N=0.0, n=0)
    # second pass needs tau2_* first -> two passes over independent data:
    # pass 1 estimates tau2_*; pass 2 (independent) estimates T_DL given tau2_*.
    # We store per-batch sums for both passes.
    def draw_studies(b):
        N = np.rint(np.exp(law["mu0"] + law["sigma0"] * rng.standard_normal(b)))
        N = np.maximum(N, law["clamp"])
        q = np.where(rng.random(b) < 1/3, 1/3, np.where(rng.random(b) < 0.5, 1/2, 2/3))
        nT = np.rint(q * N).astype(np.int64)
        nC = (N - nT).astype(np.int64)
        Th = TAU * rng.standard_normal(b)          # Delta = 0
        def arm(nflat, theta):
            starts = np.concatenate(([0], np.cumsum(nflat)[:-1]))
            tot = int(nflat.sum())
            U = (np.exp(rng.standard_normal(tot)) - math.exp(0.5)) \
                / math.sqrt(math.e * (math.e - 1))
            if theta is not None:
                U = U + np.repeat(theta, nflat)
            s1 = np.add.reduceat(U, starts)
            s2 = np.add.reduceat(U * U, starts)
            mean = s1 / nflat
            var = (s2 - nflat * mean ** 2) / (nflat - 1)
            return mean, var
        mC, sC2 = arm(nC, None)
        mT, sT2 = arm(nT, Th)
        s2p = ((nT - 1) * sT2 + (nC - 1) * sC2) / (N - 2)
        th = mT - mC
        c = 1.0 / nT + 1.0 / nC
        return N, th, s2p * c

    batches = []
    for _ in range(NBATCH):
        N, th, V = draw_studies(BATCH)
        A = 1.0 / V
        batches.append(dict(N=N, th=th, A=A, V=V))
    E_A = float(np.mean([b["A"].mean() for b in batches]))
    E_ATh = float(np.mean([b["A"].dot(b["th"]) / len(b["A"]) for b in batches]))
    E_ATh2 = float(np.mean([(b["A"] * b["th"] ** 2).mean() for b in batches]))
    mu_A = E_ATh / E_A
    q_star = E_ATh2 - 2 * mu_A * E_ATh + mu_A ** 2 * E_A
    tau2_star = max(0.0, (q_star - 1.0) / E_A)
    # pass 2: T_DL with tau2_* fixed, independent data
    rng2 = np.random.default_rng(np.random.SeedSequence(
        20260925, spawn_key=(rf.FAMILY_ID["targets"], 4, 1, seed_offset)))
    sw = swn = sn = cnt = 0.0
    tdl_batch = []
    for _ in range(NBATCH):
        N, th, V = draw_studies(BATCH)
        w = 1.0 / (V + tau2_star)
        sw += w.sum(); swn += (w * th).sum(); sn += N.sum(); cnt += w.size
        tdl_batch.append(float((w * th).sum() / w.sum()))
    T_DL = (swn / sw)            # E[w* th]/E[w*]
    E_N = sn / cnt
    # batch MCSEs
    def bmcse(vals):
        return float(np.std(vals, ddof=1) / math.sqrt(len(vals)))
    E_A_b = [b["A"].mean() for b in batches]
    E_ATh_b = [b["A"].dot(b["th"]) / len(b["A"]) for b in batches]
    E_ATh2_b = [(b["A"] * b["th"] ** 2).mean() for b in batches]
    return dict(E_A=E_A, E_ATh=E_ATh, E_ATh2=E_ATh2, mu_A=mu_A, q_star=q_star,
                tau2_star=tau2_star, T_DL=float(T_DL), E_N=float(E_N),
                T_DL_minus_Delta=float(T_DL),
                mcse=dict(E_A=bmcse(E_A_b), E_ATh=bmcse(E_ATh_b), E_ATh2=bmcse(E_ATh2_b),
                          T_DL=bmcse(tdl_batch)),
                n_studies=BATCH * NBATCH)


def main():
    t0 = time.time()
    r1 = run(0)
    r2 = run(9)   # independent stability run
    out = {"meta": {"script": "compute_targets.py", "date": time.strftime("%Y-%m-%d"),
                    "method": "general C2 by urn-level MC; S2, tau=0.5, Delta=0",
                    "batches": NBATCH, "batch_size": BATCH,
                    "seconds": time.time() - t0,
                    "note": "T2 = T3 = Delta exactly under D4/S2; T_DL - Delta is "
                            "the reported comparator-target shift. E[Theta_hat] is 0 "
                            "marginally only up to MC error; mu_A centers it."},
           "run1": r1, "run2": r2,
           "stability": dict(tau2_star_diff=r1["tau2_star"] - r2["tau2_star"],
                             T_DL_diff=r1["T_DL_minus_Delta"] - r2["T_DL_minus_Delta"])}
    # T_DL via ratio estimate: report also the ratio-of-means form
    with open(os.path.join(HERE, "..", "results", "targets", "d4_targets.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({k: out[k] for k in ("run1", "stability")}, indent=1))


if __name__ == "__main__":
    main()
