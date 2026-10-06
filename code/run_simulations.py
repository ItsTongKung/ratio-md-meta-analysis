#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Full-study simulation runner (152 cells; base grid + boundary-precision
boost).

Executes the simulation grid documented in config/simulation_design.yaml with
the six procedures of docs/METHODS.md (P1-P4, C1-C3), one shared realized
dataset per replication, the deterministic seeds of config/seeds.md, and
memory-bounded chunks.

Per-cell streaming aggregates (per-replication raw data are deterministically
regenerable from the recorded seed key + draw order; see config/seeds.md):
  every method: mean/sd of the point estimate, coverage vs own target / T2 /
  T3 / T_DL (each with MCSE), mean width, rejection vs 0 and vs T3;
  P1: variance-calibration block (E[V], E[V^2], empirical SD and variance of
  R_hat, mean rho-hat); C1: mean tau2-hat; P4: paired P4-P1 differences
  (bias, MSE, coverage - paired MCSEs from within-replication differences).

Outputs one JSON per cell under results/cell_summaries/.

CLI:
  python code/run_simulations.py --family core            # one family
  python code/run_simulations.py --family all             # full base grid
  python code/run_simulations.py --family all --boost "<override list>"
  python code/run_simulations.py --family core --dump-reps 5 --dump-cells "core,S3,1,10;core,S1,2,100;D1,S2,1,5"

The boost override list has the form 'family,S<scen>,d<dcode>,M=R;...' and
re-executes the named cells at the given total R (same seed stream: the first
R_base replications are bit-identical and the boost extends the stream). The
executed boost list is documented in config/simulation_design.yaml and
config/seeds.md; re-running the base pass followed by the boost pass
reproduces the archived cell summaries bit-for-bit (except wall time).
"""
import argparse, json, math, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import methods as rf

LOHI = {"P1": ("p1_lo", "p1_hi"), "P2": ("p2_lo", "p2_hi"), "P3": ("p3_lo", "p3_hi"),
        "P4": ("p4_lo", "p4_hi"), "C1": ("c1_lo", "c1_hi"), "C2": ("c2_lo", "c2_hi"),
        "C3": ("c3_lo", "c3_hi")}
CENTER = {"P1": "Rhat", "P2": "Rhat", "P3": "Rhat", "P4": "p4_c",
          "C1": "dl", "C2": "dl", "C3": "eq"}
OWN = {"P1": "T2", "P2": "T2", "P3": "T2", "P4": "T2",
       "C1": "T_DL", "C2": "T_DL", "C3": "T3"}
METHODS = list(LOHI)


def new_acc():
    return dict(n=0, sx=0.0, sxx=0.0, cov_T2=0, cov_T3=0, cov_T_DL=0,
                w=0.0, rej0=0, rejT3=0)


def acc_update(a, x, lo, hi, tgts, R):
    a["n"] += R
    a["sx"] += float(x.sum())
    a["sxx"] += float((x * x).sum())
    for t, v in tgts.items():
        a[f"cov_{t}"] += int(np.sum((lo <= v) & (hi >= v)))
    a["w"] += float((hi - lo).sum())
    a["rej0"] += int(np.sum((lo > 0) | (hi < 0)))
    t3 = tgts["T3"]
    a["rejT3"] += int(np.sum((lo > t3) | (hi < t3)))


def finalize_acc(a, R):
    n = a["n"]
    mean = a["sx"] / n
    var = max(a["sxx"] / n - mean ** 2, 0.0)
    sd = math.sqrt(var)
    out = dict(mean_est=mean, emp_var=var, emp_sd=sd, mcse_est=sd / math.sqrt(n),
               mean_width=a["w"] / n)
    for t in ("T2", "T3", "T_DL"):
        p = a[f"cov_{t}"] / n
        out[f"cov_{t}"] = p
        out[f"cov_{t}_mcse"] = math.sqrt(max(p * (1 - p), 1e-300) / n)
    for k in ("rej0", "rejT3"):
        p = a[k] / n
        out[k] = p
        out[f"{k}_mcse"] = math.sqrt(max(p * (1 - p), 1e-300) / n)
    return out


def run_cell(cell, chunk, tdl_d4=None, dump_reps=0, R_override=None):
    R = R_override or cell["R"]
    targets = rf.targets_for(cell, tdl_d4=tdl_d4)
    accs = {m: new_acc() for m in METHODS}
    extra = dict(P1=dict(sv=0.0, sv2=0.0, srho=0.0),
                 C1=dict(stau=0.0),
                 P4=dict(scorr=0.0, scorr2=0.0, sdd=0.0, sdd2=0.0,
                         sdmse=0.0, sdmse2=0.0, sdcov=0, sdcov2=0))
    rng = np.random.default_rng(rf.cell_seed(*cell["key"]))
    d4 = cell["branch"] == "patient"
    ch = 5_000 if d4 else chunk
    nch = math.ceil(R / ch)
    t0 = time.time()
    dumps = []
    done = 0
    for _ in range(nch):
        r = min(ch, R - done)
        if d4:
            N, th, sj2 = rf.draw_cell_d4(rng, r, cell["M"], cell)
        else:
            N, th, sj2 = rf.draw_cell(rng, r, cell["M"], cell)
        e = rf.estimators(N, th, sj2)
        if dump_reps and len(dumps) < 1:
            k = min(dump_reps, r)
            dumps.append(dict(N=N[:k].tolist(), th=th[:k].tolist(),
                              sj2=sj2[:k].tolist(),
                              Rhat=e["Rhat"][:k].tolist(), V=e["V"][:k].tolist(),
                              p1_lo=e["p1_lo"][:k].tolist(), p1_hi=e["p1_hi"][:k].tolist()))
        tgts = {t: targets[t] for t in ("T2", "T3", "T_DL")}
        for m in METHODS:
            lo, hi = e[LOHI[m][0]], e[LOHI[m][1]]
            acc_update(accs[m], e[CENTER[m]], lo, hi, tgts, r)
        extra["P1"]["sv"] += float(e["V"].sum())
        extra["P1"]["sv2"] += float((e["V"] ** 2).sum())
        extra["P1"]["srho"] += float(e["rho"].sum())
        extra["C1"]["stau"] += float(e["tau2"].sum())
        p4 = extra["P4"]
        dd = e["p4_c"] - e["Rhat"]
        dmse = (e["p4_c"] - targets["T2"]) ** 2 - (e["Rhat"] - targets["T2"]) ** 2
        dcov = ((e["p4_lo"] <= targets["T2"]) & (e["p4_hi"] >= targets["T2"])).astype(np.int64) \
             - ((e["p1_lo"] <= targets["T2"]) & (e["p1_hi"] >= targets["T2"])).astype(np.int64)
        p4["scorr"] += float(e["corr"].sum())
        p4["scorr2"] += float((e["corr"] ** 2).sum())
        p4["sdd"] += float(dd.sum())
        p4["sdd2"] += float((dd * dd).sum())
        p4["sdmse"] += float(dmse.sum())
        p4["sdmse2"] += float((dmse * dmse).sum())
        p4["sdcov"] += int(dcov.sum())
        p4["sdcov2"] += int((dcov * dcov).sum())
        done += r
    out = dict(cell=dict(cell), seed_key=list(cell["key"]), targets=targets,
               R_executed=R, chunks=nch, chunk_size=ch, methods={},
               draw_order="N -> q -> s2p -> (u) -> within-study noise"
                          + (" [D4: N -> q -> (u) -> patient draws, control block "
                             "then treatment block; s2p and arm means realized]"
                             if d4 else ""),
               wall_s=time.time() - t0)
    for m in METHODS:
        out["methods"][m] = finalize_acc(accs[m], R)
        out["methods"][m]["own_target"] = OWN[m]
        out["methods"][m]["own_cov"] = out["methods"][m][f"cov_{OWN[m]}"]
        out["methods"][m]["own_cov_mcse"] = out["methods"][m][f"cov_{OWN[m]}_mcse"]
    n = R
    out["variance_calibration"] = {
        "mean_V": extra["P1"]["sv"] / n, "mean_V2": extra["P1"]["sv2"] / n,
        "emp_sd_Rhat": out["methods"]["P1"]["emp_sd"],
        "emp_var_Rhat": out["methods"]["P1"]["emp_var"],
        "mean_rho": extra["P1"]["srho"] / n}
    out["c1_tau2_mean"] = extra["C1"]["stau"] / n
    p4 = extra["P4"]
    sd_dd = math.sqrt(max(p4["sdd2"] / n - (p4["sdd"] / n) ** 2, 0.0))
    sd_dmse = math.sqrt(max(p4["sdmse2"] / n - (p4["sdmse"] / n) ** 2, 0.0))
    mean_dcov = p4["sdcov"] / n
    sd_dcov = math.sqrt(max(p4["sdcov2"] / n - mean_dcov ** 2, 0.0))
    out["p4_paired"] = {
        "mean_correction_Khat_over_M": p4["scorr"] / n,
        "sd_correction": math.sqrt(max(p4["scorr2"] / n - (p4["scorr"] / n) ** 2, 0.0)),
        "bias_diff_P4_minus_P1": p4["sdd"] / n,
        "bias_diff_mcse_paired": sd_dd / math.sqrt(n),
        "mse_diff_P4_minus_P1": p4["sdmse"] / n,
        "mse_diff_mcse_paired": sd_dmse / math.sqrt(n),
        "cov_diff_P4_minus_P1_vs_T2": mean_dcov,
        "cov_diff_mcse_paired": sd_dcov / math.sqrt(n)}
    if dumps:
        out["replication_dump"] = dumps[0]
    return out


def cell_id(family, scen, dcode, M):
    return f"{family}_S{scen}_d{dcode}_M{M}"


def parse_key(k):
    fam, scen, dcode, M = k.split(",")
    return (fam, int(scen.lstrip("sS")), int(dcode), int(M))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--family", default="all")
    ap.add_argument("--chunk", type=int, default=10_000)
    ap.add_argument("--boost", default="", help="'family,Sscen,dcode,M=R' comma list")
    ap.add_argument("--dump-reps", type=int, default=0)
    ap.add_argument("--dump-cells", default="")
    args = ap.parse_args()
    man = rf.build_manifest()
    boosts = {}
    if args.boost:
        for item in args.boost.split(";"):
            if not item:
                continue
            k, v = item.rsplit("=", 1)
            boosts[parse_key(k)] = int(v)
    dump_cells = set()
    if args.dump_cells:
        for k in args.dump_cells.split(";"):
            dump_cells.add(parse_key(k))
    with open(os.path.join(HERE, "..", "results", "targets", "d4_targets.json")) as f:
        tdl_d4 = json.load(f)["run1"]["T_DL_minus_Delta"]
    fams = [args.family] if args.family != "all" else \
        ["core", "D1", "D2", "D3", "D4", "clamp5"]
    t0 = time.time()
    for cell in man:
        if cell["family"] not in fams:
            continue
        ro = boosts.get(cell["key"])
        cid = cell_id(*cell["key"])
        res = run_cell(cell, args.chunk, tdl_d4=tdl_d4,
                       dump_reps=args.dump_reps if cell["key"] in dump_cells else 0,
                       R_override=ro)
        path = os.path.join(HERE, "..", "results", "cell_summaries", cid + ".json")
        with open(path, "w") as f:
            json.dump(res, f, indent=1)
        p1 = res["methods"]["P1"]
        print(f"{cid}: R={res['R_executed']} wall={res['wall_s']:.1f}s "
              f"P1 own_cov={p1['own_cov']:.4f}±{p1['own_cov_mcse']:.4f}", flush=True)
    print(f"family {args.family} done in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
