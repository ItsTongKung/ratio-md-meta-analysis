#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Analysis pipeline. Produces ALL analysis tables and the summary JSON
directly from results/cell_summaries/*.json plus the committed target files
(no hand-transcribed numbers; docs/METHODS.md).

Outputs (results/analysis/):
  master_results.csv        one row per cell x method (all metrics)
  H1_primary_coverage.csv   P1 coverage for T2, core grid + classification
  H3_failure_conditions.csv all cells/methods outside or near the boundary
  estimand_map.csv          targets vs observed mean estimates (per cell)
  H4_paired_P4_P1.csv       paired P4-P1 differences per cell
  variance_calibration.csv  the five variance-calibration quantities per cell
  df_and_covariance_probe.csv  P2/P3 vs P1 paired coverage differences
  comparators.csv           C1/C2/C3 coverage for their estimands (T_DL / T3)
                            and clinical targets
  sensitivity_summary.csv   D1/D2/D3/D4/clamp5 vs core (P1, C1, C2, C3)
  summary.json              headline aggregates used by the reports
"""
import csv, json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
WS = os.path.abspath(os.path.join(HERE, "..", ".."))
RES = os.path.join(HERE, "..", "results")
ANA = os.path.join(HERE, "..", "results", "analysis")
sys.path.insert(0, HERE)
import methods as rf

NOMINAL, MARGIN = 0.95, 0.03
LOW, HIGH = NOMINAL - MARGIN, NOMINAL + MARGIN


def classify(cov, mcse):
    if cov < LOW:
        d = LOW - cov
        return ("OUTSIDE_BELOW" if d > 3 * mcse else "BOUNDARY_BELOW"), d, mcse
    if cov > HIGH:
        d = cov - HIGH
        return ("OUTSIDE_ABOVE" if d > 3 * mcse else "BOUNDARY_ABOVE"), d, mcse
    return "WITHIN", min(cov - LOW, HIGH - cov), mcse


def load_cells():
    man = rf.build_manifest()
    cells = {}
    for c in man:
        cid = f"{c['family']}_S{c['scen']}_d{c['dcode']}_M{c['M']}"
        with open(os.path.join(RES, "cell_summaries", cid + ".json")) as f:
            cells[cid] = json.load(f)
    return man, cells


def main():
    os.makedirs(ANA, exist_ok=True)
    man, cells = load_cells()
    tdl_d4 = json.load(open(os.path.join(RES, "targets", "d4_targets.json")))["run1"]

    # ---------------------------------------------------------- master ----
    mrows = []
    for c in man:
        cid = f"{c['family']}_S{c['scen']}_d{c['dcode']}_M{c['M']}"
        r = cells[cid]
        for m, md in r["methods"].items():
            row = dict(cell=cid, family=c["family"], scen=c["scen"], delta=c["delta"],
                       M=c["M"], R=r["R_executed"], method=m,
                       own_target=md["own_target"],
                       T2=r["targets"]["T2"], T3=r["targets"]["T3"], T_DL=r["targets"]["T_DL"])
            row.update({k: md[k] for k in ("mean_est", "mcse_est", "emp_var", "emp_sd",
                                           "mean_width", "rej0", "rej0_mcse",
                                           "rejT3", "rejT3_mcse")})
            for t in ("T2", "T3", "T_DL"):
                row[f"cov_{t}"] = md[f"cov_{t}"]
                row[f"cov_{t}_mcse"] = md[f"cov_{t}_mcse"]
            row["own_cov"] = md["own_cov"]
            row["own_cov_mcse"] = md["own_cov_mcse"]
            cls, dist, _ = classify(md["own_cov"], md["own_cov_mcse"])
            row["own_class"] = cls
            row["own_dist_to_boundary"] = dist
            mrows.append(row)

    def write_csv(name, rows, cols):
        with open(os.path.join(ANA, name), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)

    mcols = list(mrows[0].keys())
    write_csv("master_results.csv", mrows, mcols)

    # ------------------------------------------------------------- H1 ----
    h1 = [r for r in mrows if r["method"] == "P1" and r["family"] == "core"]
    write_csv("H1_primary_coverage.csv", h1, mcols)
    h1_counts = {}
    for r in h1:
        h1_counts[r["own_class"]] = h1_counts.get(r["own_class"], 0) + 1

    # ------------------------------------------------------------- H3 ----
    h3 = [r for r in mrows if r["own_class"] != "WITHIN"]
    h3.sort(key=lambda r: -abs(r["own_cov"] - NOMINAL))
    write_csv("H3_failure_conditions.csv", h3, mcols)

    # -------------------------------------------------------- estimand ----
    erows = []
    own_tgt_of = {"P1": "T2", "C1": "T_DL", "C2": "T_DL", "C3": "T3"}
    for r in mrows:
        if r["method"] in ("P1", "C1", "C2", "C3"):
            tgt = own_tgt_of[r["method"]]
            erows.append(dict(cell=r["cell"], family=r["family"], scen=r["scen"],
                              delta=r["delta"], M=r["M"], method=r["method"],
                              T2=r["T2"], T3=r["T3"], T_DL=r["T_DL"],
                              mean_est=r["mean_est"], mcse_est=r["mcse_est"],
                              bias_vs_own_target=r["mean_est"] - r[tgt]))
    write_csv("estimand_map.csv", erows, list(erows[0].keys()))

    # ------------------------------------------------------------- H4 ----
    h4rows = []
    for c in man:
        cid = f"{c['family']}_S{c['scen']}_d{c['dcode']}_M{c['M']}"
        r = cells[cid]
        p4 = r["p4_paired"]
        h4rows.append(dict(cell=cid, family=c["family"], scen=c["scen"],
                           delta=c["delta"], M=c["M"], R=r["R_executed"],
                           T2=r["targets"]["T2"], **p4))
    write_csv("H4_paired_P4_P1.csv", h4rows, list(h4rows[0].keys()))

    # ------------------------------------------------ variance calib ----
    vrows = []
    for c in man:
        cid = f"{c['family']}_S{c['scen']}_d{c['dcode']}_M{c['M']}"
        r = cells[cid]
        vc = r["variance_calibration"]
        vrows.append(dict(cell=cid, family=c["family"], scen=c["scen"],
                          delta=c["delta"], M=c["M"], R=r["R_executed"], **vc,
                          ratio_meanV_to_empSD=vc["mean_V"] / vc["emp_sd_Rhat"] if vc["emp_sd_Rhat"] > 0 else None,
                          ratio_meanV2_to_empVar=vc["mean_V2"] / vc["emp_var_Rhat"] if vc["emp_var_Rhat"] > 0 else None))
    write_csv("variance_calibration.csv", vrows, list(vrows[0].keys()))

    # ------------------------------------------------- P2/P3 probes ----
    prows = []
    for c in man:
        cid = f"{c['family']}_S{c['scen']}_d{c['dcode']}_M{c['M']}"
        r = cells[cid]
        for probe in ("P2", "P3"):
            p1c = r["methods"]["P1"]["cov_T2"]
            pc = r["methods"][probe]["cov_T2"]
            d = pc - p1c
            if probe == "P2":
                # Nesting: t_{M-2,.975} > t_{M-1,.975}, so the P1 (df=M-2)
                # interval is the WIDER one and P2 is nested INSIDE P1. The
                # paired difference I_P2 - I_P1 therefore takes values 0 or -1
                # with transition probability q = C_P1 - C_P2 in [0,1]; exact
                # paired variance q(1-q)/R.
                q = max(p1c - pc, 0.0)
                sdd = math.sqrt(max(q * (1 - q), 0.0))
            else:
                # P3 not nested in general -> Frechet conservative bound
                sdd = math.sqrt(max(p1c + pc - 2 * max(0.0, p1c + pc - 1.0), 0.0))
            prows.append(dict(cell=cid, family=c["family"], scen=c["scen"],
                              delta=c["delta"], M=c["M"],
                              cov_P1=p1c, cov_probe=pc, cov_diff=probe + "-P1",
                              diff=d, q_transition=(q if probe == "P2" else None),
                              diff_mcse_paired=sdd / math.sqrt(r["R_executed"]),
                              width_P1=r["methods"]["P1"]["mean_width"],
                              width_probe=r["methods"][probe]["mean_width"]))
    write_csv("df_and_covariance_probe.csv", prows, list(prows[0].keys()))

    # ---------------------------------------------------- comparators ----
    crows = [r for r in mrows if r["method"] in ("C1", "C2", "C3")]
    write_csv("comparators.csv", crows, mcols)

    # -------------------------------------------------- sensitivity ----
    srows = [r for r in mrows if r["method"] in ("P1", "C1", "C2", "C3")
             and r["family"] != "core"]
    write_csv("sensitivity_summary.csv", srows, mcols)

    # -------------------------------------------------------- summary ----
    def agg(rows, key):
        vals = [r[key] for r in rows]
        return {"min": min(vals), "max": max(vals), "n": len(vals)}

    core_p1 = h1
    by_scen = {}
    for sc in range(1, 6):
        by_scen[f"S{sc}"] = {}
        for d in (0.0, 0.5):
            sub = {r["M"]: r["own_cov"] for r in h1
                   if r["scen"] == sc and r["delta"] == d}
            by_scen[f"S{sc}"][f"d{d}"] = dict(sorted(sub.items()))
    summary = {
        "H1": {
            "counts": h1_counts,
            "own_cov_range_core_P1": agg(core_p1, "own_cov"),
            "outside_cells": [dict(cell=r["cell"], cov=r["own_cov"], mcse=r["own_cov_mcse"],
                                   cls=r["own_class"]) for r in h1 if r["own_class"] != "WITHIN"],
            "by_scenario": by_scen,
        },
        "H2": {
            "targets_core": {f"S{sc}": dict(T2_minus_T3=rf.LAWS["primary"]["t2mt3"] if sc in (3, 5) else 0.0,
                                            T_DL_minus_Delta=rf.T_DL_SHIFT["core"][sc])
                             for sc in range(1, 6)},
            "mean_est_diff_T3_vs_T2_core_M300": {
                f"S{sc}": {m: dict(est=cells[f"core_S{sc}_d1_M300"]["methods"][m]["mean_est"],
                                   T2=cells[f"core_S{sc}_d1_M300"]["targets"]["T2"],
                                   T3=cells[f"core_S{sc}_d1_M300"]["targets"]["T3"],
                                   T_DL=cells[f"core_S{sc}_d1_M300"]["targets"]["T_DL"])
                           for m in ("P1", "C1", "C3")} for sc in (3, 5)},
        },
        "H3": {
            "n_rows_nonwithin_all_methods": len(h3),
            "worst_own_cov_all_methods": min(mrows, key=lambda r: r["own_cov"])["own_cov"],
            "outside_cells_all_methods": [
                dict(cell=r["cell"], method=r["method"], cov=r["own_cov"],
                     mcse=r["own_cov_mcse"], cls=r["own_class"]) for r in h3
                if r["own_class"].startswith("OUTSIDE")][:60],
            "boundary_cells_all_methods": [
                dict(cell=r["cell"], method=r["method"], cov=r["own_cov"],
                     mcse=r["own_cov_mcse"], cls=r["own_class"]) for r in h3
                if r["own_class"].startswith("BOUNDARY")][:60],
        },
        "H4": {
            "mse_diff_sign_counts": {
                "P4_better_MSE": sum(1 for r in h4rows if r["mse_diff_P4_minus_P1"] < 0),
                "P4_worse_MSE": sum(1 for r in h4rows if r["mse_diff_P4_minus_P1"] > 0),
                "n_cells": len(h4rows)},
            "mse_diff_sig": [dict(cell=r["cell"], d=r["mse_diff_P4_minus_P1"],
                                  mcse=r["mse_diff_mcse_paired"])
                             for r in h4rows
                             if abs(r["mse_diff_P4_minus_P1"]) > 3 * r["mse_diff_mcse_paired"]],
            "cov_diff_sig": [dict(cell=r["cell"], d=r["cov_diff_P4_minus_P1_vs_T2"],
                                  mcse=r["cov_diff_mcse_paired"])
                             for r in h4rows
                             if abs(r["cov_diff_P4_minus_P1_vs_T2"]) > 3 * r["cov_diff_mcse_paired"]],
            "max_abs_bias_diff": max(abs(r["bias_diff_P4_minus_P1"]) for r in h4rows),
            "max_abs_sd_correction": max(abs(r["sd_correction"]) for r in h4rows),
            "max_abs_mean_correction": max(abs(r["mean_correction_Khat_over_M"]) for r in h4rows),
        },
        "D4_target": {"T_DL_minus_Delta": tdl_d4["T_DL_minus_Delta"],
                      "tau2_star": tdl_d4["tau2_star"],
                      "mcse_T_DL": tdl_d4["mcse"]["T_DL"]},
    }

    # ---- Deterministic full-grid additions (post-processing of the cell
    # JSONs, so every downstream number is pulled from data) ----
    p1_full = [r for r in mrows if r["method"] == "P1"]
    fg_counts = {}
    for r in p1_full:
        fg_counts[r["own_class"]] = fg_counts.get(r["own_class"], 0) + 1
    p2rows = [r for r in prows if r["cov_diff"] == "P2-P1"]
    closer = {"P1": 0, "P2": 0, "tie": 0}
    for r in p2rows:
        d1, d2 = abs(r["cov_P1"] - NOMINAL), abs(r["cov_probe"] - NOMINAL)
        closer["P1" if d1 < d2 else ("P2" if d2 < d1 else "tie")] += 1
    cls_out = {}
    for m in ("P1", "P2"):
        c2 = {}
        for r in mrows:
            if r["method"] == m:
                k = "OUTSIDE" if r["own_class"].startswith("OUTSIDE") else r["own_class"]
                c2[k] = c2.get(k, 0) + 1
        cls_out[m] = c2
    p4pos = sum(1 for r in h4rows if r["cov_diff_P4_minus_P1_vs_T2"] > 0)
    summary["H1_full_grid"] = {
        "counts": fg_counts,
        "n_cells": len(p1_full),
        "nonwithin_cells": [dict(cell=r["cell"], cov=r["own_cov"], mcse=r["own_cov_mcse"],
                                 cls=r["own_class"]) for r in p1_full if r["own_class"] != "WITHIN"],
    }
    summary["df_probe_P2_vs_P1"] = {
        "closer_to_nominal": closer,
        "decisive_OUTSIDE_counts": {m: cls_out[m].get("OUTSIDE", 0) for m in ("P1", "P2")},
        "classification_full": cls_out,
        "max_abs_paired_diff": max(abs(r["diff"]) for r in p2rows),
        "max_abs_paired_diff_cell": max(p2rows, key=lambda r: abs(r["diff"]))["cell"],
        "M3_cells_total": sum(1 for r in p2rows if r["M"] == 3),
        "M3_cells_P2_closer": sum(1 for r in p2rows if r["M"] == 3
                                  and abs(r["cov_probe"] - NOMINAL) < abs(r["cov_P1"] - NOMINAL)),
    }
    summary["P4_cov_diff_profile"] = {
        "positive_cells": p4pos,
        "max_positive": max((r["cov_diff_P4_minus_P1_vs_T2"] for r in h4rows), default=0.0),
        "significantly_positive": sum(1 for r in h4rows
                                      if r["cov_diff_P4_minus_P1_vs_T2"] > 3 * r["cov_diff_mcse_paired"]),
        "significantly_negative": sum(1 for r in h4rows
                                      if r["cov_diff_P4_minus_P1_vs_T2"] < -3 * r["cov_diff_mcse_paired"]),
    }
    acc = {"base_computation": 0, "boost_job_computation": 0,
           "unique_additional_beyond_base": 0, "final_unique_committed": 0,
           "boosted_cells": 0}
    for c in man:
        cid = f"{c['family']}_S{c['scen']}_d{c['dcode']}_M{c['M']}"
        Rb = 100000 if c["family"] == "core" else 20000
        R = cells[cid]["R_executed"]
        acc["base_computation"] += Rb
        acc["final_unique_committed"] += R
        if R != Rb:
            acc["boost_job_computation"] += R
            acc["unique_additional_beyond_base"] += R - Rb
            acc["boosted_cells"] += 1
    summary["replication_accounting"] = acc
    with open(os.path.join(ANA, "summary.json"), "w") as f:
        json.dump(summary, f, indent=1)
    print(json.dumps(summary["H1_full_grid"], indent=1))
    print(json.dumps(summary["df_probe_P2_vs_P1"], indent=1))
    print("analysis done")


if __name__ == "__main__":
    main()
