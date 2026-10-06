#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Manuscript table generation.

Builds the 11 article tables in outputs/tables/ (created when needed)
directly from results/analysis/*.csv, results/analysis/summary.json and
the M=300 cell summaries (no hand-transcribed numbers):

  Table1_core_scenarios            scenario design constants (config/simulation_design.yaml)
  Table2_coverage_classification   coverage tolerance band classification counts
  Table3_estimand_tracking_M300    large-M estimand tracking (core, Delta=0)
  Table4_secondary_probes          P4 / degrees-of-freedom / covariance probes
  TableS1_complete_master_results  master_results.csv (all 152 cells x 7 methods)
  TableS2_P1_all_152_cells         P1 rows of the master table
  TableS3_df_covariance_probes     df_and_covariance_probe.csv
  TableS4_P4_paired_results        H4_paired_P4_P1.csv
  TableS5_variance_calibration     variance_calibration.csv
  TableS6_estimand_map             estimand_map.csv
  TableS7_sensitivity_summary      sensitivity_summary.csv

Numeric formatting is deterministic: values are parsed with pandas'
default CSV float parser and printed with at most 16 decimal places
(trailing zeros stripped); integer columns stay integers;
scientific-notation values are kept as produced. Publication tables are
generated on demand and deliberately not version-controlled;
verify_reproduction.py checks their scientific content.

CLI:  python code/make_tables.py [--out DIR]
"""
import argparse, csv, json, os, sys
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
ANA = os.path.join(RES, "analysis")
OUT = os.path.join(HERE, "..", "outputs", "tables")

# Scenario design constants (config/simulation_design.yaml; prespecified design).
TABLE1_ROWS = [
    ("S1", 0.0, 0.0, "No",  "Common effect"),
    ("S2", 0.5, 0.0, "No",  "Moderate heterogeneity"),
    ("S3", 0.5, 0.1, "Yes", "Moderate heterogeneity + informative study size"),
    ("S4", 1.0, 0.0, "No",  "Strong heterogeneity"),
    ("S5", 1.0, 0.1, "Yes", "Strong heterogeneity + informative study size"),
]

TABLE2_SCOPES = [
    ("P1 core (60 cells)", "P1", lambda r: r["family"] == "core"),
    ("P1 full grid (152 cells)", "P1", lambda r: True),
    ("IVW-DL full grid", "C1", lambda r: True),
    ("HKSJ full grid", "C2", lambda r: True),
    ("Equal-weight full grid", "C3", lambda r: True),
]

TABLE3_ROWS = [  # (scenario, label, method, own target)
    ("S3", "P1 ratio", "P1", "T2"),
    ("S3", "IVW-DL / HKSJ center", "C1", "T_DL"),
    ("S3", "Equal weight", "C3", "T3"),
    ("S5", "P1 ratio", "P1", "T2"),
    ("S5", "IVW-DL / HKSJ center", "C1", "T_DL"),
    ("S5", "Equal weight", "C3", "T3"),
]


def fmt16(v):
    """Publication float formatting: <=16 decimals, trailing zeros stripped."""
    s = repr(float(v))
    if "e" in s or "E" in s or "." not in s:
        return s
    a, b = s.split(".")
    if len(b) <= 16:
        return s
    b = b[:16].rstrip("0")
    return a + "." + b if b else a


def fmt_cell(v, dtype):
    """Format one parsed pandas value the way the publication tables are formatted."""
    if pd.isna(v):
        return ""
    if not pd.api.types.is_numeric_dtype(dtype):
        return str(v)
    if pd.api.types.is_integer_dtype(dtype):
        return str(int(v))
    return fmt16(v)


def write_csv(path, header, rows):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def from_dataframe(df):
    return [[fmt_cell(df[k].iloc[i], df[k].dtype) for k in df.columns]
            for i in range(len(df))]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=OUT, help="output directory (default outputs/tables)")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    master = pd.read_csv(os.path.join(ANA, "master_results.csv"))
    summary = json.load(open(os.path.join(ANA, "summary.json")))

    # ---- Table 1: scenario design constants -----------------------------
    write_csv(os.path.join(args.out, "Table1_core_scenarios.csv"),
              ["Scenario", "tau", "gamma", "Informative study size", "Description"],
              [[s, fmt16(t), fmt16(g), inf, d] for s, t, g, inf, d in TABLE1_ROWS])

    # ---- Table 2: classification counts ---------------------------------
    rows2 = []
    for label, method, keep in TABLE2_SCOPES:
        sub = master[(master["method"] == method) & keep(master)]
        cls = sub["own_class"].value_counts()
        assert cls.get("BOUNDARY_ABOVE", 0) == 0
        rows2.append([label, len(sub),
                      int(cls.get("WITHIN", 0)), int(cls.get("OUTSIDE_ABOVE", 0)),
                      int(cls.get("BOUNDARY_BELOW", 0)), int(cls.get("OUTSIDE_BELOW", 0))])
    write_csv(os.path.join(args.out, "Table2_coverage_classification.csv"),
              ["Procedure/scope", "Cells", "Within", "Outside above",
               "Boundary below", "Outside below"], rows2)

    # ---- Table 3: estimand tracking at M=300 (core, Delta=0) ------------
    rows3 = []
    for sc, label, method, tgt in TABLE3_ROWS:
        cell = json.load(open(os.path.join(
            RES, "cell_summaries", f"core_{sc}_d1_M300.json")))
        md = cell["methods"][method]
        t, m = fmt16(cell["targets"][tgt]), fmt16(md["mean_est"])
        rows3.append([sc, label, tgt, t, m, fmt16(md["mcse_est"]),
                      repr(float(m) - float(t))])
    write_csv(os.path.join(args.out, "Table3_estimand_tracking_M300.csv"),
              ["Scenario", "Procedure", "Own target", "Target value",
               "Mean estimate", "MCSE(mean)", "Mean-target"], rows3)

    # ---- Table 4: secondary probes --------------------------------------
    n = summary["H4"]["mse_diff_sign_counts"]["n_cells"]
    worse = summary["H4"]["mse_diff_sign_counts"]["P4_worse_MSE"]
    worse_sig = len(summary["H4"]["mse_diff_sig"])
    prof = summary["P4_cov_diff_profile"]
    dfp = summary["df_probe_P2_vs_P1"]
    closer = dfp["closer_to_nominal"]
    dprobe = pd.read_csv(os.path.join(ANA, "df_and_covariance_probe.csv"))
    p3 = dprobe[(dprobe["cov_diff"] == "P3-P1") & (dprobe["delta"] == 0.5)]
    ratio = p3["width_probe"] / p3["width_P1"]
    i = ratio.idxmax()
    rows4 = [
        ["P4 correction", "MSE larger than P1", f"{worse}/{n}",
         f"{worse_sig}/{n} exceed +3 paired MCSE"],
        ["P4 correction", "Coverage numerically higher", f"{prof['positive_cells']}/{n}",
         f"{prof['significantly_positive']}/{n} significantly higher; "
         f"{prof['significantly_negative']}/{n} significantly lower"],
        ["Degrees of freedom", "Closer to 0.95",
         f"P1: {closer['P1']}; P2: {closer['P2']}; ties: {closer['tie']}",
         f"At M=3, P2 closer in {dfp['M3_cells_P2_closer']}/{dfp['M3_cells_total']}"],
        ["Covariance probe", "Maximum width inflation at \u0394\u22600",
         f"{(ratio.max() - 1.0) * 100:.1f}%", str(dprobe.loc[i, "cell"])],
    ]
    write_csv(os.path.join(args.out, "Table4_secondary_probes.csv"),
              ["Probe", "Metric", "Result", "Qualification"], rows4)

    # ---- Supplement tables S1-S7 ----------------------------------------
    s_tables = [
        ("TableS1_complete_master_results.csv", "master_results.csv", None),
        ("TableS2_P1_all_152_cells.csv", "master_results.csv",
         ["cell", "family", "scen", "delta", "M", "R", "T2", "T3", "T_DL",
          "mean_est", "mcse_est", "mean_width", "own_cov", "own_cov_mcse", "own_class"]),
        ("TableS3_df_covariance_probes.csv", "df_and_covariance_probe.csv", None),
        ("TableS4_P4_paired_results.csv", "H4_paired_P4_P1.csv", None),
        ("TableS5_variance_calibration.csv", "variance_calibration.csv", None),
        ("TableS6_estimand_map.csv", "estimand_map.csv", None),
        ("TableS7_sensitivity_summary.csv", "sensitivity_summary.csv", None),
    ]
    for out_name, src_name, cols in s_tables:
        df = pd.read_csv(os.path.join(ANA, src_name))
        if cols is not None:
            df = df.loc[df["method"] == "P1", cols]
        write_csv(os.path.join(args.out, out_name), list(df.columns),
                  from_dataframe(df))

    print("tables done:", sorted(os.listdir(args.out)))


if __name__ == "__main__":
    main()
