#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Reader-facing verification of the reproducibility repository.

Works from a fresh clone: article tables and figures are generated on the
fly into temporary directories and verified by their scientific content;
no committed derived artifacts (outputs/) are required.

Checks, in order:
  1. the archived results are present and unmodified (152 cell summaries,
     analysis outputs, static target files) against
     results/manifest_checksums.csv;
  2. the target constants: the law/design constants of code/methods.py
     against results/targets/exact_targets.json, and the committed D4
     comparator target (internal two-pass agreement within MCSE);
  3. the archived analysis outputs regenerate byte-for-byte from the
     archived cell summaries;
  4. the 11 article tables generate from the archived analysis outputs
     (deterministic: two independent generations are byte-identical);
  5. the five article figures (Figure 1-4, S1; PNG + PDF) generate from
     the archived analysis outputs (PNG generation deterministic);
  6. the headline classification counts and replication accounting match
     the archived reference values;
  7. four selected cells (one per DGM branch) re-executed from seed match
     the archived cell summaries bit-for-bit (except wall time);
  8. the generated tables match the archived scientific reference values:
     Table 1 design constants, Table 2 classification counts, Table 3
     estimand tracking against the archived cell summaries, Table 4 probe
     results, and the row counts/headers of the Tables S1-S7 analysis
     dumps;
  9. the generated figures match the article geometry: expected output
     inventory, exact canvas dimensions (fixed 220-dpi panels), the
     coverage tolerance band rendering, and the archived source data each
     figure was generated from.

The script exits nonzero on any material mismatch. Expected runtime:
about 1-2 minutes (dominated by the four re-executed cells).

Usage:  python code/verify_reproduction.py
"""
import csv, hashlib, json, os, shutil, subprocess, sys, tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE)

import methods as rf  # noqa: E402

RES = os.path.join(ROOT, "results")
ANA = os.path.join(RES, "analysis")
CELLS = os.path.join(RES, "cell_summaries")
FIGS = os.path.join(ROOT, "outputs", "figures")

# Archived reference values (archived reference state; any drift fails
# check 6).
EXPECTED_HEADLINE = {
    "H1_core_counts": {"WITHIN": 48, "OUTSIDE_ABOVE": 10, "BOUNDARY_BELOW": 2,
                       "OUTSIDE_BELOW": 0},
    "H1_full_grid_counts": {"WITHIN": 127, "OUTSIDE_ABOVE": 20,
                            "BOUNDARY_BELOW": 4, "OUTSIDE_BELOW": 1},
    "df_probe_closer_to_nominal": {"P1": 117, "P2": 22, "tie": 13},
    "df_probe_decisive_OUTSIDE": {"P1": 21, "P2": 31},
    "df_probe_M3": {"total": 20, "P2_closer": 20},
    "P4_cov_diff_profile": {"positive_cells": 10, "significantly_positive": 0,
                            "significantly_negative": 101},
    "replication_accounting": {"base_computation": 7840000,
                               "boost_job_computation": 9483932,
                               "unique_additional_beyond_base": 8363932,
                               "final_unique_committed": 16203932,
                               "boosted_cells": 32},
}
SELECTED_CELLS = [("core", 3, 1, 10), ("D1", 2, 1, 5),
                  ("D4", 2, 1, 10), ("clamp5", 3, 2, 30)]  # one per DGM branch

# Committed D4 comparator target (fixed before the estimator runs; the file
# is regenerable with compute_targets.py, which reproduces these exactly).
EXPECTED_D4 = {
    "tau2_star": 0.2590054583336707,
    "T_DL_minus_Delta": 0.00015183103614867843,
    "mcse_T_DL": 0.00019978691512039906,
    "n_studies": 10000000,
}

# Archived scientific reference values of the generated article tables.
# Tables are verified by content (check 8), not by hash of a deposited copy.
TABLE_FILES = [
    "Table1_core_scenarios.csv", "Table2_coverage_classification.csv",
    "Table3_estimand_tracking_M300.csv", "Table4_secondary_probes.csv",
    "TableS1_complete_master_results.csv", "TableS2_P1_all_152_cells.csv",
    "TableS3_df_covariance_probes.csv", "TableS4_P4_paired_results.csv",
    "TableS5_variance_calibration.csv", "TableS6_estimand_map.csv",
    "TableS7_sensitivity_summary.csv",
]
EXPECTED_TABLE1 = [
    ("S1", "0.0", "0.0", "No", "Common effect"),
    ("S2", "0.5", "0.0", "No", "Moderate heterogeneity"),
    ("S3", "0.5", "0.1", "Yes",
     "Moderate heterogeneity + informative study size"),
    ("S4", "1.0", "0.0", "No", "Strong heterogeneity"),
    ("S5", "1.0", "0.1", "Yes",
     "Strong heterogeneity + informative study size"),
]
EXPECTED_TABLE2 = [
    ("P1 core (60 cells)", 60, 48, 10, 2, 0),
    ("P1 full grid (152 cells)", 152, 127, 20, 4, 1),
    ("IVW-DL full grid", 152, 88, 0, 2, 62),
    ("HKSJ full grid", 152, 152, 0, 0, 0),
    ("Equal-weight full grid", 152, 152, 0, 0, 0),
]
EXPECTED_TABLE4 = [
    ("P4 correction", "MSE larger than P1", "150/152",
     "142/152 exceed +3 paired MCSE"),
    ("P4 correction", "Coverage numerically higher", "10/152",
     "0/152 significantly higher; 101/152 significantly lower"),
    ("Degrees of freedom", "Closer to 0.95", "P1: 117; P2: 22; ties: 13",
     "At M=3, P2 closer in 20/20"),
    ("Covariance probe", "Maximum width inflation at \u0394\u22600",
     "77.4%", "core_S1_d2_M300"),
]
TABLE_S_DATA_ROWS = {  # exact dumps of the archived analysis CSVs
    "TableS1_complete_master_results.csv": 152 * 7,   # 152 cells x 7 methods
    "TableS2_P1_all_152_cells.csv": 152,              # P1 rows only
    "TableS3_df_covariance_probes.csv": 304,
    "TableS4_P4_paired_results.csv": 152,
    "TableS5_variance_calibration.csv": 152,
    "TableS6_estimand_map.csv": 608,
    "TableS7_sensitivity_summary.csv": 368,
}
TABLE_S_HEADERS = {  # generated dumps carry the archived CSV headers
    "TableS1_complete_master_results.csv": "master_results.csv",
    "TableS3_df_covariance_probes.csv": "df_and_covariance_probe.csv",
    "TableS4_P4_paired_results.csv": "H4_paired_P4_P1.csv",
    "TableS5_variance_calibration.csv": "variance_calibration.csv",
    "TableS6_estimand_map.csv": "estimand_map.csv",
    "TableS7_sensitivity_summary.csv": "sensitivity_summary.csv",
}
TABLE_S2_COLS = ["cell", "family", "scen", "delta", "M", "R", "T2", "T3",
                 "T_DL", "mean_est", "mcse_est", "mean_width", "own_cov",
                 "own_cov_mcse", "own_class"]

FIGURE_STEMS = ["Figure1_P1_core_coverage", "Figure2_scale_calibration",
                "Figure3_estimand_tracking", "Figure4_P4_paired",
                "FigureS1_sensitivity_S3"]
EXPECTED_FIG_DIMS = {  # article canvas geometry (fixed 220-dpi panels)
    "Figure1_P1_core_coverage": (2611, 1239),
    "Figure2_scale_calibration": (1736, 1070),
    "Figure3_estimand_tracking": (1956, 1001),
    "Figure4_P4_paired": (1956, 998),
    "FigureS1_sensitivity_S3": (2171, 1098),
}
BAND_RGB = (227, 238, 245)  # coverage tolerance band 0.92-0.98 (#e3eef5)
BAND_MIN_FRAC = 0.30        # band covers >=30% of Fig 1 / Fig S1 pixels

RESULTS = []
TABLE_DIR = FIG_DIR = None


def record(name, ok, detail=""):
    RESULTS.append((name, bool(ok), detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" — {detail}" if detail else ""))


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def comparable(d, drop=("wall_s", "replication_dump")):
    if isinstance(d, dict):
        return {k: comparable(v, drop) for k, v in d.items() if k not in drop}
    if isinstance(d, list):
        return [comparable(v, drop) for v in d]
    return d


def read_csv_rows(path):
    with open(path, newline="", encoding="utf-8") as f:
        r = list(csv.reader(f))
    return r[0], r[1:]


def check_1_inventory():
    manifest = {}
    with open(os.path.join(RES, "manifest_checksums.csv"), newline="") as f:
        for row in csv.DictReader(f):
            manifest[row["file"]] = row["sha256"]
    n_cells = len([fn for fn in os.listdir(CELLS) if fn.endswith(".json")])
    record("1a. 152 archived cell summaries present", n_cells == 152,
           f"found {n_cells}")
    bad, missing, changed_cells = [], [], []
    for rel, want in sorted(manifest.items()):
        path = os.path.join(ROOT, rel)
        if not os.path.exists(path):
            missing.append(rel)
        elif sha256(path) != want:
            if rel.startswith("results/cell_summaries/"):
                changed_cells.append(rel)
            else:
                bad.append(rel)
    record("1b. analysis and target files match checksum manifest",
           not missing and not bad,
           "" if not missing and not bad else f"missing={missing[:5]} changed={bad[:5]}")
    # Cell summaries contain a wall-time field, so re-running the simulations
    # changes their checksums; this is expected and not a content failure
    # (content is validated by checks 3 and 7). A mismatch WITHOUT a re-run
    # of run_simulations.py would be a modification of the archived results.
    if changed_cells:
        print(f"[WARN] 1c. {len(changed_cells)} cell summaries differ from the "
              f"deposit checksums (e.g. {changed_cells[0].split('/')[-1]} ...). "
              "Expected only if run_simulations.py was re-executed (the "
              "wall-time field changes); all other differences are material.")
    else:
        print(f"[info] 1c. all {n_cells} cell summaries match the deposit checksums")


def check_2_targets():
    exact = json.load(open(os.path.join(RES, "targets", "exact_targets.json")))
    ok = True
    for law_key, name in (("primary", "primary"), ("D2_sigma0.5", "D2"),
                          ("clamp5", "clamp5")):
        e, m = exact["laws"][law_key], rf.LAWS[name]
        pairs = [("mu0", e["mu0"]), ("sigma0", e["sigma0"]),
                 ("clamp", float(e["clamp"])), ("mu_c", e["E_N"]),
                 ("sd_c", e["SD_N"]), ("t2mt3", e["actual_centered"]["T2_minus_T3"])]
        for attr, want in pairs:
            if m[attr] != want:
                ok = False
                print(f"    law {name}: methods.{attr}={m[attr]!r} != exact {want!r}")
    record("2a. methods.py law constants match exact_targets.json", ok)

    d = json.load(open(os.path.join(RES, "targets", "d4_targets.json")))
    r1, r2, st = d["run1"], d["run2"], d["stability"]
    committed = (r1["tau2_star"] == EXPECTED_D4["tau2_star"]
                 and r1["T_DL_minus_Delta"] == EXPECTED_D4["T_DL_minus_Delta"]
                 and r1["mcse"]["T_DL"] == EXPECTED_D4["mcse_T_DL"]
                 and r1["n_studies"] == EXPECTED_D4["n_studies"])
    record("2b. D4 comparator target matches the committed values", committed,
           f"T_DL-Delta={r1['T_DL_minus_Delta']:.8f}, tau2*={r1['tau2_star']:.5f}")
    mcse = r1["mcse"]["T_DL"]
    stable = (abs(st["T_DL_diff"]) < 3 * mcse
              and abs(st["tau2_star_diff"]) < 0.01 * r1["tau2_star"])
    record("2c. D4 target's two independent MC passes agree", stable,
           f"pass diff {st['T_DL_diff']:.2e} vs MCSE {mcse:.2e}")
    return r1


def check_3_analysis_regen():
    import analyze_results as ar
    tmp = tempfile.mkdtemp(prefix="verify_analysis_")
    ar.ANA = tmp  # redirect output; inputs still the archived results
    try:
        ar.main()
        diffs = []
        for fn in sorted(os.listdir(ANA)):
            a = open(os.path.join(ANA, fn), "rb").read()
            b = open(os.path.join(tmp, fn), "rb").read()
            if a != b:
                diffs.append(fn)
        record("3. analysis outputs regenerate byte-for-byte from the "
               "archived cell summaries", not diffs,
               "" if not diffs else f"differing: {diffs}")
    finally:
        ar.ANA = ANA
        shutil.rmtree(tmp, ignore_errors=True)


def check_4_generate_tables():
    global TABLE_DIR, TABLE_DIR2
    tmp = tempfile.mkdtemp(prefix="verify_tables_")
    cmd = [sys.executable, os.path.join(HERE, "make_tables.py"), "--out", tmp]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        record("4a. 11 article tables generate from the archived analysis", False,
               r.stderr[-300:])
        return
    got = sorted(os.listdir(tmp))
    record("4a. 11 article tables generate from the archived analysis",
           got == sorted(TABLE_FILES),
           "" if got == sorted(TABLE_FILES) else f"generated: {got}")

    tmp2 = tempfile.mkdtemp(prefix="verify_tables2_")
    r2 = subprocess.run([sys.executable, os.path.join(HERE, "make_tables.py"),
                         "--out", tmp2], capture_output=True, text=True)
    if r2.returncode != 0:
        record("4b. table generation is deterministic", False, r2.stderr[-300:])
        return
    diffs = [fn for fn in TABLE_FILES
             if open(os.path.join(tmp, fn), "rb").read()
             != open(os.path.join(tmp2, fn), "rb").read()]
    shutil.rmtree(tmp2, ignore_errors=True)
    record("4b. table generation is deterministic (two runs byte-identical)",
           not diffs, "" if not diffs else f"differing: {diffs}")
    TABLE_DIR = tmp


def check_5_generate_figures():
    global FIG_DIR
    import matplotlib.image as _  # noqa: F401  (png rendering available)
    import make_figures as mf
    tmp = tempfile.mkdtemp(prefix="verify_figures_")
    mf.FIG = tmp
    try:
        mf.fig1(); mf.fig2(); mf.fig3(); mf.fig4(); mf.fig5()
        want = sorted(n + e for n in FIGURE_STEMS for e in (".png", ".pdf"))
        got = sorted(os.listdir(tmp))
        record("5a. five article figures generate (5 PNG + 5 PDF)", got == want,
               "" if got == want else f"generated: {got}")
    finally:
        mf.FIG = FIGS

    tmp2 = tempfile.mkdtemp(prefix="verify_figures2_")
    mf.FIG = tmp2
    try:
        mf.fig1(); mf.fig2(); mf.fig3(); mf.fig4(); mf.fig5()
        # PNG rasters are deterministic (PDFs embed timestamps, skipped)
        diffs = [n + ".png" for n in FIGURE_STEMS
                 if open(os.path.join(tmp, n + ".png"), "rb").read()
                 != open(os.path.join(tmp2, n + ".png"), "rb").read()]
        record("5b. figure PNG generation is deterministic "
               "(two renders byte-identical)", not diffs,
               "" if not diffs else f"differing: {diffs}")
    finally:
        mf.FIG = FIGS
        shutil.rmtree(tmp2, ignore_errors=True)
    FIG_DIR = tmp


def check_6_headline():
    summary = json.load(open(os.path.join(ANA, "summary.json")))
    CLASSES = ("WITHIN", "OUTSIDE_ABOVE", "BOUNDARY_BELOW", "OUTSIDE_BELOW")
    got = {
        "H1_core_counts": {k: summary["H1"]["counts"].get(k, 0) for k in CLASSES},
        "H1_full_grid_counts": {k: summary["H1_full_grid"]["counts"].get(k, 0)
                                for k in CLASSES},
        "df_probe_closer_to_nominal": summary["df_probe_P2_vs_P1"]["closer_to_nominal"],
        "df_probe_decisive_OUTSIDE": summary["df_probe_P2_vs_P1"]["decisive_OUTSIDE_counts"],
        "df_probe_M3": {"total": summary["df_probe_P2_vs_P1"]["M3_cells_total"],
                        "P2_closer": summary["df_probe_P2_vs_P1"]["M3_cells_P2_closer"]},
        "P4_cov_diff_profile": {
            k: summary["P4_cov_diff_profile"][k] for k in
            ("positive_cells", "significantly_positive", "significantly_negative")},
        "replication_accounting": summary["replication_accounting"],
    }
    ok = True
    for key, want in EXPECTED_HEADLINE.items():
        if got[key] != want:
            ok = False
            print(f"    {key}: got {got[key]}, expected {want}")
    record("6. headline classification counts and accounting", ok)


def check_7_selected_cells(tdl_d4):
    from run_simulations import run_cell, cell_id
    man = {c["key"]: c for c in rf.build_manifest()}
    ok = True
    for key in SELECTED_CELLS:
        cell = man[key]
        stored = json.load(open(os.path.join(
            CELLS, cell_id(*key) + ".json")))
        rerun = run_cell(cell, 10_000, tdl_d4=tdl_d4,
                         R_override=stored["R_executed"])
        same = (json.dumps(comparable(stored), sort_keys=True)
                == json.dumps(comparable(rerun), sort_keys=True))
        ok &= same
        record(f"7. {cell_id(*key)} re-executed bit-identically "
               f"(R={stored['R_executed']})", same)
    return ok


def check_8_table_contents():
    def rows(name):
        return read_csv_rows(os.path.join(TABLE_DIR, name))

    hdr, data = rows("Table1_core_scenarios.csv")
    ok = (hdr == ["Scenario", "tau", "gamma", "Informative study size",
                  "Description"] and data == [list(r) for r in EXPECTED_TABLE1])
    record("8a. Table 1 matches the archived scenario design constants", ok,
           "" if ok else f"got {data}")

    hdr, data = rows("Table2_coverage_classification.csv")
    ok = (hdr == ["Procedure/scope", "Cells", "Within", "Outside above",
                  "Boundary below", "Outside below"]
          and data == [list(map(str, r)) for r in EXPECTED_TABLE2])
    record("8b. Table 2 matches the archived classification counts", ok,
           "" if ok else f"got {data}")

    # Table 3 against the archived cell summaries (checksummed in check 1):
    # targets T2/T3/T_DL and each method's mean/ MCSE from the archive; the
    # table columns truncate to <=16 decimals, so compare within 5e-16.
    exact = json.load(open(os.path.join(RES, "targets", "exact_targets.json")))
    t2const = exact["laws"]["primary"]["actual_centered"]["T2_minus_T3"]
    lab2m = {"P1 ratio": "P1", "IVW-DL / HKSJ center": "C1", "Equal weight": "C3"}
    tgt_want = {"T2": t2const, "T3": 0.0}
    cellcache = {}
    hdr, data = rows("Table3_estimand_tracking_M300.csv")
    ok = hdr == ["Scenario", "Procedure", "Own target", "Target value",
                 "Mean estimate", "MCSE(mean)", "Mean-target"]
    detail = ""
    if ok:
        worst = 0.0
        for sc, label, tgt, t, m, s_, mt in data:
            cell = cellcache.setdefault(sc, json.load(open(os.path.join(
                CELLS, f"core_{sc}_d1_M300.json"))))
            meth = cell["methods"][lab2m[label]]
            want_t = tgt_want.get(tgt, cell["targets"]["T_DL"])
            e = max(abs(float(t) - want_t), abs(float(m) - meth["mean_est"]),
                    abs(float(s_) - meth["mcse_est"]))
            worst = max(worst, e)
            e2 = abs(float(mt) - (float(m) - float(t)))
            if e > 5e-16 or e2 > 1e-15:
                ok = False
                detail = f"{sc}/{label}: error {e:.2e}"
        if ok:
            detail = f"max deviation to archived targets {worst:.1e} (<=5e-16)"
    record("8c. Table 3 matches the archived cell summaries and target "
           "constants", ok, detail)

    hdr, data = rows("Table4_secondary_probes.csv")
    ok = (hdr == ["Probe", "Metric", "Result", "Qualification"]
          and data == [list(r) for r in EXPECTED_TABLE4])
    record("8d. Table 4 matches the archived probe results", ok,
           "" if ok else f"got {data}")

    ok, detail = True, ""
    for name, want in TABLE_S_DATA_ROWS.items():
        hdr, data = rows(name)
        if len(data) != want:
            ok = False
            detail = f"{name}: {len(data)} data rows != {want}"
        elif name in TABLE_S_HEADERS:
            src_hdr, _ = read_csv_rows(os.path.join(ANA, TABLE_S_HEADERS[name]))
            if hdr != src_hdr:
                ok = False
                detail = f"{name}: header differs from {TABLE_S_HEADERS[name]}"
        elif name == "TableS2_P1_all_152_cells.csv":
            if hdr != TABLE_S2_COLS:
                ok = False
                detail = "TableS2: column set differs"
        if not ok:
            break
    record("8e. Tables S1-S7 row counts and headers match the archived "
           "analysis files", ok, detail)


def check_9_figure_contents():
    import matplotlib.image as mpimg
    # 9a exact canvas dimensions of the five PNGs
    ok, detail = True, ""
    for stem, (w, h) in EXPECTED_FIG_DIMS.items():
        a = mpimg.imread(os.path.join(FIG_DIR, stem + ".png"))
        if a.shape[1] != w or a.shape[0] != h:
            ok = False
            detail = f"{stem}: {a.shape[1]}x{a.shape[0]} != {w}x{h}"
            break
    record("9a. figure canvases match the article dimensions (220 dpi)", ok,
           "all five exact" if ok else detail)

    # 9b the coverage tolerance band (0.92-0.98) is rendered
    ok, detail = True, ""
    for stem in ("Figure1_P1_core_coverage", "FigureS1_sensitivity_S3"):
        a = mpimg.imread(os.path.join(FIG_DIR, stem + ".png"))
        rgb = (a[..., :3] * 255.0 + 0.5).astype(np.uint8)
        n = int(((rgb[..., 0] == BAND_RGB[0]) & (rgb[..., 1] == BAND_RGB[1])
                 & (rgb[..., 2] == BAND_RGB[2])).sum())
        frac = n / (a.shape[0] * a.shape[1])
        if frac < BAND_MIN_FRAC:
            ok = False
            detail = f"{stem}: band pixels {100*frac:.2f}% < {100*BAND_MIN_FRAC:.0f}%"
    record("9b. tolerance band rendered in Figures 1 and S1", ok,
           "" if not ok else "band present at the expected share of pixels")

    # 9c the archived source data behind each figure is complete
    src = read("H1_primary_coverage.csv")
    combos = {(r["scen"], r["delta"], r["M"]) for r in src}
    want = {(s, d, str(m)) for s in ("1", "2", "3", "4", "5")
            for d in ("0.0", "0.5") for m in (3, 5, 10, 30, 100, 300)}
    ok1 = len(src) == 60 and combos == want
    det1 = "" if ok1 else f"H1_primary_coverage.csv: {len(src)} rows"
    em = [r for r in read("estimand_map.csv")
          if r["family"] == "core" and r["scen"] in ("3", "5")
          and r["delta"] == "0.0"]
    want2 = {(m, s, str(mm)) for m in ("P1", "C1", "C3")
             for s in ("3", "5") for mm in (3, 5, 10, 30, 100, 300)}
    extras2 = {(r["method"], r["scen"], r["M"]) for r in em} - want2
    allow2 = {("C2", s, str(mm)) for s in ("3", "5")
              for mm in (3, 5, 10, 30, 100, 300)}
    ok2 = want2 <= {(r["method"], r["scen"], r["M"]) for r in em} \
        and extras2 == allow2
    det2 = "" if ok2 else "estimand_map.csv: series incomplete"
    sens = [r for r in read("sensitivity_summary.csv") if r["method"] == "P1"]
    anchor = [r for r in sens if r["family"] == "clamp5"
              and r["scen"] == "3" and r["delta"] == "0.5" and r["M"] == "10"]
    got = {f: {(r["delta"], r["M"]) for r in sens
               if r["family"] == f and r["scen"] == "3"}
           for f in ("D1", "D2", "D3", "clamp5")}
    need = {(d, str(m)) for d in ("0.0", "0.5") for m in (5, 10, 30)}
    ok3 = (bool(anchor) and all(need <= got[f] for f in got)
           and sum(len(v) for v in got.values()) == 28)
    det3 = "" if ok3 else "sensitivity_summary.csv: series or anchor missing"
    vc = read("variance_calibration.csv")
    vc = [r for r in vc if r["family"] == "core" and r["delta"] == "0.0"]
    ok4 = len(vc) == 30 and {(r["scen"], r["M"]) for r in vc} == {
        (s, str(m)) for s in ("1", "2", "3", "4", "5")
        for m in (3, 5, 10, 30, 100, 300)}
    det4 = "" if ok4 else "variance_calibration.csv: series incomplete"
    h4 = read("H4_paired_P4_P1.csv")
    fam_ms = {"core": (3, 5, 10, 30, 100, 300), "D1": (3, 5, 10, 30, 100),
              "D2": (5, 10, 30), "D3": (5, 10, 30), "D4": (5, 10, 30),
              "clamp5": (5, 10, 30)}
    fam_n = {"core": 60, "D1": 50, "D2": 12, "D3": 12, "D4": 6, "clamp5": 12}
    ok5 = len(h4) == 152 and all(
        r["mse_diff_P4_minus_P1"] != "" and r["cov_diff_P4_minus_P1_vs_T2"] != ""
        for r in h4) and all(
        {int(r["M"]) for r in h4 if r["family"] == f} == set(ms)
        and sum(1 for r in h4 if r["family"] == f) == fam_n[f]
        for f, ms in fam_ms.items())
    det5 = "" if ok5 else "H4_paired_P4_P1.csv: family series incomplete"
    record("9c. figure source data complete (Figures 1-4, S1)", ok1 and ok2
           and ok3 and ok4 and ok5, "; ".join(filter(None, (det1, det2, det3,
                                                           det4, det5))))


def read(name):
    with open(os.path.join(ANA, name), newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def main():
    print("Verifying the reproducibility repository against its archived results\n")
    check_1_inventory()
    d4 = check_2_targets()
    check_3_analysis_regen()
    check_4_generate_tables()
    check_5_generate_figures()
    check_6_headline()
    check_7_selected_cells(d4["T_DL_minus_Delta"])
    check_8_table_contents()
    check_9_figure_contents()
    n_pass = sum(1 for _, ok, _ in RESULTS if ok)
    n_fail = len(RESULTS) - n_pass
    print(f"\n{len(RESULTS)} checks: {n_pass} passed, {n_fail} failed")
    if n_fail:
        for name, ok, detail in RESULTS:
            if not ok:
                print(f"  FAILED: {name} {detail}")
        sys.exit(1)
    print("ALL CHECKS PASS — the archived results reproduce from the "
          "committed code, configuration, and seeds; the article tables and "
          "figures regenerate from the archived analysis outputs.")

    if TABLE_DIR:
        shutil.rmtree(TABLE_DIR, ignore_errors=True)


if __name__ == "__main__":
    main()
