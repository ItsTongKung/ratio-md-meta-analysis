#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Article figures, generated exclusively from results/analysis/*.csv (no
hand-transcribed numbers). Generates the five publication figures (PNG +
PDF) into outputs/figures/, creating that directory when needed:

  Figure1_P1_core_coverage     published ratio procedure, coverage of T2,
                               60-cell core grid
  Figure2_scale_calibration    finite-sample interval-scale calibration
                               (core, Delta = 0)
  Figure3_estimand_tracking    distinct estimands under informative study size
  Figure4_P4_paired            prespecified plug-in bias correction relative
                               to the published ratio procedure
  FigureS1_sensitivity_S3      sensitivity behavior in informative-size
                               scenario S3

The figures use the article's public terminology and the article panel
geometry (fixed 220-dpi canvas). Generation is deterministic: identical
inputs produce byte-identical PNG rasters under the pinned environment
(verify_reproduction.py checks the inventory, dimensions, tolerance-band
rendering, and source data). Derived tables and figures are deliberately
not version-controlled; this script recreates them on demand.
"""
import csv, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ANA = os.path.join(HERE, "..", "results", "analysis")
FIG = os.path.join(HERE, "..", "outputs", "figures")

DPI = 220
BAND = "#e3eef5"          # coverage tolerance band 0.92-0.98
C0, C1, C2, C3, C4, C5 = [f"C{i}" for i in range(6)]


def read(name):
    with open(os.path.join(ANA, name)) as f:
        return list(csv.DictReader(f))


def fnum(x):
    return float(x) if x not in ("", None) else float("nan")


def save(fig, name):
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(os.path.join(FIG, name + ".png"), dpi=DPI)
    fig.savefig(os.path.join(FIG, name + ".pdf"))
    plt.close(fig)


def px_ax(fig, W, H, x0, x1, y0, y1):
    """Axes placed at pixel coordinates (y measured from the top)."""
    return fig.add_axes([x0 / W, (H - y1) / H, (x1 - x0) / W, (y1 - y0) / H])


def new_fig(w_px, h_px):
    return plt.figure(figsize=(w_px / DPI, h_px / DPI))


MS = [3, 5, 10, 30, 100, 300]
XL = (0.3456, 412.9)      # Figure 1 log-x window
XL2 = (2.38236, 377.626)  # Figure 2 log-x window (default margins)
XL3 = (2.36797, 377.86)   # Figure 3 log-x window
XL4 = (2.37936, 377.176)   # Figure 4 log-x window
XL5 = (4.5633, 32.77)     # Figure S1 log-x window


def fig1():
    rows = [r for r in read("H1_primary_coverage.csv")]
    W, H = 2611, 1239
    fig = new_fig(W, H)
    lefts = [130.5, 631.5, 1133.0, 1634.5, 2135.5]
    rights = [584.5, 1085.5, 1586.5, 2087.0, 2588.5]
    rowsy = [(221.5, 636.5), (723.5, 1138.5)]
    for i, scen in enumerate(["1", "2", "3", "4", "5"]):
        for j, d in enumerate(["0.0", "0.5"]):
            ax = px_ax(fig, W, H, lefts[i], rights[i], rowsy[j][0], rowsy[j][1])
            sub = {int(r["M"]): r for r in rows if r["scen"] == scen and r["delta"] == d}
            covs = [fnum(sub[M]["own_cov"]) for M in MS]
            err = [3 * fnum(sub[M]["own_cov_mcse"]) for M in MS]
            ax.axhspan(0.92, 0.98, color=BAND, zorder=0)
            ax.axhline(0.95, color=C0, lw=1.2, ls="--", zorder=1)
            ax.errorbar(MS, covs, yerr=err, marker="o", ms=5.8, capsize=2,
                        color=C0, zorder=3)
            ax.set_xscale("log")
            ax.set_xlim(*XL)
            ax.set_ylim(0.915079, 1.001989)
            ax.set_yticks([0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00])
            ax.set_title(f"S{scen}, \u0394={'0' if d == '0.0' else '0.5'}",
                         fontsize=9, pad=5.65)
            ax.set_xticks(MS)
            if j == 1:
                ax.set_xticklabels([str(m) for m in MS])
            else:
                ax.set_xticklabels([])
            ax.tick_params(labelsize=7)
            if j == 1:
                ax.set_xlabel("Number of studies (M)", fontsize=8.0)
            if i == 0:
                ax.set_ylabel("Coverage for $T_2$", fontsize=8, labelpad=3.0)
            else:
                ax.set_yticklabels([])
    fig.text(0.5, 1.0 - 22 / H, "Published ratio procedure: coverage of $T_2$ "
             "across the 60-cell core grid", ha="center", va="top", fontsize=11)
    save(fig, "Figure1_P1_core_coverage")


def fig2():
    rows = [r for r in read("variance_calibration.csv")
            if r["family"] == "core" and r["delta"] == "0.0"]
    W, H = 1736, 1070
    fig = new_fig(W, H)
    ax = px_ax(fig, W, H, 166.5, 1710.5, 70.5, 955.5)
    ax.axhline(1.0, color=C0, lw=1.2, ls="--", zorder=1)
    for scen in ["1", "2", "3", "4", "5"]:
        sub = {int(r["M"]): fnum(r["ratio_meanV_to_empSD"]) for r in rows
               if r["scen"] == scen}
        ax.plot(MS, [sub[M] for M in MS], marker="o", ms=5.8, label=f"S{scen}")
    ax.set_xscale("log")
    ax.set_xlim(*XL2)
    ax.set_ylim(0.712, 1.011)
    ax.set_yticks([0.75, 0.80, 0.85, 0.90, 0.95, 1.00])
    ax.set_xticks(MS)
    ax.set_xticklabels([str(m) for m in MS])
    ax.tick_params(labelsize=10)
    ax.set_xlabel("Number of studies (M)", fontsize=10)
    ax.set_ylabel("Mean reported scale / empirical SD of $\\hat{R}$", fontsize=10, labelpad=3.0)
    ax.legend(loc="lower right", ncol=5, fontsize=8, frameon=True)
    fig.text(938.5 / W, 1.0 - 24 / H, "Finite-sample interval-scale calibration of the "
             "published ratio procedure (core, \u0394 = 0)", ha="center",
             va="top", fontsize=12)
    save(fig, "Figure2_scale_calibration")


def fig3():
    rows = [r for r in read("estimand_map.csv")
            if r["scen"] in ("3", "5") and r["family"] == "core"]
    W, H = 1956, 1001
    fig = new_fig(W, H)
    panels = [("3", 151.5, 1021.5, 0.08361854460463167, 0.018306227057702698),
              ("5", 1064.5, 1933.5, 0.08361854460463167, 0.0073860233323441915)]
    for scen, x0, x1, t2, tdl in panels:
        ax = px_ax(fig, W, H, x0, x1, 238.5, 725.5)
        ax.axhline(t2, color=C0, lw=1.2, ls="-", zorder=1)
        ax.axhline(tdl, color=C0, lw=1.2, ls="--", zorder=1)
        ax.axhline(0.0, color=C0, lw=1.2, ls=":", zorder=1)
        for m, color in (("P1", C0), ("C1", C1), ("C3", C2)):
            sub = {int(r["M"]): fnum(r["mean_est"]) for r in rows
                   if r["scen"] == scen and r["method"] == m and r["delta"] == "0.0"}
            ax.plot(MS, [sub[M] for M in MS], marker="o", ms=5.8, color=color)
        ax.set_xscale("log")
        ax.set_xlim(*XL3)
        ax.set_ylim(-0.0061391, 0.0879665)
        ax.set_yticks([0.00, 0.02, 0.04, 0.06, 0.08])
        ax.set_xticks(MS)
        ax.set_xticklabels([str(m) for m in MS])
        ax.tick_params(labelsize=10)
        ax.set_xlabel("Number of studies (M)", fontsize=10)
        ax.set_title(f"Informative size: S{scen}, \u0394 = 0", fontsize=11.8, pad=6.4)
    fig.axes[0].set_ylabel("Mean point estimate", fontsize=10, labelpad=3.0)
    handles = [
        plt.Line2D([], [], color=C0, marker="o", lw=1.5,
                   label="Published ratio procedure"),
        plt.Line2D([], [], color=C1, marker="o", lw=1.5,
                   label="IVW-DL / Hartung\u2013Knapp center"),
        plt.Line2D([], [], color=C2, marker="o", lw=1.5, label="Equal weight"),
        plt.Line2D([], [], color=C0, lw=1.2, label="$T_2$"),
        plt.Line2D([], [], color=C0, lw=1.2, ls="--", label="$T_{DL}$"),
        plt.Line2D([], [], color=C0, lw=1.2, ls=":", label="$T_3$"),
    ]
    fig.legend(handles=handles, ncol=3, fontsize=7, frameon=True, columnspacing=2.1,
               loc="upper center", bbox_to_anchor=(0.5, 1.0 - 893 / H))
    fig.text(0.5, 1.0 - 22 / H, "Distinct estimands under informative study size",
             ha="center", va="top", fontsize=12)
    save(fig, "Figure3_estimand_tracking")


def fig4():
    rows = read("H4_paired_P4_P1.csv")
    fams = [("core", "Core", C0), ("D1", "D1 scale mixture", C1),
            ("D2", "D2 reduced skew", C2), ("D3", "D3 imbalanced allocation", C3),
            ("D4", "D4 non-normal outcomes", C4), ("clamp5", "Lower clamp", C5)]
    W, H = 1956, 998
    fig = new_fig(W, H)
    panels = [("mse_diff_P4_minus_P1", 151.5, 938.5, "MSE difference",
               "Bias-corrected - published ratio MSE",
               (-0.005025, 0.104289), [0.00, 0.02, 0.04, 0.06, 0.08, 0.10]),
              ("cov_diff_P4_minus_P1_vs_T2", 1147.5, 1933.5, "Coverage difference",
               "Bias-corrected - published ratio coverage",
               (-0.009509, 0.000924), [0.000, -0.002, -0.004, -0.006, -0.008])]
    for key, x0, x1, title, ylab, ylim, yticks in panels:
        ax = px_ax(fig, W, H, x0, x1, 232.5, 678.5)
        ax.axhline(0.0, color=C0, lw=1.2, ls="--", zorder=1)
        for fam, label, color in fams:
            sub = [(int(r["M"]), fnum(r[key])) for r in rows if r["family"] == fam]
            ax.plot([p[0] for p in sub], [p[1] for p in sub], "o", ms=5.5,
                    color=color, linestyle="none")
        ax.set_xscale("log")
        ax.set_xlim(*XL4)
        ax.set_ylim(*ylim)
        ax.set_yticks(yticks)
        ax.set_xticks(MS)
        ax.set_xticklabels([str(m) for m in MS])
        ax.tick_params(labelsize=10)
        ax.set_xlabel("Number of studies (M)", fontsize=10)
        ax.set_ylabel(ylab, fontsize=9.8, labelpad=3.0)
        ax.set_title(title, fontsize=12, pad=5.7)
    handles = [plt.Line2D([], [], color=c, marker="o", lw=0, label=lb)
               for _, lb, c in fams]
    fig.legend(handles=handles, ncol=3, fontsize=7, frameon=True, columnspacing=1.9,
               loc="upper center", bbox_to_anchor=(0.5, 1.0 - 899 / H))
    fig.text(0.5, 1.0 - 22 / H, "Prespecified plug-in bias correction relative "
             "to the published ratio procedure", ha="center", va="top",
             fontsize=12)
    save(fig, "Figure4_P4_paired")


def fig5():
    rows = [r for r in read("sensitivity_summary.csv") if r["method"] == "P1"]
    core = read("H1_primary_coverage.csv")
    fams = [("core", "Core", C0), ("D1", "D1 scale mixture", C1),
            ("D2", "D2 reduced skew", C2), ("D3", "D3 imbalanced allocation", C3),
            ("clamp5", "Lower clamp", C4)]
    W, H = 2171, 1098
    fig = new_fig(W, H)
    MS5 = [5, 10, 30]
    for j, d in enumerate(["0.0", "0.5"]):
        x0, x1 = (155.5, 1130.5) if j == 0 else (1174.5, 2148.5)
        ax = px_ax(fig, W, H, x0, x1, 252.5, 795.5)
        ax.axhspan(0.92, 0.98, color=BAND, zorder=0)
        ax.axhline(0.95, color=C0, lw=1.2, ls="--", zorder=1)
        for fam, label, color in fams:
            src = [r for r in (core if fam == "core" else rows)
                   if r["scen"] == "3" and r["delta"] == d
                   and r["family"] == fam]
            sub = {int(r["M"]): fnum(r["own_cov"]) for r in src if int(r["M"]) in MS5}
            ax.plot(MS5, [sub[M] for M in MS5], marker="o", ms=5.8, color=color,
                    label=label)
        if d == "0.5":
            clamp = [r for r in rows if r["family"] == "clamp5"
                     and r["delta"] == "0.5" and r["M"] == "10"]
            ax.annotate("Below 0.92 by >3 MCSE",
                        xy=(10, fnum(clamp[0]["own_cov"])),
                        xytext=(8.5, 0.9407), fontsize=10,
                        arrowprops=dict(arrowstyle="->", color="k", lw=1.0))
        ax.set_xscale("log")
        ax.set_xlim(*XL5)
        ax.set_ylim(0.914038, 0.982981)
        ax.set_yticks([0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98])
        ax.set_xticks([5, 10, 30])
        ax.set_xticklabels(["5", "10", "30"])
        ax.tick_params(labelsize=10)
        ax.set_xlabel("Number of studies (M)", fontsize=10)
        ax.set_title(f"S3, \u0394={'0' if d == '0.0' else '0.5'}", fontsize=12, pad=5.7)
        if j == 0:
            ax.set_ylabel("Coverage for $T_2$", fontsize=10, labelpad=3.0)
    handles = [plt.Line2D([], [], color=c, marker="o", lw=1.5, label=lb)
               for _, lb, c in fams]
    fig.legend(handles=handles, ncol=3, fontsize=8, frameon=True,
               loc="upper center", bbox_to_anchor=(0.4986, 1.0 - 984 / H))
    fig.text(0.5, 1.0 - 22 / H, "Sensitivity behavior in informative-size "
             "scenario S3", ha="center", va="top", fontsize=12)
    save(fig, "FigureS1_sensitivity_S3")


if __name__ == "__main__":
    fig1()
    fig2()
    fig3()
    fig4()
    fig5()
    print("figures done:", sorted(f for f in os.listdir(FIG)))
