# Finite-Sample Operating Characteristics and Estimand Characterization of a Ratio-Estimation Procedure for Mean-Difference Meta-Analysis

This is the reproducibility repository for the simulation study behind the
article of the same title.

**Purpose.** A published interval procedure for mean-difference
meta-analysis based on ratio estimation (the "P1" procedure) was evaluated
in a 152-cell simulation grid (16,203,932 unique replications) together
with three prespecified variants (P2–P4) and three standard comparator
procedures (random-effects inverse-variance weighting with
DerSimonian–Laird moment estimation and a z interval, conventional
Hartung–Knapp, and equal weighting). The study tracks three population
targets — the study-size-weighted superpopulation estimand, the
superpopulation mean of study effects, and the inverse-variance
probability limit — and describes where the published procedure is
calibrated, where it is not, and the interval-scale behavior that was
quantitatively consistent with the observed finite-sample coverage
pattern. This repository contains the simulation and analysis code,
configuration, deterministic seed specification, archived results, and
generation scripts needed to reproduce the article's numerical results,
tables, and publication figures.

Final article tables and figures are generated from the archived analysis
outputs and are not stored as version-controlled derived artifacts.

## Repository contents

```
ratio-md-meta-analysis/
├── README.md, CITATION.cff, requirements.txt,
│   LICENSE-CODE, LICENSE-CONTENT, .gitignore
├── code/                       seven Python scripts (see below)
├── config/
│   ├── simulation_design.yaml  full design as executed (scenarios, DGMs, boost)
│   └── seeds.md                master seed, SeedSequence keys, draw order
├── results/
│   ├── cell_summaries/         152 archived per-cell result JSON files
│   ├── analysis/               analysis CSVs + summary.json (archived)
│   ├── targets/                target files (exact, comparator, D4 Monte Carlo)
│   └── manifest_checksums.csv  SHA-256 of the archived results
└── docs/
    ├── METHODS.md              computational methods (estimators, targets, DGMs)
    └── REPRODUCIBILITY.md      how to reproduce (full and results-only)
```

`outputs/` is a generated working directory (gitignored, not part of the
repository contents): running the table/figure commands below creates it.

`code/` scripts:

| Script | Role |
|---|---|
| `methods.py` | data-generating models, estimators, design constants |
| `run_simulations.py` | simulation runner (152 cells; base grid + boost) |
| `compute_targets.py` | D4 comparator probability-limit target (Monte Carlo) |
| `analyze_results.py` | analysis pipeline (analysis CSVs + summary) |
| `make_tables.py` | article tables (11) → `outputs/tables/` |
| `make_figures.py` | article figures (5) → `outputs/figures/` |
| `verify_reproduction.py` | reader-facing verification of the archived results and of the generated tables/figures |

Archived in the repository: the simulation and analysis code,
configuration, seed documentation, the 152 cell summaries, the final
analysis datasets, the population targets, and the checksum manifest.
Generated on demand: the manuscript tables
(`python code/make_tables.py`) and manuscript figures
(`python code/make_figures.py`).

## Quick reproduction (results-only, no simulations)

Regenerates, from the archived cell summaries, the analysis outputs, all
11 tables, and all 5 article figures, then verifies everything:

```bash
python code/analyze_results.py
python code/make_tables.py
python code/make_figures.py
python code/verify_reproduction.py
```

Generation steps take a few seconds each; `verify_reproduction.py`
(additionally re-executing four cells from seed) takes ≈ 1 min. After
the generation steps, the article tables and figures appear under
`outputs/tables/` and `outputs/figures/`.

## Full reproduction (all simulations)

Re-runs the complete study from seed: 152 cells (base grid 7,840,000
replications plus the boundary-precision boost pass over the 32 cells at
the coverage tolerance band boundary), then the results-only steps above.

Run the D4 comparator-target Monte Carlo and the base grid:

```bash
python code/compute_targets.py
python code/run_simulations.py --family all
```

Then run the boundary-precision boost, passing the exact 32-cell list from
[docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md) as the `--boost`
argument (do not substitute a placeholder or shorten it), and finish with
the results-only commands above:

```bash
python code/analyze_results.py
python code/make_tables.py
python code/make_figures.py
python code/verify_reproduction.py
```

### Expected runtime (measured on the development machine)

D4 target Monte Carlo ≈ 45 s; base pass ≈ 105 s (7.84M replications);
boundary-precision boost pass ≈ 119 s; analysis ≈ 2 s, tables ≈ 3 s,
figures ≈ 4 s — full reproduction ≈ 5 min total (Python 3.12, single
process, NumPy vectorized). Runtimes will differ on other machines; the
results do not (see below).

## Reproducibility

Master seed 20260925; each cell uses
`numpy.random.SeedSequence(20260925, spawn_key=(family, scenario, delta_code, M))`
(draw order and stream conventions in [config/seeds.md](config/seeds.md)).
Under the pinned environment (below), regeneration is deterministic: the
archived cell summaries and analysis outputs and the generated tables
reproduce byte-for-bit, four re-executed cells match the archived
summaries exactly, and repeated figure renders reproduce the PNG rasters
byte-identically (PDF files carry identical vector content and differ
only in embedded creation timestamps). `verify_reproduction.py` reports
this.

## Article outputs (generated on demand)

| Article item | Generated by |
|---|---|
| Table 1 — core scenario design | `make_tables.py` → `outputs/tables/Table1_core_scenarios.csv` |
| Table 2 — coverage classification | `make_tables.py` → `outputs/tables/Table2_coverage_classification.csv` |
| Table 3 — large-M estimand tracking | `make_tables.py` → `outputs/tables/Table3_estimand_tracking_M300.csv` |
| Table 4 — secondary probes | `make_tables.py` → `outputs/tables/Table4_secondary_probes.csv` |
| Table S1 — complete master results | `make_tables.py` → `outputs/tables/TableS1_complete_master_results.csv` |
| Table S2 — all 152 cells, published procedure | `make_tables.py` → `outputs/tables/TableS2_P1_all_152_cells.csv` |
| Table S3 — df / covariance probes | `make_tables.py` → `outputs/tables/TableS3_df_covariance_probes.csv` |
| Table S4 — P4 paired results | `make_tables.py` → `outputs/tables/TableS4_P4_paired_results.csv` |
| Table S5 — variance calibration | `make_tables.py` → `outputs/tables/TableS5_variance_calibration.csv` |
| Table S6 — estimand map | `make_tables.py` → `outputs/tables/TableS6_estimand_map.csv` |
| Table S7 — sensitivity summary | `make_tables.py` → `outputs/tables/TableS7_sensitivity_summary.csv` |
| Figure 1 — core-grid coverage | `make_figures.py` → `outputs/figures/Figure1_P1_core_coverage.png` / `.pdf` |
| Figure 2 — scale calibration | `make_figures.py` → `outputs/figures/Figure2_scale_calibration.png` / `.pdf` |
| Figure 3 — estimand tracking | `make_figures.py` → `outputs/figures/Figure3_estimand_tracking.png` / `.pdf` |
| Figure 4 — P4 paired comparison | `make_figures.py` → `outputs/figures/Figure4_P4_paired.png` / `.pdf` |
| Figure S1 — sensitivity (S3) | `make_figures.py` → `outputs/figures/FigureS1_sensitivity_S3.png` / `.pdf` |

The machine-readable source data behind every table and figure is
`results/analysis/` (one row per cell × method) and
`results/analysis/summary.json`.

## Software environment

Python 3.12 (verified on 3.12.10) with numpy 2.4.6, scipy 1.17.1,
matplotlib 3.10.9, pandas 3.0.3 — pinned in
[requirements.txt](requirements.txt). Pure Python/NumPy; no compilation
is required:

```bash
python -m pip install -r requirements.txt
```

## License

Source code in `code/` is licensed under the MIT License (see
[LICENSE-CODE](LICENSE-CODE)). Documentation, configuration, and the
archived scientific results — including the tables and figures generated
from them — are licensed under the Creative Commons Attribution 4.0
International License (CC BY 4.0), unless otherwise noted (see
[LICENSE-CONTENT](LICENSE-CONTENT)).

## Citation

See [CITATION.cff](CITATION.cff). If you use this repository, please
cite the article and the repository version.

## DOI

Repository DOI: pending. A persistent identifier will be added after
archival release; this section and `CITATION.cff` will be updated then.
No placeholder DOI is used.

## Contact

Phongsakon Sriwicha —
[ORCID 0009-0001-0561-1076](https://orcid.org/0009-0001-0561-1076),
King Mongkut's University of Technology Thonburi, Bangkok, Thailand.
