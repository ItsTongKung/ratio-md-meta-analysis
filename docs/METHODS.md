# Computational methods

This document describes the computational content of the simulation study
"Finite-Sample Operating Characteristics and Estimand Characterization of a
Ratio-Estimation Procedure for Mean-Difference Meta-Analysis": the
data-generating mechanisms, the estimators, the population targets, and the
analysis quantities. It uses the same terminology as the article; the
statistical background and interpretation are in the article and its
supplement, not here.

## 1. Setting

A meta-analysis of \(M\) studies comparing treatment with control reports,
per study, a mean difference \(\hat\theta_i\) with variance
\(s_i^2 = s_{p,i}^2\,c_i\), where \(c_i = 1/n_{T,i} + 1/n_{C,i}\) and
\(s_{p,i}^2\) is the pooled variance. The response is the ratio estimator of
a mean difference on the patient scale,

\[
\widehat R \;=\; \frac{\sum_{i=1}^M N_i\hat\theta_i}{\sum_{i=1}^M N_i},
\]

with \(\{N_i\}\) the study sizes. The simulation study evaluates the
finite-sample operating characteristics of the published interval procedure
for \(\widehat R\) (procedure **P1** below) and of prespecified variants and
comparators, under a common data-generating mechanism (DGM) family and
several sensitivity mechanisms.

## 2. Data-generating mechanisms

All summary-statistic DGMs share, per study: a size draw, an allocation
draw, a pooled-variance draw, an effect draw, and a within-study noise draw
(exact order in `config/seeds.md`).

**Study size.** \(N = \max(\text{clamp}, \operatorname{round}(e^{\mu_0 +
\sigma_0 Z}))\), \(Z \sim N(0,1)\). Primary law: \(\mu_0 = 3\),
\(\sigma_0 = 0.75\), clamp = 10. Sensitivity laws: reduced study-size skew
(D2, \(\sigma_0 = 0.5\)) and the lower-clamp sensitivity with minimum study
size 5 (clamp-5).

**Allocation.** Treatment fraction \(q \in \{1/3, 1/2, 2/3\}\) with equal
probability (sensitivity D3: \(P(q = 1/3) = 0.6\), \(P(q = 2/3) = 0.4\)).
Then \(n_T = \operatorname{round}(qN)\), \(n_C = N - n_T\).

**Pooled variance.** \(s_p^2 = \chi^2_{N-2}/(N-2)\) (with \(N-2\) replaced
by \(\max(N-2, 1)\) in the draw).

**Effect and noise.** With \(\Delta\) the true mean difference and
scenarios S1–S5 defined by \((\tau, \gamma)\) (see
`config/simulation_design.yaml`), every scenario effect law uses one
standard-normal latent variable per study, \(Z_i \sim N(0,1)\),
independent of \((N_i, q_i)\) and of the pooled-variance and within-study
noise draws — in the implementation all study-level components are drawn
in sequence from one stream, so \(Z_i \perp (N_i, q_i)\) and the other
per-study draws (see `config/seeds.md`):

- S1 (common effect): \(\Theta_i = \Delta\)
- S2/S4 (heterogeneity \(\tau = 0.5/1.0\)): \(\Theta_i = \Delta + \tau Z_i\)
- S3/S5 (heterogeneity + informative size): \(\Theta_i = \Delta +
  \gamma\,(N_i - \mathrm{E}[N])/\mathrm{SD}(N) + \tau Z_i\) (actual-centered)

This is the repository's computational notation: the heterogeneity term
\(\tau Z_i\) is \(N(0, \tau^2)\), which the manuscript writes as
\(u_i \sim N(0, \tau^2)\). The standard-normal latent-variable form is used
here because it mirrors the implementation exactly.

Branches:

- **normal** (core, D2, D3, clamp-5): \(\hat\theta_i = \Theta_i +
  \sqrt{c_i}\,\varepsilon_i\), \(\varepsilon_i\) independent of \(s^2_{p,i}\)
  (corrected normal-outcome representation);
- **mixture** (D1, the variance-coupled normal scale-mixture sensitivity): \(\hat\theta_i = \Theta_i +
  \sqrt{s^2_{p,i} c_i}\,\varepsilon_i\) — a variance-coupled normal scale
  mixture;
- **patient-level** (D4, the non-normal patient-level-outcome sensitivity):
  patient-level outcomes are standardized lognormal,
  \(U = (e^Z - e^{1/2})/\sqrt{e(e-1)}\); arm means and \(s_p^2\) are
  **realized** from the patient data (control-arm block drawn first, then
  treatment-arm block); \(\hat\theta\) and \(s^2\) are their realized
  two-sample statistics.

## 3. Population targets

Three targets are tracked per cell (with \(\Delta\) added where the cell has
a nonzero \(\Delta\)):

- **T2** — the study-size-weighted superpopulation estimand,
  \(T_2 = \mathrm{E}[N\Theta]/\mathrm{E}[N]\), the probability limit of
  \(\widehat R\). Because \(T_3 = \mathrm{E}[\Theta]\),
  \(T_2 - T_3 = \operatorname{Cov}(N,\Theta)/\mathrm{E}[N]\) in general. It
  is zero for S1/S2/S4 (for S1 \(\Theta\) is a constant; for S2/S4
  \(\Theta\) is independent of \(N\)), and for the actual-centered
  informative-size scenarios (S3, S5) it reduces to
  \(T_2 - T_3 = \gamma\,\mathrm{SD}(N)/\mathrm{E}[N]\). With the
  prespecified \(\gamma = 0.1\) this equals \(\gamma\,\mathrm{CV}(N)\),
  exactly the archived law constants in
  `results/targets/exact_targets.json` (primary law: 0.08361854; D2:
  0.05215928; clamp-5: 0.08665149).
- **T3** — the superpopulation mean of study-level effects,
  \(\mathrm{E}[\Theta]\) (= \(\Delta\) in all scenarios here).
- **T_DL** — the probability limit of the inverse-variance (DerSimonian–Laird
  / Hartung–Knapp) center \(\mathrm{E}[w\hat\theta]/\mathrm{E}[w]\) with
  \(w = 1/(s^2 + \tau^2_*)\). Closed forms per family/scenario are in
  `results/targets/comparator_targets.json`; for D4 (realized arm statistics,
  no closed form) it was computed by Monte Carlo under the general C2 form at
the study level (`code/compute_targets.py` →
`results/targets/d4_targets.json`: \(\tau^2_* = 0.25901\),
\(T_{DL} - \Delta = 0.00015183\) (the committed run-1 value used
downstream), MCSE 0.00020).

## 4. Estimation procedures

All six procedures are evaluated on the same realized data in every
replication. Implementation: `code/methods.py` (`estimators`).

**P1 — published ratio procedure** (the published ratio procedure for the
mean difference, retained verbatim): With \(\bar y, \bar z\) the study-level
means of \(Y_i = N_i\hat\theta_i\) and \(N_i\), \(s_Y, s_Z\) their sample
SDs, and \(\hat\rho\) the sample correlation of \(Y_i\) and \(N_i\):

1. center: \(\widehat R = \sum N_i \hat\theta_i / \sum N_i\);
2. scale (delta-method, covariance term retained):
   \(\widehat V^2 = \frac{1}{M}\left(\frac{s_Y^2}{\bar z^2} +
   \frac{(\bar y s_Z)^2}{\bar z^4} - 2\hat\rho\,\frac{\bar y\, s_Y s_Z}{\bar z^3}\right)\);
3. interval: \(\widehat R \pm t_{M-2,0.975}\,\widehat V\).

**P2 — \(M-1\) degrees-of-freedom variant.** P1 with \(t_{M-1,0.975}\). Since
\(t_{M-2,.975} > t_{M-1,.975}\), the P2 interval is nested inside P1.

**P3 — no-covariance variant.** P1 with the covariance term dropped from the scale (the \(-2\hat\rho\)
term).

**P4 — plug-in bias-corrected variant.** Classical plug-in \(O(M^{-1})\)
ratio-bias correction of the P1 center (Cochran form), retaining P1's scale
and critical value:
\(\widehat K = \widehat R\, s_Z^2/\bar z^2 - \hat\rho\, s_Y s_Z/\bar z^2\),
center \(\widehat R - \widehat K/M\).

**C1 — IVW-DL (z) comparator**: DerSimonian–Laird moment estimator
\(\hat\tau^2 = \max(0, (Q - (M-1))/c^*)\), rerandomized weights
\(w_i = 1/(s_i^2 + \hat\tau^2)\), center \(\sum w_i \hat\theta_i / \sum
w_i\), interval \(\pm z_{0.975}\sqrt{1/\sum w_i}\).

**C2 — Hartung–Knapp comparator (conventional HKSJ)**: C1's center and rerandomized weights with the
Hartung–Knapp scale \(\widehat{se}^2 = \sum w_i(\hat\theta_i - \hat\mu)^2 /
((M-1)\sum w_i)\) and \(t_{M-1,0.975}\).

**C3 — equal-weight estimator**: unweighted mean \(\bar{\hat\theta}\) with
\(t_{M-1,0.975}\, s_{\hat\theta}/\sqrt M\).

Each interval is evaluated against the procedure's corresponding estimand:
P1–P4 against \(T_2\); C1/C2 against \(T_{DL}\); C3 against \(T_3\).

## 5. Simulation and analysis quantities

- Replications per cell: base 100,000 (core) / 20,000 (sensitivity),
  extended by the boundary-precision boost in the 32 cells at the coverage
  tolerance band boundary — the boost rule (trigger and plug-in form) was
  fixed before the boost pass, and per-cell boost sizes used the rule's
  stated plug-in form with the observed worst-boundary \(q(1-q)\) from the
  base run, a refinement that did not change the 0.92/0.98 classification
  boundary (see `config/simulation_design.yaml`); final unique total
  16,203,932.
- Per method: mean and SD of the point estimate (with MCSE), coverage of
  \(T_2, T_3, T_{DL}\) (each with MCSE), mean interval width, and rejection
  rates against 0 and \(T_3\).
- **Coverage tolerance band and classification** (operational Monte Carlo
  rule, not a hypothesis test): the prespecified coverage tolerance band is
  0.92–0.98 around nominal 0.95; it is neither a rejection region nor an
  externally validated adequacy criterion. Coverage for \(T_2\) inside the
  band is labeled WITHIN. Outside the band, the departure is labeled OUTSIDE
  if its distance to the violated boundary exceeds 3 MCSEs at the executed R,
  and BOUNDARY otherwise (the estimate lies outside 0.92/0.98 but within
  three MCSEs of the boundary). These labels are data values in the archived
  analysis files.
- **Paired comparisons** (plug-in variant vs published procedure;
  \(M-1\)/no-covariance variants vs P1): within-replication differences.
  For P2 (nested inside P1) the paired difference \(I_{P2} -
  I_{P1} \in \{0, -1\}\) has transition probability \(q = C_{P1} - C_{P2}\)
  and exact paired variance \(q(1-q)/R\); for P3 (not nested) the
  Fréchet bound is used.
- **P4 paired profile**: paired bias difference, MSE difference, and coverage
  difference, each with paired MCSE from within-replication differences.
- **P1 variance calibration**: \(\mathrm{E}[\widehat V]\),
  \(\mathrm{E}[\widehat V^2]\), the empirical SD and variance of
  \(\widehat R\), and mean \(\hat\rho\), reported as
  \(\mathrm{E}[\widehat V]/\mathrm{SD}(\widehat R)\) and
  \(\mathrm{E}[\widehat V^2]/\mathrm{Var}(\widehat R)\) ratios.
- **D4 comparator target**: two independent Monte Carlo passes of
  \(2\times 10^7\) studies each (`code/compute_targets.py`), batch-based
  MCSEs, with a stability check between passes.

## 6. Pipeline

```
compute_targets.py  ->  results/targets/d4_targets.json
run_simulations.py  ->  results/cell_summaries/<cell>.json   (152 files)
analyze_results.py  ->  results/analysis/*.csv + summary.json
make_tables.py      ->  outputs/tables/*.csv                 (11 tables)
make_figures.py     ->  outputs/figures/Figure*.png|pdf    (Figures 1-4, S1)
verify_reproduction.py  ->  reader-facing integrity checks
```

The table and figure outputs are generated on demand from the archived
analysis outputs (a generated working directory; not stored in the
repository).

## 7. Software

Python 3.12, NumPy 2.4.6, SciPy 1.17.1 (random streams: PCG64 via
`numpy.random.default_rng`; see `config/seeds.md`). Environment pins:
`requirements.txt`.
