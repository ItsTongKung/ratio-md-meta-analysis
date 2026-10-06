# Deterministic seeds and draw order

## Seed allocation

Every simulation cell draws from exactly one generator:

```python
rng = np.random.default_rng(
    np.random.SeedSequence(20260925,
                           spawn_key=(family_id, scenario_id, delta_id, M_id)))
```

- master seed: **20260925**
- `family_id`: core = 1, D1 = 2, D2 = 3, D3 = 4, D4 = 5, clamp5 = 6,
  targets = 7, verification = 8
- `scenario_id` ∈ {1, ..., 5} (S1–S5)
- `delta_id`: 1 ↔ Δ = 0, 2 ↔ Δ = 0.5
- `M_id` = number of studies per meta-analysis

The `(family_id, scenario_id, delta_id, M_id)` keys are unique across all 152
cells and across the target/verification families. The D4 comparator-target
computation (`code/compute_targets.py`) uses the targets family (id 7) with
spawn keys `(7, 4, 0, 0)` and `(7, 4, 1, 0)`, so the target pass is fully
separated from every simulation cell.

The underlying bit generator is PCG64 (the NumPy `default_rng` stream).
Re-executing a cell with the committed code and its seed key reproduces the
archived cell summary **bit-for-bit** (except the recorded wall time); this
was verified on the prespecified computational environment listed in
`docs/REPRODUCIBILITY.md`.

## Draw order within a cell (fixed)

Summary-statistic DGMs (core normal, D1 mixture, D2, D3, clamp-5), per chunk
of replications, per study:

1. `N = max(clamp, round(exp(mu0 + sigma0*Z)))` — one standard-normal draw
2. `q` — one draw per study (uniform choice over {1/3, 1/2, 2/3}; D3:
   uniform < 0.6 → 1/3, else 2/3); `nT = round(q*N)`, `nC = N - nT`
3. `s2p = chi2_{N-2}/(N-2)` — one gamma draw
4. effect noise `u` (S1: none; S2/S4: `tau * normal`; S3/S5: `tau * normal`)
5. within-study noise — one normal draw:
   `theta_hat = Theta + sqrt(c)*Z` (normal branch) or
   `theta_hat = Theta + sqrt(s2p*c)*Z` (D1 mixture branch),
   with `c = 1/nT + 1/nC`

D4 (patient-level) consumes N → q → effect noise → patient draws
(standardized-lognormal noise, control-arm block then treatment-arm block
within a chunk); `s2p` and the arm means are **realized** from the patient
data (no `s2p` draw).

## Chunking and stream conventions

Replications are generated in sequential chunks (10,000 replications per
chunk; 5,000 for D4). Every per-element draw consumes the underlying bit
stream in a fixed order, so results do not depend on the chunk size; the
archived runs used the chunk sizes above.

The boundary-precision boost pass (documented in `simulation_design.yaml`)
re-executes a cell at a larger total R **on the same generator**: the first
R_base replications are bit-identical to the base run and the boost extends
the stream.

## Reproducibility notes

- Aggregate statistics (means, variances, coverage counts, paired
  differences) are accumulated exactly across chunks — floating-point
  accumulation order is fixed by the chunk sequence, so chunked runs are
  bit-identical to single-pass runs at the same chunk boundaries.
- Bit-identity of re-execution was verified on Python 3.12.10 with
  numpy 2.4.6 / scipy 1.17.1 (see `docs/REPRODUCIBILITY.md`). Other NumPy
  versions can change distribution sampling; use the prespecified computational environment for
  exact bit-level reproduction, or expect results equal up to Monte Carlo
  stream differences otherwise.
- No other randomness enters the study: the analysis, table, and figure code
  are deterministic transformations of the archived cell summaries.
