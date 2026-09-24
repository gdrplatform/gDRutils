# drc parity probe

Prototype, not part of the package API. Branch material for GDR-3627.

`gDRutils` fits dose-response curves through `drc::drm()` with `drc::LL.4()` and
`drc::LL.3u()`. Everything around the solver — the F-test against the flat model, r2,
`xc50`, the AOC over a predicted grid, capping, priors, bounds and the constant-fit
fallback — is implemented here. drc is used as one thing: a bounded least-squares
solver for a 3- or 4-parameter log-logistic.

These scripts ask what happens if that one thing is replaced by
`stats::optim(method = "L-BFGS-B")` on the same objective, and what the difference looks
like where it matters.

## Scripts

| File | What it does |
|---|---|
| `drc_parity.R` | Fits every curve in the three public workshops datasets with both backends and compares fit type, `xc50`, parameters and residuals |
| `native_fit_fn.R` | Registers the native backend as an ordinary `fit_fn` through `apply_fit()`, writing into the same `Metrics` assay under its own `fit_source` |
| `results_summary.txt` | Console output of `drc_parity.R` from the run described below |

Both scripts take the same model, start values and bounds that `.prepareFitModel()`
builds, so the comparison isolates the solver.

## Running

`drc_parity.R` needs a [gDRworkshops](https://github.com/gdrplatform/gDRworkshops)
checkout for its input data:

```bash
git clone https://github.com/gdrplatform/gDRworkshops
GDRWORKSHOPS_EXAMPLES=gDRworkshops/examples Rscript drc_parity.R
```

`native_fit_fn.R` runs on the synthetic MAE shipped with `gDRtestData` and needs nothing
else.

## What the run showed

1229 curves — PRISM 774 (RV), Goetz 215 RV + 200 GR, Zhou 20 RV + 20 GR — through the
3-parameter fixed-`x_0` path.

- **The parameterisation decides the outcome.** Optimising `ec50` on the linear scale,
  the space drc declares it in, changes `fit_type` on 140 of 1229 curves. Optimising it
  on `log10`, nothing else changed, brings that to 13.
- **Those 13 are mostly the native fit doing better.** On 9 of them drc returns a
  constant fit and the native one a sigmoidal fit with higher r2 and lower F p-value
  (worst case r2 0.221 → 0.694). Two go the other way; two are curves where
  `drc::drm()` errors out.
- **Where both fit a curve** (696 of them), `xc50` agrees to a median relative
  difference of 1.3e-05, 86.9% within 0.1% and 94.1% within 1%.
- **Neither solver dominates on residuals**: native better on 343, tie on 106, drc
  better on 247.
- **Speed**: PRISM 774 curves in 1.7 s against 17.6 s. Not a like-for-like measurement —
  drc also builds a full model object.

Two optimisers on one objective, so this is a statement about agreement in the reported
metrics, not an independent validation of the fitting itself. Parameter-level tolerances
over all curves are not meaningful here: on flat curves the parameters are
unidentifiable and the two backends settle on opposite bounds with near-identical
residuals.

`native_fit_fn.R` carries its own caveat. Its demo fit function hardcodes the
relative-viability bound `x_inf >= 0` for every slice, where `fit_curves()` uses
`x_inf >= -1` for growth rate, and that single difference accounts for all 19 of its
disagreements on the synthetic data. The configuration around the solver matters as much
as the solver.
