# Method provenance and implementation decisions

## Local source material inspected

- `Rcode_Revise/MAL_function.R`: binary median adaptive selection,
  gamma/lambda grid, wAMD tuning, logistic refit and survey-weighted effect.
- `Rcode_Revise/MAL_function_gamma5.R` and `MAL_function_gamma10.R`:
  sensitivity versions using convergence constants 5 and 10.
- `GOAL-main/R/GOAL.R`: continuous Gaussian treatment selection, npCBPS,
  weighted correlation tuning and known-covariate logic. This older method
  uses a mean outcome regression; GMAL uses a median outcome regression.
- Final local manuscript: median regression, continuous DWC and binary wAMD;
  real-data application uses bootstrap intervals and requires prescreening.

The old GOAL package is authored by Qian Gao and is not relabeled or copied
wholesale. Functions here were rewritten around an explicit data interface.
The manuscript's method and all paper authors are credited in CITATION.

## Explicit changes from the supplied scripts

1. `gmal()` reads its own data arguments; it never looks up global Y or tau,
   sources absolute paths, clears the workspace, or attaches libraries.
2. The outcome quantile is fixed to 0.5 (the published median method).
3. Binary and continuous workflows are exposed explicitly, not inferred.
4. Standardization uses sample standard deviations. Known covariates are
   standardized too and receive zero penalty; this extends the incomplete
   binary script's unused `covar` argument.
5. Near-zero median coefficients have a documented beta_floor of 1e-8 for
   finite adaptive weights. Five-decimal selection matches the binary script;
   the legacy continuous known-covariate branch used eight decimals.
6. `lqa::lqa.default` receives an explicitly constructed intercept/design
   matrix, avoiding the archived formula method's non-list contrasts warning.
   The LQA solver and n^exponent penalty parameterization are retained.
7. Continuous tuning uses weighted Pearson correlation implemented with
   stats::cov.wt rather than adding the wCorr dependency.
8. Failed candidates are recorded with errors and excluded with a warning.
   Non-finite weights, constant columns and rank-deficient full models stop
   explicitly. Binary propensity clipping is opt-in and changes weights.
9. Empty adjustment sets retain unit weights, matching the supplied scripts.
10. Survey SEs are labeled conditional; no bootstrap or full published
    simulation reproduction is claimed. Synthetic examples use a smaller,
    transparent log-normal error design.

## Validation scope

Run `R CMD build GMAL` and `R CMD check --no-manual GMAL_0.1.0.tar.gz`.
The check includes independent binary reference calculations, continuous
balance calculations, forced inclusion, input errors, and response invariance.
See the adjacent delivery `validation` folder for actual logs and environment.
Mac validation does not by itself establish Windows/Linux compatibility;
GitHub Actions exercises these platforms after publication.
