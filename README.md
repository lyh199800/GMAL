# GMAL

**Generalized Median Adaptive Lasso for Causal Inference**

An R implementation of the method introduced in:

Yahang Liu, Qian Gao, Kecheng Wei, Chen Huang, Ce Wang, Yongfu Yu,
Guoyou Qin, and Tong Wang (2024).
[High-dimensional generalized median adaptive lasso with application to omics data](https://doi.org/10.1093/bib/bbae059).
*Briefings in Bioinformatics*, 25(2), bbae059.

GMAL uses median outcome regression to construct adaptive penalty weights,
selects covariates by covariate balance, and estimates a weighted marginal
mean treatment effect. It supports **binary and continuous treatments**.

This is a newly packaged research implementation developed from local research
materials. It is not on CRAN and does not claim to reproduce every published
simulation or the restricted ADNI application. See [implementation notes](DEVELOPMENT.md).

## Installation

R >= 4.1 is required. `lqa` is archived on CRAN and must be installed first:

```r
install.packages(c("quantreg", "survey", "CBPS", "remotes"))
install.packages(
  "https://cran.r-project.org/src/contrib/Archive/lqa/lqa_1.0-3.tar.gz",
  repos = NULL, type = "source"
)
# After this repository has been published:
remotes::install_github("lyh199800/GMAL")
```

For a downloaded source archive:

```r
install.packages("GMAL_0.1.0.tar.gz", repos = NULL, type = "source")
```

## Binary treatment

```r
library(GMAL)
set.seed(2026)
d <- simulate_gmal(n = 300, p = 10)
fit <- gmal(d, covariates = paste0("X", 1:10), treatment = "A", outcome = "Y")
fit
fit$selected
fit$estimate
fit$balance
fit$tuning
predict(fit, c(0, 1))
```

## Continuous treatment

```r
set.seed(2026)
d <- simulate_gmal(n = 300, p = 10, treatment_type = "continuous")
fit <- gmal(d, paste0("X", 1:10), "A", "Y", treatment_type = "continuous")
fit
```

Continuous balancing uses `CBPS::npCBPS`, matching the method's balancing
approach. Covariates known in advance can be forced into every adjustment set:

```r
fit <- gmal(d, paste0("X", 2:10), "A", "Y",
            treatment_type = "continuous", known = "X1")
```

## Inputs and interpretation

- Provide complete numeric covariates and outcome; code binary exposure as 0/1.
- Expand categorical covariates to numeric indicator columns before fitting.
- The unpenalized median outcome model requires a full-rank design and
  `n > number of covariates + 2`. Screen larger candidate sets first.
- Inspect covariate balance, weights, effective sample size and failed tuning
  candidates. `ps_clip` is an explicit optional binary propensity threshold;
  it is disabled by default.
- Median regression constructs selection weights. The final estimand is a
  **mean effect or linear dose-response slope**, not a quantile effect.
- `std_error` treats selected covariates and estimated weights as fixed.
  It excludes selection and weight-estimation uncertainty. For inference,
  consider bootstrapping the entire pipeline and report failed resamples.
- Causal interpretation requires consistency, no unmeasured confounding,
  positivity and an appropriate marginal response model.

## Citation and provenance

Run `citation("GMAL")` for the full citation. GMAL builds on outcome-adaptive
lasso and generalized outcome-adaptive lasso (GOAL). The supplied legacy GOAL
package credits Qian Gao; that package and third-party comparison methods are
not redistributed as new GMAL code here. All paper authors are credited in
the citation. This repository contains method code and simulated examples only.

## Contact

Yahang Liu · [yhliu22@m.fudan.edu.cn](mailto:yhliu22@m.fudan.edu.cn)

License: GPL-3. Package version: 0.1.0.
