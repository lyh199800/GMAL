\name{gmal}
\alias{gmal}
\alias{print.gmal}
\alias{coef.gmal}
\alias{predict.gmal}
\title{Generalized Median Adaptive Lasso for Causal Inference}
\usage{
gmal(data, covariates, treatment, outcome,
  treatment_type = c("binary", "continuous"), known = NULL,
  lambda_exponents = c(-10, -5, -2, -1, -.75, -.5, -.25, .25, .49),
  gamma_convergence = 2, beta_floor = 1e-8,
  selection_digits = 5, ps_clip = NULL)
\method{print}{gmal}(x, ...)
\method{coef}{gmal}(object, ...)
\method{predict}{gmal}(object, newdata, ...)
}
\arguments{
\item{data}{Data frame with complete, finite numeric analysis columns.}
\item{covariates}{Character vector of candidate covariate names.}
\item{treatment, outcome}{Single treatment and outcome column names. Binary treatment is coded 0/1.}
\item{treatment_type}{Explicit treatment type; continuous treatment requires CBPS.}
\item{known}{Optional covariate names forced into all adjustment sets with zero penalty. Must not overlap candidates.}
\item{lambda_exponents}{Grid of exponents \eqn{e} defining \eqn{\lambda=n^e}.}
\item{gamma_convergence}{Defines adaptive exponent \eqn{\gamma=2(c-e+1)}.}
\item{beta_floor}{Positive lower bound on absolute median coefficients when constructing penalty weights.}
\item{selection_digits}{Coefficients that round to zero at this many decimal places are excluded.}
\item{ps_clip}{Optional explicit binary propensity clipping threshold. Default NULL performs no clipping.}
\item{x, object}{A fitted gmal object.}
\item{newdata}{Numeric vector of treatment values for predicting the linear marginal response.}
\item{...}{Reserved additional arguments.}
}
\description{Uses median regression outcome coefficients to weight adaptive lasso treatment-model penalties and selects tuning by covariate balance.}
\details{
Covariates, including known covariates, are standardized using sample standard deviations.
An unpenalized median outcome regression includes treatment and all covariates.
Treatment-model selection uses lqa with a binomial family for binary treatment or Gaussian family for continuous treatment.
Selected binary covariates enter an unpenalized logistic refit for inverse probability weights. If no covariates are selected, unit weights are used, as in the supplied research scripts.
Continuous treatment weights use \code{CBPS::npCBPS} with \code{corprior=.1/n}.
The binary score is the sum of absolute weighted covariate mean differences times absolute median outcome coefficients.
The continuous score uses absolute weighted Pearson correlations with treatment.
The first minimum in the supplied tuning grid wins. Failed candidates are recorded and excluded, with a warning; if all fail the function stops.
The reported standard error is conditional on selection and estimated weights. It does not include selection or weight-estimation uncertainty. For research inference consider resampling the entire pipeline; published application intervals used bootstrap resampling.
The estimand is a mean treatment effect (binary) or a linear marginal dose-response slope (continuous), not a median treatment effect.
The full median model requires more observations than covariates plus two and a full-rank design. For larger candidate sets, apply and document a separate screening procedure first.
Causal interpretation additionally requires consistency, no unmeasured confounding, adequate overlap and an appropriate marginal response model.
}
\value{A gmal object containing \code{estimate}, conditional \code{std_error}, selected covariate names in \code{selected}, selected \code{lambda} and \code{lambda_exponent}, observation \code{weights}, binary \code{propensity} where available, covariate \code{balance}, all candidate \code{tuning} scores and errors, candidate \code{path}, weighted survey \code{fit}, median \code{outcome_model}, covariate \code{centers} and \code{scales}, and \code{effective_sample_size}. Coef returns intercept and treatment slope; predict returns linear marginal response predictions.}
\references{Liu Y et al. (2024). High-dimensional generalized median adaptive lasso with application to omics data. Briefings in Bioinformatics 25(2), bbae059. \doi{10.1093/bib/bbae059}.}
\examples{
set.seed(2026)
d <- simulate_gmal(n = 200, p = 6)
f <- gmal(d, paste0("X", 1:6), "A", "Y")
f
f$balance
predict(f, c(0, 1))
}
