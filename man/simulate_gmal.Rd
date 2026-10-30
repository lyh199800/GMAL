\name{simulate_gmal}
\alias{simulate_gmal}
\title{Generate a Small Simulated GMAL Example}
\usage{
simulate_gmal(n = 300L, p = 10L,
  treatment_type = c("binary", "continuous"), effect = 2)
}
\arguments{
\item{n}{Integer sample size at least 20.}
\item{p}{Integer number of covariates at least 4.}
\item{treatment_type}{Binary or continuous.}
\item{effect}{True additive treatment effect.}
}
\description{Generates independent normal covariates with X1 a confounder, X2 prognostic, X3 an instrument and remaining covariates spurious. Errors are centered log-normal. This is an illustrative example, not the full published simulation design. The caller controls the random seed.}
\value{A data frame with columns X1 through Xp, A (treatment), and Y (outcome).}
\examples{
set.seed(2026)
d <- simulate_gmal()
head(d)
}
