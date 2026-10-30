simulate_gmal <- function(n = 300L, p = 10L,
                          treatment_type = c("binary", "continuous"), effect = 2) {
  treatment_type <- match.arg(treatment_type)
  if (length(n) != 1L || !is.finite(n) || n < 20 || n %% 1 != 0 ||
      length(p) != 1L || !is.finite(p) || p < 4 || p %% 1 != 0 ||
      length(effect) != 1L || !is.finite(effect)) stop("Use integer n >= 20, p >= 4 and finite effect.")
  X <- matrix(stats::rnorm(n * p), n, p)
  a <- if (treatment_type == "binary") {
    stats::rbinom(n, 1, stats::plogis(.6 * X[, 1] + .6 * X[, 3]))
  } else .6 * X[, 1] + .6 * X[, 3] + stats::rnorm(n)
  # Centered log-normal errors demonstrate skewed outcomes.
  y <- effect * a + X[, 1] + X[, 2] + stats::rlnorm(n, 0, .5) - exp(.125)
  d <- as.data.frame(X); names(d) <- paste0("X", seq_len(p))
  d$A <- a; d$Y <- y
  d
}
