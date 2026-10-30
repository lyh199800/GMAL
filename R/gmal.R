# Algorithm: Liu et al. (2024), doi:10.1093/bib/bbae059.
# Source mapping and implementation changes are documented in DEVELOPMENT.md.

#' Generalized median adaptive lasso
#' See the installed help page ?gmal for arguments and statistical limitations.
gmal <- function(data, covariates, treatment, outcome,
                 treatment_type = c("binary", "continuous"), known = NULL,
                 lambda_exponents = c(-10, -5, -2, -1, -.75, -.5, -.25, .25, .49),
                 gamma_convergence = 2, beta_floor = 1e-8,
                 selection_digits = 5, ps_clip = NULL) {
  call <- match.call()
  treatment_type <- match.arg(treatment_type)
  if (!is.data.frame(data)) stop("data must be a data.frame.")
  if (!is.character(covariates) || !length(covariates) || anyDuplicated(covariates))
    stop("covariates must be distinct column names.")
  if (!is.null(known) && (!is.character(known) || anyDuplicated(known)))
    stop("known must be distinct column names.")
  if (!is.character(treatment) || length(treatment) != 1L ||
      !is.character(outcome) || length(outcome) != 1L)
    stop("treatment and outcome must each name one column.")
  cols <- c(known, covariates, treatment, outcome)
  if (anyDuplicated(cols) || anyDuplicated(names(data)) || !all(cols %in% names(data)))
    stop("Use existing, unique, non-overlapping column names.")
  if (!all(vapply(data[cols], is.numeric, logical(1))) ||
      any(!is.finite(as.matrix(data[cols]))))
    stop("Analysis columns must be numeric and finite; handle missing values explicitly.")
  if (!is.numeric(lambda_exponents) || !length(lambda_exponents) ||
      any(!is.finite(lambda_exponents)) || anyDuplicated(lambda_exponents))
    stop("lambda_exponents must contain distinct finite numbers.")
  if (length(gamma_convergence) != 1L || !is.finite(gamma_convergence) ||
      any(2 * (gamma_convergence - lambda_exponents + 1) <= 0))
    stop("gamma_convergence must yield positive adaptive exponents.")
  if (length(beta_floor) != 1L || !is.finite(beta_floor) || beta_floor <= 0)
    stop("beta_floor must be positive and finite.")
  if (length(selection_digits) != 1L || !is.finite(selection_digits) ||
      selection_digits < 0 || selection_digits > 12 || selection_digits %% 1 != 0)
    stop("selection_digits must be an integer from 0 to 12.")
  if (!is.null(ps_clip) && (length(ps_clip) != 1L || !is.finite(ps_clip) ||
                           ps_clip <= 0 || ps_clip >= .5))
    stop("ps_clip must be NULL or a number strictly between 0 and 0.5.")
  if (treatment_type == "continuous" && !is.null(ps_clip))
    stop("ps_clip applies only to binary treatments.")
  xnames <- c(known, covariates)
  X <- as.matrix(data[xnames]); A <- data[[treatment]]; Y <- data[[outcome]]
  n <- nrow(X)
  if (n <= ncol(X) + 2L)
    stop("The unpenalized median model requires n > number of covariates + 2; screen first.")
  if (treatment_type == "binary" && !setequal(unique(A), c(0, 1)))
    stop("Binary treatment must contain both 0 and 1.")
  if (length(unique(A)) < 2L) stop("Treatment must vary.")
  scales <- apply(X, 2, stats::sd)
  if (any(scales <= 0)) stop("Remove constant covariates.")
  centers <- colMeans(X)
  X <- scale(X, center = centers, scale = scales)
  if (qr(cbind(1, A, X))$rank < ncol(X) + 2L)
    stop("Outcome design is rank deficient; remove collinear covariates.")
  # Internal names allow non-syntactic user column names without formula injection.
  internal <- paste0("x", seq_len(ncol(X)))
  d <- as.data.frame(X); names(d) <- internal; d$a <- A; d$y <- Y
  median_fit <- quantreg::rq(stats::reformulate(c("a", internal), "y"),
                             data = d, tau = .5)
  beta <- stats::coef(median_fit)[internal]
  if (any(!is.finite(beta))) stop("Median regression produced non-finite coefficients.")
  names(beta) <- xnames
  k <- length(known)
  fam <- if (treatment_type == "binary") stats::binomial() else stats::gaussian()
  design_matrix <- stats::model.matrix(stats::reformulate(internal), data = d)
  if (treatment_type == "continuous" && !requireNamespace("CBPS", quietly = TRUE))
    stop("Install CBPS for continuous treatment: install.packages('CBPS').")
  paths <- vector("list", length(lambda_exponents))
  rows <- vector("list", length(paths))
  for (i in seq_along(paths)) {
    exponent <- lambda_exponents[i]
    gamma <- 2 * (gamma_convergence - exponent + 1)
    lambda <- n^exponent
    result <- tryCatch({
      penalty_weights <- exp(-gamma * log(pmax(abs(beta), beta_floor)))
      if (k) penalty_weights[seq_len(k)] <- 0
      if (any(!is.finite(penalty_weights)) || !is.finite(lambda) || lambda <= 0)
        stop("Non-finite adaptive penalty; adjust the tuning grid or beta_floor.")
      penalty <- lqa::adaptive.lasso(lambda = lambda, al.weights = penalty_weights)
      selection_fit <- lqa::lqa.default(x = design_matrix, y = A,
        penalty = penalty, family = fam)
      b <- stats::coef(selection_fit)[internal]
      if (any(!is.finite(b))) stop("Non-finite treatment coefficients.")
      keep <- which(round(b, selection_digits) != 0)
      if (k) keep <- sort(unique(c(seq_len(k), keep)))
      if (!length(keep)) {
        w <- rep(1, n); ps <- NULL
      } else if (treatment_type == "binary") {
        ps_fit <- stats::glm(stats::reformulate(internal[keep], "a"), data = d,
                             family = stats::binomial())
        if (!isTRUE(ps_fit$converged)) stop("Propensity refit did not converge.")
        ps <- stats::predict(ps_fit, type = "response")
        if (!is.null(ps_clip)) ps <- pmin(1 - ps_clip, pmax(ps_clip, ps))
        if (any(!is.finite(ps)) || any(ps <= 0 | ps >= 1))
          stop("Boundary propensity scores; inspect overlap or explicitly set ps_clip.")
        w <- ifelse(A == 1, 1 / ps, 1 / (1 - ps))
      } else {
        cb <- CBPS::npCBPS(stats::reformulate(internal[keep], "a"), data = d,
                            corprior = .1 / n, print.level = 0)
        w <- as.numeric(cb$weights); ps <- NULL
      }
      if (length(w) != n || any(!is.finite(w)) || any(w <= 0))
        stop("Invalid balancing weights.")
      balance <- .gmal_balance(X, A, w, treatment_type)
      score <- sum(abs(beta) * balance)
      if (!is.finite(score)) stop("Non-finite balance score.")
      list(selected = xnames[keep], weights = w, propensity = ps,
           balance = balance, score = score, treatment_coef = stats::setNames(b, xnames))
    }, error = function(e) list(error = conditionMessage(e)))
    paths[[i]] <- result
    rows[[i]] <- data.frame(exponent = exponent, lambda = lambda, gamma = gamma,
      score = if (is.null(result$error)) result$score else Inf,
      n_selected = if (is.null(result$error)) length(result$selected) else NA_integer_,
      error = if (is.null(result$error)) "" else result$error)
  }
  tuning <- do.call(rbind, rows)
  good <- is.finite(tuning$score)
  if (!any(good)) stop("All tuning candidates failed: ", paste(unique(tuning$error), collapse = "; "))
  if (any(!good)) warning(sum(!good), " tuning candidates failed; inspect $tuning.")
  best <- which.min(tuning$score); chosen <- paths[[best]]
  weighted_data <- data.frame(y = Y, a = A, w = chosen$weights)
  design <- survey::svydesign(ids = ~1, weights = ~w, data = weighted_data)
  fit <- survey::svyglm(y ~ a, design = design)
  result <- list(call = call, treatment_type = treatment_type,
    estimate = unname(stats::coef(fit)["a"]),
    std_error = sqrt(stats::vcov(fit)["a", "a"]), selected = chosen$selected,
    lambda_exponent = tuning$exponent[best], lambda = tuning$lambda[best],
    weights = chosen$weights, propensity = chosen$propensity,
    balance = data.frame(covariate = xnames,
      unweighted = as.numeric(.gmal_balance(X, A, rep(1, n), treatment_type)),
      weighted = as.numeric(chosen$balance), outcome_coefficient = as.numeric(beta)),
    tuning = tuning, path = paths, fit = fit, outcome_model = median_fit,
    centers = centers, scales = scales, ps_clip = ps_clip,
    effective_sample_size = sum(chosen$weights)^2 / sum(chosen$weights^2))
  class(result) <- "gmal"
  result
}

.gmal_balance <- function(X, A, w, type) {
  if (type == "binary") {
    apply(X, 2, function(x) abs(stats::weighted.mean(x[A == 1], w[A == 1]) -
                               stats::weighted.mean(x[A == 0], w[A == 0])))
  } else {
    apply(X, 2, function(x) abs(stats::cov.wt(cbind(x, A), wt = w, cor = TRUE)$cor[1, 2]))
  }
}

print.gmal <- function(x, ...) {
  cat("GMAL (", x$treatment_type, " treatment)\n", sep = "")
  cat("Effect estimate:", format(x$estimate, digits = 5),
      "  Conditional SE:", format(x$std_error, digits = 4), "\n")
  cat("Selected:", if (length(x$selected)) paste(x$selected, collapse = ", ") else "(none)", "\n")
  cat("Lambda exponent:", x$lambda_exponent,
      "  Effective sample size:", round(x$effective_sample_size, 1), "\n")
  invisible(x)
}

coef.gmal <- function(object, ...) stats::coef(object$fit)

predict.gmal <- function(object, newdata, ...) {
  if (!is.numeric(newdata) || any(!is.finite(newdata)))
    stop("newdata must be a finite numeric vector of treatment values.")
  b <- stats::coef(object$fit)
  unname(b[1] + b[2] * newdata)
}
