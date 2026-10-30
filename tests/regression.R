library(GMAL)
expect_error <- function(expr) {
  failed <- inherits(try(force(expr), silent = TRUE), "try-error")
  stopifnot(failed)
}
set.seed(2026)
d <- simulate_gmal(n = 240, p = 6)
xn <- paste0("X", 1:6)
grid <- c(-2, -.5, .25)
f <- gmal(d, xn, "A", "Y", lambda_exponents = grid)
stopifnot(is.finite(f$estimate), all(f$weights > 0),
          length(f$weights) == nrow(d), nrow(f$tuning) == length(grid),
          abs(diff(predict(f, c(0, 1))) - f$estimate) < 1e-12)

# Independently reconstruct the supplied binary algorithm on a safe data set.
X <- scale(as.matrix(d[xn])); dd <- as.data.frame(X); dd$A <- d$A; dd$Y <- d$Y
med <- quantreg::rq(stats::reformulate(c("A", xn), "Y"), data = dd, tau = .5)
beta <- coef(med)[xn]
ref_scores <- numeric(length(grid)); ref_sets <- vector("list", length(grid))
ref_weights <- vector("list", length(grid))
for (i in seq_along(grid)) {
  penalty <- lqa::adaptive.lasso(lambda = nrow(d)^grid[i],
    al.weights = abs(beta)^(-2 * (2 - grid[i] + 1)))
  # Formula interface matches the source; suppress its known contrasts warning.
  sel <- suppressWarnings(lqa::lqa.formula(stats::reformulate(xn, "A"),
                         data = dd, penalty = penalty, family = binomial()))
  chosen <- xn[round(coef(sel)[xn], 5) != 0]
  if (length(chosen)) {
    ps <- predict(glm(stats::reformulate(chosen, "A"), data = dd,
                      family = binomial()), type = "response")
    w <- ifelse(dd$A == 1, 1 / ps, 1 / (1 - ps))
  } else w <- rep(1, nrow(d))
  diffs <- vapply(xn, function(j) {
    abs(sum(dd[dd$A == 1, j] * w[dd$A == 1]) / sum(w[dd$A == 1]) -
        sum(dd[dd$A == 0, j] * w[dd$A == 0]) / sum(w[dd$A == 0]))
  }, numeric(1))
  ref_scores[i] <- sum(diffs * abs(beta)); ref_sets[[i]] <- chosen
  ref_weights[[i]] <- w
}
ib <- which.min(ref_scores)
stopifnot(isTRUE(all.equal(f$tuning$score, ref_scores, tolerance = 1e-8)),
          identical(f$selected, ref_sets[[ib]]),
          isTRUE(all.equal(f$weights, unname(ref_weights[[ib]]), tolerance = 1e-8)))
refd <- data.frame(y = d$Y, a = d$A, w = ref_weights[[ib]])
refdesign <- survey::svydesign(ids = ~1, weights = ~w, data = refd)
ref_effect <- coef(survey::svyglm(y ~ a, design = refdesign))["a"]
stopifnot(abs(f$estimate - ref_effect) < 1e-8)

# Known variables stay selected; non-syntactic names work; globals cannot leak.
known <- gmal(d, xn[-1], "A", "Y", known = "X1", lambda_exponents = grid)
stopifnot("X1" %in% known$selected)
renamed <- d; names(renamed)[1] <- "age in years"
g <- gmal(renamed, c("age in years", xn[-1]), "A", "Y", lambda_exponents = grid)
stopifnot(abs(f$estimate - g$estimate) < 1e-8)
Y <- rep(-999, nrow(d)); tau <- .1
fg <- gmal(d, xn, "A", "Y", lambda_exponents = grid)
stopifnot(abs(f$estimate - fg$estimate) < 1e-12)
dshift <- d; dshift$Y <- dshift$Y + 10
fs <- gmal(dshift, xn, "A", "Y", lambda_exponents = grid)
stopifnot(abs(f$estimate - fs$estimate) < 1e-6)

bad <- d; bad$X1[1] <- NA
expect_error(gmal(bad, xn, "A", "Y"))
bad <- d; bad$X1 <- 1
expect_error(gmal(bad, xn, "A", "Y"))
bad <- d; bad$X2 <- bad$X1
expect_error(gmal(bad, xn, "A", "Y"))
bad <- d; bad$A <- 1
expect_error(gmal(bad, xn, "A", "Y"))
expect_error(gmal(d[1:7, ], xn, "A", "Y"))
expect_error(gmal(d, xn, "A", "Y", known = "X1"))
expect_error(gmal(d, xn, "A", "Y", ps_clip = .7))

if (requireNamespace("CBPS", quietly = TRUE)) {
  set.seed(42)
  dc <- simulate_gmal(n = 240, p = 6, treatment_type = "continuous")
  fc <- gmal(dc, xn, "A", "Y", treatment_type = "continuous", lambda_exponents = grid)
  stopifnot(is.finite(fc$estimate), all(fc$weights > 0), all(is.finite(fc$balance$weighted)))
  # Explicit weighted correlation calculation independent of cov.wt.
  xx <- scale(as.matrix(dc[xn])); a <- dc$A; w <- fc$weights / sum(fc$weights)
  ac <- a - sum(w * a)
  corr <- apply(xx, 2, function(x) {
    xc <- x - sum(w * x)
    abs(sum(w * xc * ac) / sqrt(sum(w * xc^2) * sum(w * ac^2)))
  })
  stopifnot(max(abs(unname(corr) - fc$balance$weighted)) < 1e-10)
  fck <- gmal(dc, xn[-1], "A", "Y", treatment_type = "continuous", known = "X1", lambda_exponents = grid)
  stopifnot("X1" %in% fck$selected)
}
cat("All GMAL regression checks passed.\n")
