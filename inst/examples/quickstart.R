library(GMAL)
set.seed(2026)
d <- simulate_gmal()
fit <- gmal(d, paste0("X", 1:10), "A", "Y")
print(fit)
fit$balance
fit$tuning
predict(fit, c(0, 1))

# Continuous treatment requires CBPS.
set.seed(2026)
dc <- simulate_gmal(treatment_type = "continuous")
fc <- gmal(dc, paste0("X", 1:10), "A", "Y", treatment_type = "continuous")
print(fc)
