# ==============================================================================
# 02 - STATIONARITY, LAG SELECTION & COINTEGRATION
# ==============================================================================
# source("R/00_setup_data.R") must be run first.
# ==============================================================================

# ------------------------------------------------------------------------
# 1. Augmented Dickey-Fuller unit root tests (levels, log levels)
# ------------------------------------------------------------------------
run_adf <- function(series, name, type = "drift", lags = 3) {
  test <- ur.df(na.omit(series), type = type, lags = lags, selectlags = "Fixed")
  cat(sprintf("ADF (%s): %-20s stat = %6.3f\n", type, name, test@teststat[1]))
  test
}

adf_gasoline <- run_adf(log_gasoline, "Log Gasoline")
adf_diesel   <- run_adf(log_diesel,   "Log Diesel")
adf_brent    <- run_adf(log_brent_3M, "Log Brent 3M")

# All three series are expected to be I(1): non-stationary in levels,
# stationary in first differences - a prerequisite for testing cointegration.
adf_d_gasoline <- run_adf(d_log_gasoline, "D Log Gasoline")
adf_d_diesel   <- run_adf(d_log_diesel,   "D Log Diesel")
adf_d_brent    <- run_adf(d_log_brent_3M, "D Log Brent 3M")

# ------------------------------------------------------------------------
# 2. VAR lag-order selection (AIC / BIC) for each oil-fuel system
# ------------------------------------------------------------------------
select_var_lags <- function(y, x, lag_max = 12) {
  sel <- VARselect(cbind(x, y), lag.max = lag_max, type = "const")
  data.frame(
    Lag = 1:lag_max,
    AIC = sel$criteria["AIC(n)", ],
    BIC = sel$criteria["SC(n)", ]
  )
}

lag_selection_gasoline <- select_var_lags(d_log_gasoline, d_log_brent_3M)
lag_selection_diesel   <- select_var_lags(d_log_diesel,   d_log_brent_3M)

print(lag_selection_gasoline)
print(lag_selection_diesel)

# ------------------------------------------------------------------------
# 3. Engle-Granger cointegration test: is there a stable long-run
#    relationship between oil and each retail fuel price?
# ------------------------------------------------------------------------
run_cointegration_test <- function(log_fuel, log_oil, fuel_name, lags = 3) {
  coint_reg <- lm(log_fuel ~ log_oil)
  coeftest(coint_reg, vcov = vcovHAC(coint_reg))

  resid_reg <- residuals(coint_reg)
  adf_resid <- ur.df(resid_reg, type = "drift", lags = lags, selectlags = "Fixed")

  cat(sprintf("\n=== %s COINTEGRATION TEST ===\n", toupper(fuel_name)))
  cat("ADF test statistic:", adf_resid@teststat[1], "\n")
  cat("Engle-Granger critical values: 1% = -3.90 | 5% = -3.34 | 10% = -3.04\n")
  cat("Long-run elasticity (beta):", coef(coint_reg)[2], "\n")

  list(model = coint_reg, residuals = resid_reg, adf = adf_resid,
       beta = coef(coint_reg)[2], r_squared = summary(coint_reg)$r.squared)
}

coint_gasoline <- run_cointegration_test(log_gasoline, log_brent_3M, "Gasoline")
coint_diesel   <- run_cointegration_test(log_diesel,   log_brent_3M, "Diesel")

# Error correction terms (equilibrium deviations), used as ECT(t-1) in the
# VECM / asymmetric ECM specifications in the next scripts.
ect_gasoline <- coint_gasoline$residuals
ect_diesel   <- coint_diesel$residuals

# ------------------------------------------------------------------------
# 4. Summary table + residual diagnostics plot
# ------------------------------------------------------------------------
cointegration_summary <- data.frame(
  Relationship = c("Oil -> Gasoline", "Oil -> Diesel"),
  Beta         = c(coint_gasoline$beta, coint_diesel$beta),
  R_squared    = c(coint_gasoline$r_squared, coint_diesel$r_squared),
  ADF_stat     = c(coint_gasoline$adf@teststat[1], coint_diesel$adf@teststat[1]),
  Cointegrated = c(
    ifelse(coint_gasoline$adf@teststat[1] < -3.34, "Yes (5%)", "No"),
    ifelse(coint_diesel$adf@teststat[1]   < -3.34, "Yes (5%)", "No")
  )
)
print(cointegration_summary)

par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
plot(ect_gasoline, type = "l", main = "Gasoline: Cointegration Residuals",
     ylab = "Residuals", xlab = "Time"); abline(h = 0, col = "red", lty = 2)
acf(ect_gasoline, main = "Gasoline: ACF of Residuals")
plot(ect_diesel, type = "l", main = "Diesel: Cointegration Residuals",
     ylab = "Residuals", xlab = "Time"); abline(h = 0, col = "red", lty = 2)
acf(ect_diesel, main = "Diesel: ACF of Residuals")
par(mfrow = c(1, 1))

cat("\nCointegration residuals saved as 'ect_gasoline' and 'ect_diesel'\n")
