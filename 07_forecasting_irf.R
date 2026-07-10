# ==============================================================================
# 03 - VECTOR ERROR CORRECTION MODEL (VECM)
# ==============================================================================
# source("R/00_setup_data.R") and R/02_stationarity_cointegration.R must be
# run first (this script needs ect_gasoline / ect_diesel).
#
# The original analysis was run separately with 3 and 5 lags of the
# differenced series as a robustness check. Both are produced here from a
# single function instead of two near-identical scripts.
# ==============================================================================

# ------------------------------------------------------------------------
# 1. VECM estimator: delta_fuel = f(delta_oil, lags, ECT(t-1))
# ------------------------------------------------------------------------
run_vecm <- function(d_log_fuel, d_log_oil, ect, fuel_name, n_lags) {

  ect_lag <- ect[1:(length(ect) - 1)]  # align ECT(t-1) with differenced series

  lags_oil  <- create_lags(d_log_oil,  n_lags)
  lags_fuel <- create_lags(d_log_fuel, n_lags)

  model_data <- na.omit(data.frame(
    d_log_fuel = d_log_fuel,
    d_log_oil  = d_log_oil,
    lags_oil,
    lags_fuel,
    ect_lag    = ect_lag
  ))
  colnames(model_data)[3:(2 + n_lags)]                 <- paste0("oil_lag",  1:n_lags)
  colnames(model_data)[(3 + n_lags):(2 + 2 * n_lags)]   <- paste0("fuel_lag", 1:n_lags)

  formula <- as.formula(paste(
    "d_log_fuel ~ d_log_oil +",
    paste(paste0("oil_lag", 1:n_lags), collapse = " + "), "+",
    paste(paste0("fuel_lag", 1:n_lags), collapse = " + "), "+ ect_lag"
  ))

  model  <- lm(formula, data = model_data)
  robust <- coeftest(model, vcov = vcovHAC(model))

  alpha       <- coef(model)["ect_lag"]
  speed_adj   <- abs(alpha) * 100
  half_life   <- log(0.5) / log(1 + alpha)
  sr_elast    <- coef(model)["d_log_oil"]

  cat(sprintf("\n=== %s VECM (%d lags) ===\n", toupper(fuel_name), n_lags))
  print(robust)
  cat(sprintf("Error correction coefficient (alpha): %.4f\n", alpha))
  cat(sprintf("Speed of adjustment: %.2f%% per week\n", speed_adj))
  cat(sprintf("Half-life: %.1f weeks\n", half_life))
  cat(sprintf("R-squared: %.4f | Adj. R-squared: %.4f\n",
              summary(model)$r.squared, summary(model)$adj.r.squared))

  list(model = model, robust = robust, alpha = alpha, speed_adj = speed_adj,
       half_life = half_life, sr_elasticity = sr_elast, n_obs = nrow(model_data))
}

# ------------------------------------------------------------------------
# 2. Estimate for both lag specifications (3 and 5), both fuels
# ------------------------------------------------------------------------
vecm_gasoline_3 <- run_vecm(d_log_gasoline, d_log_brent_3M, ect_gasoline, "Gasoline", 3)
vecm_diesel_3   <- run_vecm(d_log_diesel,   d_log_brent_3M, ect_diesel,   "Diesel",   3)

vecm_gasoline_5 <- run_vecm(d_log_gasoline, d_log_brent_3M, ect_gasoline, "Gasoline", 5)
vecm_diesel_5   <- run_vecm(d_log_diesel,   d_log_brent_3M, ect_diesel,   "Diesel",   5)

# ------------------------------------------------------------------------
# 3. Summary table (short-run vs long-run elasticity, speed of adjustment)
# ------------------------------------------------------------------------
build_vecm_summary <- function(res_gas, res_diesel, n_lags) {
  data.frame(
    Fuel                  = c("Gasoline", "Diesel"),
    Lags                  = n_lags,
    Short_Run_Elasticity  = round(c(res_gas$sr_elasticity, res_diesel$sr_elasticity), 4),
    Long_Run_Elasticity   = round(c(coint_gasoline$beta, coint_diesel$beta), 4),
    ECM_Coefficient       = round(c(res_gas$alpha, res_diesel$alpha), 4),
    Speed_of_Adjustment_% = round(c(res_gas$speed_adj, res_diesel$speed_adj), 2),
    Half_Life_Weeks       = round(c(res_gas$half_life, res_diesel$half_life), 1),
    R_squared             = round(c(summary(res_gas$model)$r.squared,
                                     summary(res_diesel$model)$r.squared), 4)
  )
}

vecm_summary_3lag <- build_vecm_summary(vecm_gasoline_3, vecm_diesel_3, 3)
vecm_summary_5lag <- build_vecm_summary(vecm_gasoline_5, vecm_diesel_5, 5)

print(vecm_summary_3lag)
print(vecm_summary_5lag)

# ------------------------------------------------------------------------
# 4. Export LaTeX tables
# ------------------------------------------------------------------------
export_vecm_table <- function(summary_df, n_lags) {
  print(xtable(summary_df,
               caption = sprintf("Vector Error Correction Model Results (%d Lags)", n_lags),
               label = sprintf("tab:vecm_summary_%dlag", n_lags),
               align = c("l", "l", "c", rep("r", 5)),
               digits = c(0, 0, 0, 4, 4, 4, 2, 1)),
        include.rownames = FALSE, caption.placement = "top", booktabs = TRUE,
        file = file.path("output", sprintf("table_vecm_summary_%dlag.tex", n_lags)))
}

export_vecm_table(vecm_summary_3lag, 3)
export_vecm_table(vecm_summary_5lag, 5)

cat("\nVECM tables exported to output/table_vecm_summary_3lag.tex and _5lag.tex\n")
