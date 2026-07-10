# ==============================================================================
# 07 - VAR, IMPULSE RESPONSE FUNCTIONS & FORECASTING
# ==============================================================================
# source("R/00_setup_data.R") and R/02_stationarity_cointegration.R must be
# run first.
#
# Traces out how a 1% shock to Brent crude propagates to retail gasoline and
# diesel prices over a 20-week horizon, using a VAR(5) estimated on the
# stationary (differenced) series - the lag order chosen in
# R/02_stationarity_cointegration.R via AIC/BIC.
# ==============================================================================

# ------------------------------------------------------------------------
# 1. Estimate VAR(5) systems and compute bootstrapped IRFs
# ------------------------------------------------------------------------
run_irf <- function(d_log_fuel, d_log_oil, fuel_name, n_ahead = 20, n_lags = 5) {
  var_data <- na.omit(cbind(oil = d_log_oil, fuel = d_log_fuel))
  var_model <- VAR(var_data, p = n_lags, type = "const")

  irf_result <- irf(var_model, impulse = "oil", response = "fuel",
                     n.ahead = n_ahead, boot = TRUE, runs = 500, ci = 0.95)

  cat(sprintf("VAR(%d) estimated for Oil -> %s, IRF computed over %d weeks\n",
              n_lags, fuel_name, n_ahead))

  list(var = var_model, irf = irf_result, n_obs = nrow(var_data))
}

irf_gas    <- run_irf(d_log_gasoline, d_log_brent_3M, "Gasoline")
irf_diesel <- run_irf(d_log_diesel,   d_log_brent_3M, "Diesel")

# ------------------------------------------------------------------------
# 2. Extract responses at key horizons (immediate, 1 month, 1 quarter, ...)
# ------------------------------------------------------------------------
extract_irf_values <- function(irf_obj, weeks) {
  idx <- weeks + 1  # R is 1-indexed; week 0 = impact period
  data.frame(
    week = weeks,
    response = irf_obj$irf[[1]][idx],
    lower    = irf_obj$Lower[[1]][idx],
    upper    = irf_obj$Upper[[1]][idx]
  ) |> transform(significant = ifelse(lower > 0 | upper < 0, "Yes", "No"))
}

key_weeks <- c(0, 1, 2, 4, 8, 12, 20)
gas_key_values    <- extract_irf_values(irf_gas$irf,    key_weeks)
diesel_key_values <- extract_irf_values(irf_diesel$irf, key_weeks)

cat("\nGasoline response to a 1% oil shock:\n"); print(gas_key_values)
cat("\nDiesel response to a 1% oil shock:\n");   print(diesel_key_values)

# ------------------------------------------------------------------------
# 3. Plot IRFs
# ------------------------------------------------------------------------
plot(irf_gas$irf,    main = "Response of Gasoline Prices to a Brent Oil Shock")
plot(irf_diesel$irf, main = "Response of Diesel Prices to a Brent Oil Shock")

# ------------------------------------------------------------------------
# 4. Summary metrics (peak response, cumulative effect, weeks significant)
# ------------------------------------------------------------------------
summarize_irf <- function(irf_obj, weeks_full = 0:20) {
  vals <- extract_irf_values(irf_obj, weeks_full)
  cum_response <- cumsum(vals$response)
  list(
    peak_response     = max(vals$response),
    week_of_peak      = vals$week[which.max(vals$response)],
    immediate         = vals$response[vals$week == 1],
    short_run_w4      = vals$response[vals$week == 4],
    medium_run_w12    = vals$response[vals$week == 12],
    cumulative_w20    = cum_response[length(cum_response)],
    weeks_significant = sum(vals$significant == "Yes")
  )
}

irf_summary_gas    <- summarize_irf(irf_gas$irf)
irf_summary_diesel <- summarize_irf(irf_diesel$irf)

cat("\n=== IRF SUMMARY: GASOLINE ===\n"); str(irf_summary_gas)
cat("\n=== IRF SUMMARY: DIESEL ===\n");   str(irf_summary_diesel)

# ------------------------------------------------------------------------
# 5. Export summary table to LaTeX
# ------------------------------------------------------------------------
irf_summary_table <- data.frame(
  Metric = c("Peak response", "Week of peak", "Immediate (week 1)",
             "Short-run (week 4)", "Medium-run (week 12)",
             "Cumulative (week 20)", "Weeks significant"),
  Gasoline = round(unlist(irf_summary_gas[c("peak_response", "week_of_peak", "immediate",
                                             "short_run_w4", "medium_run_w12",
                                             "cumulative_w20", "weeks_significant")]), 4),
  Diesel   = round(unlist(irf_summary_diesel[c("peak_response", "week_of_peak", "immediate",
                                                "short_run_w4", "medium_run_w12",
                                                "cumulative_w20", "weeks_significant")]), 4)
)

print(xtable(irf_summary_table,
             caption = "Summary Statistics of Impulse Response Functions",
             label = "tab:irf_summary"),
      include.rownames = FALSE, caption.placement = "top", booktabs = TRUE,
      file = file.path("output", "table_irf_summary.tex"))

cat("\nTable exported to output/table_irf_summary.tex\n")
