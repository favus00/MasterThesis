# ==============================================================================
# 04 - ASYMMETRIC PASS-THROUGH: NARDL & ASYMMETRIC ERROR CORRECTION MODEL
# ==============================================================================
# source("R/00_setup_data.R") and R/02_stationarity_cointegration.R must be
# run first.
#
# Both models test whether retail fuel prices respond differently to oil
# price increases than to decreases ("rockets and feathers" pass-through).
#   - NARDL (Shin, Yu & Greenwood-Nimmo, 2014): decomposes the oil price
#     level itself into cumulative positive/negative partial sums.
#   - Asymmetric ECM (Granger & Lee, 1989): decomposes both the oil price
#     change AND the cointegration error term into positive/negative parts.
# ==============================================================================

# ------------------------------------------------------------------------
# 1. Shared helper: split a series into positive/negative components
# ------------------------------------------------------------------------
split_pos_neg <- function(x) {
  list(pos = ifelse(x > 0, x, 0), neg = ifelse(x < 0, x, 0))
}

# ==============================================================================
# PART A - NARDL
# ==============================================================================
run_nardl <- function(log_fuel, log_oil, fuel_name, p = 4, q = 4, h_max = 12) {

  delta_oil <- diff(log_oil)
  oil_split <- split_pos_neg(delta_oil)
  oil_pos <- cumsum(oil_split$pos)   # cumulative partial sum of oil increases
  oil_neg <- cumsum(oil_split$neg)   # cumulative partial sum of oil decreases

  df <- data.frame(y = log_fuel[-1], oil_pos = oil_pos, oil_neg = oil_neg)

  max_lag <- max(p, q)
  model_data <- data.frame(
    y = df$y,
    create_lags(df$y, p),
    create_lags(df$oil_pos, q),
    create_lags(df$oil_neg, q)
  )
  colnames(model_data) <- c("y",
                             paste0("y_lag", 1:p),
                             paste0("oil_pos_lag", 1:q),
                             paste0("oil_neg_lag", 1:q))
  model_data <- model_data[-(1:max_lag), ]

  formula <- as.formula(paste("y ~", paste(colnames(model_data)[-1], collapse = " + ")))
  model <- lm(formula, data = model_data)

  coefs    <- coef(model)
  phi      <- coefs[grep("y_lag", names(coefs))]
  beta_pos <- coefs[grep("oil_pos_lag", names(coefs))]
  beta_neg <- coefs[grep("oil_neg_lag", names(coefs))]

  # Long-run multipliers: sum of short-run coefficients / (1 - sum of AR coefficients)
  lrm_pos <- sum(beta_pos) / (1 - sum(phi))
  lrm_neg <- sum(beta_neg) / (1 - sum(phi))

  cat(sprintf("\n=== %s NARDL (p=%d, q=%d) ===\n", toupper(fuel_name), p, q))
  print(summary(model))
  cat(sprintf("Long-run multiplier, oil increases: %.4f\n", lrm_pos))
  cat(sprintf("Long-run multiplier, oil decreases: %.4f\n", lrm_neg))

  # Cumulative dynamic multipliers (how the effect builds up over h weeks)
  cdm_pos <- cumsum(beta_pos[1:min(h_max, length(beta_pos))])
  cdm_neg <- cumsum(beta_neg[1:min(h_max, length(beta_neg))])

  list(model = model, lrm_pos = lrm_pos, lrm_neg = lrm_neg,
       cdm_pos = cdm_pos, cdm_neg = cdm_neg)
}

nardl_gasoline <- run_nardl(log_gasoline, log_brent_3M, "Gasoline")
nardl_diesel   <- run_nardl(log_diesel,   log_brent_3M, "Diesel")

plot_cumulative_multipliers <- function(nardl_result, fuel_name) {
  h_max <- length(nardl_result$cdm_pos)
  df_long <- pivot_longer(
    data.frame(h = 1:h_max, Positive = nardl_result$cdm_pos, Negative = nardl_result$cdm_neg),
    cols = c("Positive", "Negative"), names_to = "Shock", values_to = "Multiplier"
  )
  ggplot(df_long, aes(x = h, y = Multiplier, color = Shock)) +
    geom_line(linewidth = 1.2) + geom_point(size = 2) +
    labs(title = paste("Cumulative Dynamic Multipliers -", fuel_name),
         x = "Weeks after shock", y = "Multiplier") +
    theme_minimal() +
    scale_color_manual(values = c(Positive = "steelblue", Negative = "firebrick"))
}

print(plot_cumulative_multipliers(nardl_gasoline, "Gasoline"))
print(plot_cumulative_multipliers(nardl_diesel,   "Diesel"))

# ==============================================================================
# PART B - ASYMMETRIC ERROR CORRECTION MODEL (Granger & Lee, 1989)
# ==============================================================================
run_asymmetric_ecm <- function(d_log_fuel, d_log_oil, ect, fuel_name) {

  oil_split <- split_pos_neg(c(NA, d_log_oil))
  oil_pos      <- oil_split$pos[-1]
  oil_neg      <- oil_split$neg[-1]
  oil_pos_lag1 <- c(NA, oil_split$pos[-length(oil_split$pos)])[-1]
  oil_neg_lag1 <- c(NA, oil_split$neg[-length(oil_split$neg)])[-1]

  ect_lag <- ect[1:(length(ect) - 1)]
  ect_split <- split_pos_neg(ect_lag)

  model <- lm(d_log_fuel ~ oil_pos + oil_neg + oil_pos_lag1 + oil_neg_lag1 +
                ect_split$pos + ect_split$neg)
  robust <- coeftest(model, vcov = vcovHAC(model))

  cat(sprintf("\n=== %s ASYMMETRIC ECM ===\n", toupper(fuel_name)))
  print(robust)

  coefs <- coef(model)
  cat("Short-run response to oil increases (beta+):", round(coefs["oil_pos"], 4), "\n")
  cat("Short-run response to oil decreases (beta-):", round(coefs["oil_neg"], 4), "\n")
  cat("Adjustment when above equilibrium (lambda+):", round(coefs[5], 4), "\n")
  cat("Adjustment when below equilibrium (lambda-):", round(coefs[6], 4), "\n")
  cat("(lambda+ more negative than lambda- => faster pass-through of increases,",
      "consistent with 'rockets and feathers')\n")

  list(model = model, robust = robust)
}

aecm_gasoline <- run_asymmetric_ecm(d_log_gasoline, d_log_brent_3M, ect_gasoline, "Gasoline")
aecm_diesel   <- run_asymmetric_ecm(d_log_diesel,   d_log_brent_3M, ect_diesel,   "Diesel")
