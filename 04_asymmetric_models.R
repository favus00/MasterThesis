# ==============================================================================
# 06 - GRANGER CAUSALITY
# ==============================================================================
# source("R/00_setup_data.R") must be run first.
#
# Tests the direction of price transmission: does oil Granger-cause retail
# fuel prices (and not vice-versa, as economic theory would predict), and is
# there any cross-fuel dynamic between gasoline and diesel?
# ==============================================================================

# ------------------------------------------------------------------------
# 1. Reusable Granger causality test (unrestricted vs. restricted model,
#    HAC-robust Wald test) - replaces six separate near-identical blocks
# ------------------------------------------------------------------------
granger_test_hac <- function(y, x, lags = 5, desc = "X -> Y") {

  n <- length(y)
  data <- data.frame(y = y, x = x)
  for (i in 1:lags) {
    data[[paste0("y_lag", i)]] <- c(rep(NA, i), y[1:(n - i)])
    data[[paste0("x_lag", i)]] <- c(rep(NA, i), x[1:(n - i)])
  }
  data <- na.omit(data)

  f_unres <- as.formula(paste("y ~", paste0("y_lag", 1:lags, collapse = " + "),
                               "+", paste0("x_lag", 1:lags, collapse = " + ")))
  f_res   <- as.formula(paste("y ~", paste0("y_lag", 1:lags, collapse = " + ")))

  model_unres <- lm(f_unres, data = data)
  model_res   <- lm(f_res,   data = data)

  wald <- waldtest(model_unres, model_res, vcov = vcovHAC(model_unres))
  f_stat  <- wald$F[2]
  p_value <- wald$`Pr(>F)`[2]

  cat(sprintf("%-24s F = %7.3f,  p = %s\n", desc, f_stat,
              ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))))

  list(f_stat = f_stat, p_value = p_value, desc = desc, n_obs = nrow(data))
}

# ------------------------------------------------------------------------
# 2. Run all six pairwise tests (5 lags, matching the VAR lag selection)
# ------------------------------------------------------------------------
cat("Granger causality tests (5 lags, HAC-robust):\n\n")

gc_oil_gas     <- granger_test_hac(d_log_gasoline, d_log_brent_3M, 5, "Oil -> Gasoline")
gc_gas_oil     <- granger_test_hac(d_log_brent_3M, d_log_gasoline, 5, "Gasoline -> Oil")
gc_oil_diesel  <- granger_test_hac(d_log_diesel,   d_log_brent_3M, 5, "Oil -> Diesel")
gc_diesel_oil  <- granger_test_hac(d_log_brent_3M, d_log_diesel,   5, "Diesel -> Oil")
gc_gas_diesel  <- granger_test_hac(d_log_diesel,   d_log_gasoline, 5, "Gasoline -> Diesel")
gc_diesel_gas  <- granger_test_hac(d_log_gasoline, d_log_diesel,   5, "Diesel -> Gasoline")

# ------------------------------------------------------------------------
# 3. Summary table + interpretation
# ------------------------------------------------------------------------
granger_results <- data.frame(
  Direction = c(gc_oil_gas$desc, gc_gas_oil$desc, gc_oil_diesel$desc,
                gc_diesel_oil$desc, gc_gas_diesel$desc, gc_diesel_gas$desc),
  F_stat  = round(c(gc_oil_gas$f_stat, gc_gas_oil$f_stat, gc_oil_diesel$f_stat,
                     gc_diesel_oil$f_stat, gc_gas_diesel$f_stat, gc_diesel_gas$f_stat), 3),
  p_value = round(c(gc_oil_gas$p_value, gc_gas_oil$p_value, gc_oil_diesel$p_value,
                     gc_diesel_oil$p_value, gc_gas_diesel$p_value, gc_diesel_gas$p_value), 4)
)
print(granger_results)

cat("\nInterpretation:\n")
cat("- Oil -> Gasoline / Oil -> Diesel: expected to be strongly significant,",
    "consistent with oil driving retail prices.\n")
cat("- Reverse directions (fuel -> oil): a small country's retail fuel demand",
    "should not move world oil prices; significance here would point to a",
    "spurious/confounded relationship rather than genuine causality.\n")
cat("- Gasoline <-> Diesel: tests whether the two products move together",
    "beyond their shared dependence on oil.\n")

# ------------------------------------------------------------------------
# 4. Export LaTeX table
# ------------------------------------------------------------------------
format_pval <- function(p) if (p < 0.001) "$<$0.001" else sprintf("%.3f", p)
get_sig <- function(p) if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else ""

granger_results$Sig <- sapply(granger_results$p_value, get_sig)
print(xtable(granger_results,
             caption = "Granger Causality Test Results: Direction of Price Transmission",
             label = "tab:granger_causality"),
      include.rownames = FALSE, caption.placement = "top", booktabs = TRUE,
      file = file.path("output", "table_granger_causality.tex"))

cat("\nTable exported to output/table_granger_causality.tex\n")
