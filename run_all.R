# ==============================================================================
# 01 - DESCRIPTIVE STATISTICS & ACF/PACF DIAGNOSTICS
# ==============================================================================
# source("R/00_setup_data.R") must be run first.
# ==============================================================================

library(moments)   # skewness, kurtosis
library(gridExtra)
library(RColorBrewer)

# ------------------------------------------------------------------------
# 1. Descriptive statistics table (levels, logs, first differences)
# ------------------------------------------------------------------------
calc_descriptive_stats <- function(ts_data, series_name) {
  clean_data <- na.omit(as.numeric(ts_data))
  jb_test <- jarque.bera.test(clean_data)

  data.frame(
    Series    = series_name,
    N_Obs     = length(clean_data),
    Mean      = round(mean(clean_data), 4),
    Median    = round(median(clean_data), 4),
    Std_Dev   = round(sd(clean_data), 4),
    Min       = round(min(clean_data), 4),
    Max       = round(max(clean_data), 4),
    Skewness  = round(skewness(clean_data), 4),
    Kurtosis  = round(kurtosis(clean_data), 4),
    CV        = round(sd(clean_data) / mean(clean_data), 4),
    JB_Stat   = round(as.numeric(jb_test$statistic), 4),
    JB_pvalue = round(jb_test$p.value, 4),
    Normal_95 = ifelse(jb_test$p.value > 0.05, "Yes", "No")
  )
}

series_list <- list(
  "Gasoline (levels)"   = ts_gasoline,
  "Diesel (levels)"     = ts_diesel,
  "Brent 3M (levels)"   = ts_brent_3M,
  "Brent SPOT (levels)" = ts_brent_SPOT,
  "Log Gasoline"        = log_gasoline,
  "Log Diesel"          = log_diesel,
  "Log Brent 3M"        = log_brent_3M,
  "Log Brent SPOT"      = log_brent_SPOT,
  "D Log Gasoline"      = d_log_gasoline,
  "D Log Diesel"        = d_log_diesel,
  "D Log Brent 3M"      = d_log_brent_3M
)

descriptive_table <- do.call(rbind, Map(calc_descriptive_stats, series_list, names(series_list)))
print(descriptive_table)

knitr::kable(descriptive_table,
             caption = "Descriptive Statistics for Fuel Price Series (2005-2024)",
             digits = 4)

# ------------------------------------------------------------------------
# 2. ACF / PACF plots
# ------------------------------------------------------------------------
create_acf_pacf_plots <- function(series, series_name, max_lags = 24) {
  acf_result  <- acf(series, lag.max = max_lags, plot = FALSE)
  pacf_result <- pacf(series, lag.max = max_lags, plot = FALSE)
  ci_line <- qnorm((1 + 0.95) / 2) / sqrt(length(series))

  acf_data  <- data.frame(Lag = 1:max_lags, ACF = as.numeric(acf_result$acf[-1]),
                           Upper_CI = ci_line, Lower_CI = -ci_line)
  pacf_data <- data.frame(Lag = 1:max_lags, PACF = as.numeric(pacf_result$acf),
                           Upper_CI = ci_line, Lower_CI = -ci_line)

  p1 <- ggplot(acf_data, aes(x = Lag, y = ACF)) +
    geom_hline(yintercept = 0) +
    geom_hline(aes(yintercept = Upper_CI), linetype = "dashed", color = "blue") +
    geom_hline(aes(yintercept = Lower_CI), linetype = "dashed", color = "blue") +
    geom_segment(aes(xend = Lag, yend = 0)) +
    labs(title = paste("ACF -", series_name)) +
    theme_minimal()

  p2 <- ggplot(pacf_data, aes(x = Lag, y = PACF)) +
    geom_hline(yintercept = 0) +
    geom_hline(aes(yintercept = Upper_CI), linetype = "dashed", color = "blue") +
    geom_hline(aes(yintercept = Lower_CI), linetype = "dashed", color = "blue") +
    geom_segment(aes(xend = Lag, yend = 0)) +
    labs(title = paste("PACF -", series_name)) +
    theme_minimal()

  grid.arrange(p1, p2, ncol = 2)
}

create_acf_pacf_plots(log_gasoline, "Log Gasoline")
create_acf_pacf_plots(log_diesel,   "Log Diesel")
create_acf_pacf_plots(log_brent_3M, "Log Brent 3M")
