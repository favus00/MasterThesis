# =============================================================================
# Asymmetric price transmission in the Italian fuel market
# Single end-to-end script (data -> descriptive stats -> unit roots ->
# cointegration -> VECM -> Granger -> IRF -> NARDL)
#
# Order follows the thesis chapters:
#   1  Data                                   (Ch. 3.1-3.2, 3.4.1)
#   2  Descriptive statistics and figures     (Ch. 3.3)
#   3  Unit root tests                        (Ch. 3.4.2)
#   4  VAR lag selection                      (Ch. 3.5.1)
#   5  Engle-Granger cointegration            (Ch. 3.5.2)
#   6  VECM (5 lags main, 3 lags robustness)  (Ch. 3.6)
#   7  VECM residual diagnostics              (Ch. 3.6)
#   8  Granger causality                      (Ch. 3.6.1)
#   9  Impulse response functions             (Ch. 3.6.2)
#  10  Forward vs spot oil prices             (Ch. 3.3.4, optional)
#  11  NARDL(6,6,6)                           (Ch. 4)
#  12  Check against the values reported in the thesis (optional)
#
# Usage: set the working directory to the folder containing the data file,
# then source("thesis_pipeline.R"). Tables go to output/tables, figures to
# output/figures.
# =============================================================================


# ----------------------------------------------------------------------------
# 0. SETUP
# ----------------------------------------------------------------------------


packages <- c("readxl", "urca", "tseries", "vars", "lmtest", "sandwich",
              "car", "dynlm", "xtable", "ggplot2", "moments")
to_install <- setdiff(packages, rownames(installed.packages()))
if (length(to_install) > 0) {
  stop("Please install the missing packages first: ",
       paste(to_install, collapse = ", "))
}
invisible(lapply(packages, library, character.only = TRUE))



setwd("C:/Users/ricky/Desktop/GitHub thesis")
list.files()   # devi vedere sia thesis_pipeline.R sia Dati_completo.xlsx

# --- Configuration -----------------------------------------------------------

DATA_FILE          <- "Dati_completo.xlsx"
OUT_DIR            <- "output"
VECM_LAGS_MAIN     <- 5     # main VECM specification (also used for Granger/IRF)
VECM_LAGS_ROBUST   <- 3     # robustness check (lag-selection benchmark)
NARDL_LAGS         <- 6     # NARDL(6,6,6)
IRF_HORIZON        <- 20    # weeks
IRF_RUNS           <- 500   # bootstrap replications
DM_HORIZON         <- 50    # weeks, dynamic multipliers
DM_SHOCK           <- 0.01  # size of the permanent shock (1%)
RUN_FORWARD_VS_SPOT      <- TRUE   # section 10
RUN_NARDL_LAG_SELECTION  <- TRUE   # AIC/BIC over NARDL(p,p,p), section 11
RUN_THESIS_CHECK         <- TRUE   # section 12

set.seed(2026)  # bootstrap confidence bands are reproducible

dir.create(file.path(OUT_DIR, "tables"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)

# --- Helpers -----------------------------------------------------------------
section <- function(txt) {
  cat("\n", strrep("=", 72), "\n", txt, "\n", strrep("=", 72), "\n", sep = "")
}

save_table <- function(df, name, caption = name, digits = 4) {
  write.csv(df, file.path(OUT_DIR, "tables", paste0(name, ".csv")),
            row.names = FALSE)
  print(xtable(df, caption = caption, digits = digits),
        file = file.path(OUT_DIR, "tables", paste0(name, ".tex")),
        include.rownames = FALSE, booktabs = TRUE,
        caption.placement = "top")
  invisible(df)
}

save_plot <- function(p, name, width = 10, height = 6) {
  ggsave(file.path(OUT_DIR, "figures", paste0(name, ".png")), plot = p,
         width = width, height = height, dpi = 300)
}

stars <- function(p) {
  ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
}

fuel_names <- c(gasoline = "Gasoline", diesel = "Diesel")
fuel_colors <- c("Gasoline" = "#E74C3C", "Diesel" = "#3498DB")


# ----------------------------------------------------------------------------
# 1. DATA
# ----------------------------------------------------------------------------
# Sources: retail pre-tax prices (EUR per 1,000 litres, EC Weekly Oil Bulletin),
# Brent 3rd continuous futures and Brent spot (USD/bbl) and exchange rates
# (LSEG). Brent is converted to EUR/bbl with the forward rate of matching
# maturity for futures and with the spot rate for spot prices.

section("1. DATA")

raw <- read_excel(DATA_FILE)

COL <- list(
  date   = "Date",
  gas    = "Gasoline",
  dsl    = "Diesel",
  fut    = "ICE - Brent  Continuous 3rd future Settlment Price",  # double space is in the source file
  fx3m   = "EXCHANGE RATE 3M FWD",
  spot   = "Brent Spot USD/BBL",
  fxspot = "USD TO EUR (LSEG DS) - EXCHANGE RATE"
)
missing_cols <- setdiff(unlist(COL), names(raw))
if (length(missing_cols) > 0) {
  stop("Columns not found in ", DATA_FILE, ": ", paste(missing_cols, collapse = " | "))
}

dat_all <- data.frame(
  Date       = as.Date(raw[[COL$date]]),
  Gasoline   = as.numeric(raw[[COL$gas]]),
  Diesel     = as.numeric(raw[[COL$dsl]]),
  Brent_3M   = as.numeric(raw[[COL$fut]])  / as.numeric(raw[[COL$fx3m]]),
  Brent_SPOT = as.numeric(raw[[COL$spot]]) / as.numeric(raw[[COL$fxspot]])
)

# Missing weeks are dropped (listwise deletion), not interpolated
dat <- na.omit(dat_all)
cat("Observations in file:      ", nrow(dat_all), "\n")
cat("Observations after na.omit:", nrow(dat), "\n")

make_ts <- function(x) ts(x, start = c(2005, 1), frequency = 52)

ts_gasoline   <- make_ts(dat$Gasoline)
ts_diesel     <- make_ts(dat$Diesel)
ts_brent_3M   <- make_ts(dat$Brent_3M)
ts_brent_spot <- make_ts(dat$Brent_SPOT)

log_gasoline   <- log(ts_gasoline)
log_diesel     <- log(ts_diesel)
log_brent_3M   <- log(ts_brent_3M)
log_brent_spot <- log(ts_brent_spot)

d_log_gasoline   <- diff(log_gasoline)
d_log_diesel     <- diff(log_diesel)
d_log_brent_3M   <- diff(log_brent_3M)
d_log_brent_spot <- diff(log_brent_spot)

dates <- dat$Date
log_fuel <- list(gasoline = log_gasoline, diesel = log_diesel)
d_log_fuel <- list(gasoline = d_log_gasoline, diesel = d_log_diesel)


# ----------------------------------------------------------------------------
# 2. DESCRIPTIVE STATISTICS AND FIGURES
# ----------------------------------------------------------------------------

section("2. DESCRIPTIVE STATISTICS AND FIGURES")

desc_stats <- function(x, label) {
  x <- na.omit(as.numeric(x))
  jb <- jarque.bera.test(x)
  data.frame(
    Series    = label,
    N_Obs     = length(x),
    Mean      = mean(x),
    Median    = median(x),
    Std_Dev   = sd(x),
    Min       = min(x),
    Max       = max(x),
    Skewness  = skewness(x),
    Kurtosis  = kurtosis(x),
    CV        = sd(x) / mean(x),
    JB_Stat   = as.numeric(jb$statistic),
    JB_pvalue = jb$p.value,
    Normal_95 = ifelse(jb$p.value > 0.05, "Yes", "No"),
    stringsAsFactors = FALSE
  )
}

base_series <- list("Gasoline" = ts_gasoline, "Diesel" = ts_diesel,
                    "Brent 3M" = ts_brent_3M, "Brent SPOT" = ts_brent_spot)
transforms <- list(
  list(prefix = "",         suffix = " (levels)", f = identity),
  list(prefix = "Log ",     suffix = "",          f = log),
  list(prefix = "D ",       suffix = "",          f = diff),
  list(prefix = "D Log ",   suffix = "",          f = function(x) diff(log(x)))
)
desc_rows <- list()
for (tr in transforms) {
  for (nm in names(base_series)) {
    desc_rows[[length(desc_rows) + 1]] <-
      desc_stats(tr$f(base_series[[nm]]), paste0(tr$prefix, nm, tr$suffix))
  }
}
descriptive_table <- do.call(rbind, desc_rows)
print(descriptive_table)
save_table(descriptive_table, "table_descriptive_statistics",
           "Descriptive statistics (2005-2024)", digits = 4)

# --- Figures -----------------------------------------------------------------
events <- data.frame(
  label = c("2008 crisis", "COVID-19", "Ukraine war"),
  start = as.Date(c("2008-06-01", "2020-03-01", "2022-02-24")),
  end   = as.Date(c("2009-06-30", "2020-10-31", "2024-01-01")),
  fill  = c("red", "orange", "purple"),
  stringsAsFactors = FALSE
)
add_events <- function(p) {
  for (i in seq_len(nrow(events))) {
    p <- p + annotate("rect", xmin = events$start[i], xmax = events$end[i],
                      ymin = -Inf, ymax = Inf, alpha = 0.1, fill = events$fill[i])
  }
  p
}

df_levels <- data.frame(Date = dates, Gasoline = dat$Gasoline, Diesel = dat$Diesel,
                        Brent_3M = dat$Brent_3M, Brent_SPOT = dat$Brent_SPOT)

p_fuel <- add_events(
  ggplot(df_levels, aes(x = Date)) +
    geom_line(aes(y = Gasoline, color = "Gasoline"), linewidth = 0.7) +
    geom_line(aes(y = Diesel, color = "Diesel"), linewidth = 0.7) +
    scale_color_manual(values = fuel_colors) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    labs(title = "Italian fuel prices (2005-2024)", subtitle = "Pre-tax retail prices",
         x = "Year", y = "EUR per 1,000 litres", color = "Fuel") +
    theme_minimal() + theme(legend.position = "bottom",
                            axis.text.x = element_text(angle = 45, hjust = 1))
)
save_plot(p_fuel, "fig_fuel_prices_levels")

p_brent <- add_events(
  ggplot(df_levels, aes(x = Date)) +
    geom_line(aes(y = Brent_3M, color = "3M futures"), linewidth = 0.7) +
    geom_line(aes(y = Brent_SPOT, color = "Spot"), linewidth = 0.7, alpha = 0.8) +
    scale_color_manual(values = c("3M futures" = "#2ECC71", "Spot" = "#F39C12")) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    labs(title = "Brent: spot vs 3-month futures (2005-2024)",
         subtitle = "Converted to EUR per barrel", x = "Year",
         y = "EUR per barrel", color = "Price") +
    theme_minimal() + theme(legend.position = "bottom",
                            axis.text.x = element_text(angle = 45, hjust = 1))
)
save_plot(p_brent, "fig_brent_spot_vs_futures")

df_chg <- data.frame(
  Date   = rep(dates[-1], 3),
  value  = c(as.numeric(d_log_gasoline), as.numeric(d_log_diesel),
             as.numeric(d_log_brent_3M)),
  series = rep(c("Gasoline", "Diesel", "Brent (3M futures)"),
               each = length(d_log_gasoline))
)
p_chg <- add_events(
  ggplot(df_chg, aes(x = Date, y = value)) +
    geom_line(linewidth = 0.3) +
    facet_wrap(~ series, ncol = 1, scales = "free_y") +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    labs(title = "Weekly log-differences", x = "Year", y = "Log-difference") +
    theme_minimal()
)
save_plot(p_chg, "fig_log_differences", height = 8)

# ACF / PACF of log-levels and of first differences
acf_pacf_png <- function(series_list, file, title) {
  png(file.path(OUT_DIR, "figures", file), width = 1200, height = 1000, res = 150)
  par(mfrow = c(length(series_list), 2), mar = c(4, 4, 3, 1))
  for (nm in names(series_list)) {
    acf(as.numeric(series_list[[nm]]),  lag.max = 24, main = paste(nm, "- ACF"))
    pacf(as.numeric(series_list[[nm]]), lag.max = 24, main = paste(nm, "- PACF"))
  }
  dev.off()
}
acf_pacf_png(list(Gasoline = log_gasoline, Diesel = log_diesel, "Brent 3M" = log_brent_3M),
             "fig_acf_pacf_log_levels.png", "log levels")
acf_pacf_png(list(Gasoline = d_log_gasoline, Diesel = d_log_diesel, "Brent 3M" = d_log_brent_3M),
             "fig_acf_pacf_first_differences.png", "first differences")


# ----------------------------------------------------------------------------
# 3. UNIT ROOT TESTS (constant-only specifications)
# ----------------------------------------------------------------------------
# ADF (drift, AIC), Phillips-Perron (Z-tau, constant), KPSS (level) and
# Zivot-Andrews (break in intercept). All figures in the table are computed
# here; nothing is typed by hand.

section("3. UNIT ROOT TESTS")

unit_root_objects <- function(x) {
  list(
    adf = ur.df(x, type = "drift", selectlags = "AIC"),
    pp  = ur.pp(x, type = "Z-tau", model = "constant"),
    kpss = ur.kpss(x, type = "mu"),
    za  = ur.za(x, model = "intercept", lag = NULL)
  )
}

ur_series <- list("Gasoline price"  = log_gasoline,
                  "Diesel price"    = log_diesel,
                  "Brent 3M future" = log_brent_3M)
ur_objs <- lapply(ur_series, unit_root_objects)

unit_root_table <- data.frame(
  Variable = names(ur_series),
  ADF  = sapply(ur_objs, function(o) o$adf@teststat[1]),
  PP   = sapply(ur_objs, function(o) as.numeric(o$pp@teststat)),
  KPSS = sapply(ur_objs, function(o) as.numeric(o$kpss@teststat)),
  ZA   = sapply(ur_objs, function(o) as.numeric(o$za@teststat)),
  row.names = NULL
)

# Critical values are read from the test objects (they do not depend on the
# series). Note the KPSS columns are 10%, 5%, 2.5%, 1%: select them by name.
o1 <- ur_objs[[1]]
cv_rows <- data.frame(
  Variable = c("Critical value 1%", "Critical value 5%", "Critical value 10%"),
  ADF  = o1$adf@cval["tau2", c("1pct", "5pct", "10pct")],
  PP   = o1$pp@cval[1, c("1pct", "5pct", "10pct")],
  KPSS = o1$kpss@cval[1, c("1pct", "5pct", "10pct")],
  ZA  = unname(o1$za@cval[1:3]),
  row.names = NULL
)
unit_root_table <- rbind(unit_root_table, cv_rows)
print(unit_root_table)
save_table(unit_root_table, "table_unit_root_tests",
           "Unit root tests (H0: unit root for ADF, PP, ZA; H0: stationarity for KPSS)",
           digits = 3)

# First differences should be stationary (series are I(1))
ur_diff <- data.frame(
  Series = c("D Log Gasoline", "D Log Diesel", "D Log Brent 3M"),
  ADF_stat = sapply(list(d_log_gasoline, d_log_diesel, d_log_brent_3M),
                    function(x) ur.df(x, type = "drift", selectlags = "AIC")@teststat[1]),
  Critical_5pct = ur_objs[[1]]$adf@cval["tau2", "5pct"]
)
ur_diff$Stationary_5pct <- ur_diff$ADF_stat < ur_diff$Critical_5pct
print(ur_diff)
save_table(ur_diff, "table_unit_root_first_differences",
           "ADF test on first differences", digits = 3)


# ----------------------------------------------------------------------------
# 4. VAR LAG SELECTION (on first differences)
# ----------------------------------------------------------------------------

section("4. LAG SELECTION")

lag_selection <- function(d_fuel) {
  sel <- VARselect(cbind(d_brent = d_log_brent_3M, d_fuel = d_fuel),
                   lag.max = 12, type = "const")
  data.frame(Lag = 1:12,
             AIC = as.numeric(sel$criteria["AIC(n)", ]),
             BIC = as.numeric(sel$criteria["SC(n)", ]))
}
lag_gas <- lag_selection(d_log_gasoline)
lag_dsl <- lag_selection(d_log_diesel)
lag_table <- data.frame(Lag = 1:12,
                        AIC_Gasoline = lag_gas$AIC, BIC_Gasoline = lag_gas$BIC,
                        AIC_Diesel   = lag_dsl$AIC, BIC_Diesel   = lag_dsl$BIC)
print(lag_table)
cat("Lag minimising AIC / BIC - Gasoline:", which.min(lag_gas$AIC), "/", which.min(lag_gas$BIC), "\n")
cat("Lag minimising AIC / BIC - Diesel:  ", which.min(lag_dsl$AIC), "/", which.min(lag_dsl$BIC), "\n")
save_table(lag_table, "table_var_lag_selection",
           "VAR lag selection: AIC and BIC", digits = 4)


# ----------------------------------------------------------------------------
# 5. ENGLE-GRANGER COINTEGRATION
# ----------------------------------------------------------------------------
# Step 1: static regression log(fuel) = a + b * log(Brent 3M) + u
# Step 2: ADF on the residuals (drift, 3 lags). Engle-Granger critical values.

section("5. ENGLE-GRANGER COINTEGRATION")

EG_CV <- c("1%" = -3.90, "5%" = -3.34, "10%" = -3.04)

engle_granger <- function(log_y, log_x) {
  d <- data.frame(y = as.numeric(log_y), x = as.numeric(log_x))
  lr <- lm(y ~ x, data = d)
  adf <- ur.df(residuals(lr), type = "drift", lags = 3, selectlags = "Fixed")
  list(model = lr, ect = residuals(lr), adf = adf,
       beta = unname(coef(lr)["x"]), r2 = summary(lr)$r.squared)
}
eg <- lapply(log_fuel, engle_granger, log_x = log_brent_3M)

for (f in names(eg)) {
  cat("\n---", fuel_names[f], ": long-run regression (HAC standard errors) ---\n")
  print(coeftest(eg[[f]]$model, vcov = vcovHAC(eg[[f]]$model)))
}

eg_table <- data.frame(
  Relationship = paste("Oil ->", fuel_names),
  Beta         = sapply(eg, function(e) e$beta),
  R_squared    = sapply(eg, function(e) e$r2),
  ADF_stat     = sapply(eg, function(e) e$adf@teststat[1]),
  Cointegrated_5pct = ifelse(sapply(eg, function(e) e$adf@teststat[1]) < EG_CV["5%"],
                             "Yes", "No"),
  row.names = NULL
)
print(eg_table)
cat("Engle-Granger critical values: 1% = -3.90, 5% = -3.34, 10% = -3.04\n")
save_table(eg_table, "table_engle_granger", "Engle-Granger cointegration test", digits = 4)

# Residuals (equilibrium errors)
png(file.path(OUT_DIR, "figures", "fig_cointegration_residuals.png"),
    width = 1200, height = 900, res = 150)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (f in names(eg)) {
  plot(eg[[f]]$ect, type = "l", main = paste(fuel_names[f], "- cointegration residuals"),
       ylab = "Residual", xlab = "Time")
  abline(h = 0, col = "red", lty = 2)
  acf(eg[[f]]$ect, main = paste(fuel_names[f], "- ACF of residuals"))
}
dev.off()


# ----------------------------------------------------------------------------
# 6. VECM (single-equation error-correction form, estimated by OLS)
# ----------------------------------------------------------------------------
#   d(fuel_t) = c + g0*d(brent_t) + sum_i g_i*d(brent_{t-i})
#               + sum_i h_i*d(fuel_{t-i}) + alpha*ECT_{t-1} + e_t
# ECT = residuals of the Engle-Granger regression of section 5.

section("6. VECM")

make_lags <- function(x, k) {
  n <- length(x)
  m <- sapply(seq_len(k), function(i) c(rep(NA, i), x[seq_len(n - i)]))
  matrix(m, nrow = n)
}

estimate_vecm <- function(d_fuel, d_brent, ect, lags) {
  d_fuel <- as.numeric(d_fuel); d_brent <- as.numeric(d_brent)
  n <- length(d_fuel)
  B <- make_lags(d_brent, lags); colnames(B) <- paste0("d_brent_L", seq_len(lags))
  H <- make_lags(d_fuel, lags);  colnames(H) <- paste0("d_fuel_L",  seq_len(lags))
  df <- data.frame(d_fuel = d_fuel, d_brent = d_brent, B, H,
                   ect_lag = ect[seq_len(n)])   # ECT(t-1) aligned with the differences
  df <- na.omit(df)
  rhs <- c("d_brent", colnames(B), colnames(H), "ect_lag")
  fit <- lm(reformulate(rhs, response = "d_fuel"), data = df)
  alpha <- unname(coef(fit)["ect_lag"])
  list(model = fit, data = df, lags = lags, n = nrow(df),
       robust = coeftest(fit, vcov = vcovHAC(fit)),
       alpha = alpha,
       speed_pct = abs(alpha) * 100,
       half_life = log(0.5) / log(1 + alpha),
       sr_elast = unname(coef(fit)["d_brent"]),
       r2 = summary(fit)$r.squared, adj_r2 = summary(fit)$adj.r.squared)
}

vecm <- list()
for (k in c(VECM_LAGS_ROBUST, VECM_LAGS_MAIN)) {
  for (f in names(fuel_names)) {
    vecm[[paste0(f, "_", k)]] <-
      estimate_vecm(d_log_fuel[[f]], d_log_brent_3M, eg[[f]]$ect, lags = k)
  }
}

for (f in names(fuel_names)) {
  cat("\n--- VECM", fuel_names[f], "(", VECM_LAGS_MAIN, "lags, HAC standard errors ) ---\n")
  print(vecm[[paste0(f, "_", VECM_LAGS_MAIN)]]$robust)
}

vecm_summary <- do.call(rbind, lapply(names(vecm), function(key) {
  v <- vecm[[key]]
  f <- sub("_[0-9]+$", "", key)
  data.frame(Fuel = fuel_names[f], Lags = v$lags, N = v$n,
             SR_Elasticity = v$sr_elast, LR_Elasticity = eg[[f]]$beta,
             ECM_Coef = v$alpha, Speed_Adj_Pct = v$speed_pct,
             Half_Life_Weeks = v$half_life, R_squared = v$r2, Adj_R_squared = v$adj_r2,
             row.names = NULL)
}))
print(vecm_summary)
save_table(vecm_summary, "table_vecm_summary", "VECM results (3 and 5 lags)", digits = 4)

# Full coefficient tables (main specification)
for (f in names(fuel_names)) {
  cf <- vecm[[paste0(f, "_", VECM_LAGS_MAIN)]]$robust
  cf_df <- data.frame(Term = rownames(cf), Estimate = cf[, 1], HAC_SE = cf[, 2],
                      t_value = cf[, 3], p_value = cf[, 4], Sig = stars(cf[, 4]),
                      row.names = NULL)
  save_table(cf_df, paste0("table_vecm_coefficients_", f),
             paste("VECM coefficients -", fuel_names[f]), digits = 4)
}


# ----------------------------------------------------------------------------
# 7. VECM RESIDUAL DIAGNOSTICS
# ----------------------------------------------------------------------------

section("7. VECM RESIDUAL DIAGNOSTICS")

vecm_diagnostics <- do.call(rbind, lapply(names(vecm), function(key) {
  v <- vecm[[key]]; f <- sub("_[0-9]+$", "", key); e <- residuals(v$model)
  lb10 <- Box.test(e, lag = 10, type = "Ljung-Box")
  lb20 <- Box.test(e, lag = 20, type = "Ljung-Box")
  bg4  <- bgtest(v$model, order = 4)
  data.frame(Fuel = fuel_names[f], Lags = v$lags,
             LjungBox10_stat = unname(lb10$statistic), LjungBox10_p = lb10$p.value,
             LjungBox20_stat = unname(lb20$statistic), LjungBox20_p = lb20$p.value,
             BG4_stat = unname(bg4$statistic), BG4_p = bg4$p.value, row.names = NULL)
}))
print(vecm_diagnostics)
save_table(vecm_diagnostics, "table_vecm_diagnostics",
           "VECM residual diagnostics (Ljung-Box, Breusch-Godfrey)", digits = 3)


# ----------------------------------------------------------------------------
# 8. GRANGER CAUSALITY (pairwise, first differences, HAC, 5 lags)
# ----------------------------------------------------------------------------

section("8. GRANGER CAUSALITY")

granger_test <- function(y, x, lags, label) {
  y <- as.numeric(y); x <- as.numeric(x); n <- length(y)
  d <- data.frame(y = y, x = x)
  for (i in seq_len(lags)) {
    d[[paste0("y_L", i)]] <- c(rep(NA, i), y[seq_len(n - i)])
    d[[paste0("x_L", i)]] <- c(rep(NA, i), x[seq_len(n - i)])
  }
  d <- na.omit(d)
  unres <- lm(reformulate(c(paste0("y_L", seq_len(lags)), paste0("x_L", seq_len(lags))),
                          response = "y"), data = d)
  res   <- lm(reformulate(paste0("y_L", seq_len(lags)), response = "y"), data = d)
  w <- waldtest(unres, res, vcov = vcovHAC(unres))
  data.frame(Direction = label, F_stat = w$F[2], p_value = w$`Pr(>F)`[2],
             Sig = stars(w$`Pr(>F)`[2]), Lags = lags, N = nrow(d))
}

L <- VECM_LAGS_MAIN
granger_table <- rbind(
  granger_test(d_log_gasoline, d_log_brent_3M, L, "Oil -> Gasoline"),
  granger_test(d_log_brent_3M, d_log_gasoline, L, "Gasoline -> Oil"),
  granger_test(d_log_diesel,   d_log_brent_3M, L, "Oil -> Diesel"),
  granger_test(d_log_brent_3M, d_log_diesel,   L, "Diesel -> Oil"),
  granger_test(d_log_diesel,   d_log_gasoline, L, "Gasoline -> Diesel"),
  granger_test(d_log_gasoline, d_log_diesel,   L, "Diesel -> Gasoline")
)
print(granger_table)
save_table(granger_table, "table_granger_causality",
           paste0("Granger causality tests (", L, " lags, HAC)"), digits = 3)


# ----------------------------------------------------------------------------
# 9. IMPULSE RESPONSE FUNCTIONS (VAR in differences + lagged ECT as exogenous)
# ----------------------------------------------------------------------------
# Brent is ordered first (Cholesky). The shock is a one-standard-deviation
# orthogonalised innovation to d(log Brent); responses are in log-differences.

section("9. IMPULSE RESPONSE FUNCTIONS")

irf_for_fuel <- function(f) {
  d <- vecm[[paste0(f, "_", VECM_LAGS_MAIN)]]$data
  Y <- data.frame(d_brent = d$d_brent, d_fuel = d$d_fuel)
  X <- matrix(d$ect_lag, ncol = 1, dimnames = list(NULL, "ect_lag"))
  v <- do.call(VAR, list(y = Y, p = VECM_LAGS_MAIN, type = "const", exogen = X))
  run <- function(cum) irf(v, impulse = "d_brent", response = "d_fuel",
                           n.ahead = IRF_HORIZON, boot = TRUE, ci = 0.95,
                           runs = IRF_RUNS, cumulative = cum)
  to_df <- function(o, label) data.frame(
    weeks = 0:IRF_HORIZON, irf = c(o$irf[[1]]),
    lower = c(o$Lower[[1]]), upper = c(o$Upper[[1]]),
    fuel = label, stringsAsFactors = FALSE)
  list(normal = to_df(run(FALSE), fuel_names[f]),
       cumulative = to_df(run(TRUE), fuel_names[f]))
}
irf_res <- lapply(names(fuel_names), irf_for_fuel)
df_irf     <- do.call(rbind, lapply(irf_res, `[[`, "normal"))
df_irf_cum <- do.call(rbind, lapply(irf_res, `[[`, "cumulative"))

plot_irf <- function(df, title, ylab) {
  ggplot(df, aes(x = weeks)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = fuel), alpha = 0.2) +
    geom_line(aes(y = irf, color = fuel), linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    scale_color_manual(values = fuel_colors) + scale_fill_manual(values = fuel_colors) +
    labs(title = title,
         subtitle = "Orthogonalised Brent shock, 95% bootstrap confidence band",
         x = "Weeks after shock", y = ylab, color = "Fuel", fill = "Fuel") +
    theme_minimal() + theme(legend.position = "bottom")
}
save_plot(plot_irf(df_irf, "Impulse response: Brent -> retail fuel prices", "Response"),
          "fig_irf_normal")
save_plot(plot_irf(df_irf_cum, "Cumulative impulse response", "Cumulative response"),
          "fig_irf_cumulative")

key_weeks <- c(0, 1, 2, 4, 8, 12, 20)
irf_key <- rbind(transform(df_irf[df_irf$weeks %in% key_weeks, ], type = "Normal"),
                 transform(df_irf_cum[df_irf_cum$weeks %in% key_weeks, ], type = "Cumulative"))
save_table(irf_key, "table_irf_key_weeks", "IRF at selected horizons", digits = 4)


# ----------------------------------------------------------------------------
# 10. FORWARD VS SPOT OIL PRICES (robustness of the cost proxy)
# ----------------------------------------------------------------------------
# Bivariate VAR in log-levels (4 lags): response of retail prices to a Brent
# shock when Brent is measured with 3M futures vs spot prices.

if (RUN_FORWARD_VS_SPOT) {
  section("10. FORWARD VS SPOT")
  irf_levels <- function(log_fuel_series, log_brent_series, fuel_label, brent_label) {
    v <- VAR(cbind(brent = log_brent_series, fuel = log_fuel_series),
             p = 4, type = "const")
    o <- irf(v, impulse = "brent", response = "fuel", n.ahead = IRF_HORIZON,
             boot = TRUE, ci = 0.95, runs = IRF_RUNS)
    data.frame(weeks = 0:IRF_HORIZON, irf = c(o$irf[[1]]),
               lower = c(o$Lower[[1]]), upper = c(o$Upper[[1]]),
               fuel = fuel_label, brent = brent_label)
  }
  df_fs <- rbind(
    irf_levels(log_gasoline, log_brent_3M,   "Gasoline", "3M futures"),
    irf_levels(log_gasoline, log_brent_spot, "Gasoline", "Spot"),
    irf_levels(log_diesel,   log_brent_3M,   "Diesel",   "3M futures"),
    irf_levels(log_diesel,   log_brent_spot, "Diesel",   "Spot")
  )
  p_fs <- ggplot(df_fs, aes(x = weeks, y = irf, color = brent, fill = brent)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    facet_wrap(~ fuel) +
    scale_color_manual(values = c("3M futures" = "#2ECC71", "Spot" = "#F39C12")) +
    scale_fill_manual(values  = c("3M futures" = "#2ECC71", "Spot" = "#F39C12")) +
    labs(title = "Response of retail prices to a Brent shock: futures vs spot",
         x = "Weeks after shock", y = "Response (log-levels)", color = "Oil price",
         fill = "Oil price") +
    theme_minimal() + theme(legend.position = "bottom")
  save_plot(p_fs, "fig_forward_vs_spot_irf", width = 11)
  save_table(df_fs[df_fs$weeks %in% key_weeks, ], "table_forward_vs_spot_irf",
             "IRF with futures vs spot oil prices", digits = 4)
}


# ----------------------------------------------------------------------------
# 11. NARDL(6,6,6)  (Shin, Yu and Greenwood-Nimmo, 2014)
# ----------------------------------------------------------------------------
#   d(y_t) = c + phi*y_{t-1} + th+*x+_{t-1} + th-*x-_{t-1}
#            + sum_{i=1..p} rho_i d(y_{t-i})
#            + sum_{i=0..q} (pi+_i d(x+_{t-i}) + pi-_i d(x-_{t-i})) + e_t
# x+ / x- are the partial sums of positive / negative changes of log Brent (3M).

section("11. NARDL")

delta_brent <- as.numeric(diff(log_brent_3M))
brent_pos <- c(0, cumsum(pmax(delta_brent, 0)))
brent_neg <- c(0, cumsum(pmin(delta_brent, 0)))

nardl_data <- ts(cbind(gasoline = as.numeric(log_gasoline),
                       diesel   = as.numeric(log_diesel),
                       brent_pos = brent_pos, brent_neg = brent_neg),
                 start = c(2005, 1), frequency = 52)

fit_nardl <- function(dep, p, q, start = NULL) {
  fml <- as.formula(sprintf(
    "d(%s) ~ L(%s, 1) + L(brent_pos, 1) + L(brent_neg, 1) + L(d(%s), 1:%d) + L(d(brent_pos), 0:%d) + L(d(brent_neg), 0:%d)",
    dep, dep, dep, p, q, q))
  if (is.null(start)) dynlm(fml, data = nardl_data) else dynlm(fml, data = nardl_data, start = start)
}

# Optional: information criteria over NARDL(p,p,p), p = 1..12, on a common
# sample (documents the choice of 6 lags).
if (RUN_NARDL_LAG_SELECTION) {
  max_p <- 12
  common_start <- time(nardl_data)[max_p + 2]
  nardl_ic <- do.call(rbind, lapply(1:max_p, function(p) {
    g <- fit_nardl("gasoline", p, p, start = common_start)
    d <- fit_nardl("diesel",   p, p, start = common_start)
    data.frame(Lag = p, AIC_Gasoline = AIC(g), BIC_Gasoline = BIC(g),
               AIC_Diesel = AIC(d), BIC_Diesel = BIC(d))
  }))
  print(nardl_ic)
  cat("NARDL lag minimising AIC / BIC - Gasoline:", nardl_ic$Lag[which.min(nardl_ic$AIC_Gasoline)],
      "/", nardl_ic$Lag[which.min(nardl_ic$BIC_Gasoline)], "\n")
  cat("NARDL lag minimising AIC / BIC - Diesel:  ", nardl_ic$Lag[which.min(nardl_ic$AIC_Diesel)],
      "/", nardl_ic$Lag[which.min(nardl_ic$BIC_Diesel)], "\n")
  save_table(nardl_ic, "table_nardl_lag_selection", "NARDL(p,p,p) information criteria", digits = 2)
}

# --- 11.1 Engle-Granger check on the ECT (no drift, 4 lags) -----------------
NARDL_EG_CV5 <- -3.37
nardl_eg <- data.frame(
  Model = fuel_names,
  ADF_stat = sapply(eg, function(e) ur.df(e$ect, type = "none", lags = 4)@teststat[1]),
  Critical_5pct = NARDL_EG_CV5, row.names = NULL
)
nardl_eg$Cointegrated <- nardl_eg$ADF_stat < nardl_eg$Critical_5pct
print(nardl_eg)
save_table(nardl_eg, "table_nardl_engle_granger", "Engle-Granger test on the ECT", digits = 3)

# --- 11.2 Estimation ---------------------------------------------------------
nardl <- list(gasoline = fit_nardl("gasoline", NARDL_LAGS, NARDL_LAGS),
              diesel   = fit_nardl("diesel",   NARDL_LAGS, NARDL_LAGS))
for (f in names(nardl)) {
  cat("\n---", fuel_names[f], ": NARDL(", NARDL_LAGS, ",", NARDL_LAGS, ",", NARDL_LAGS, ") ---\n")
  print(summary(nardl[[f]]))
}

# Coefficients with classical and HAC standard errors
for (f in names(nardl)) {
  cl <- coeftest(nardl[[f]]); hc <- coeftest(nardl[[f]], vcov = vcovHAC(nardl[[f]]))
  tab <- data.frame(Term = rownames(cl), Estimate = cl[, 1],
                    SE = cl[, 2], p_value = cl[, 4],
                    HAC_SE = hc[, 2], HAC_p_value = hc[, 4], row.names = NULL)
  save_table(tab, paste0("table_nardl_coefficients_", f),
             paste("NARDL coefficients -", fuel_names[f]), digits = 4)
}

# --- 11.3 Bounds test and asymmetry (Wald) tests ----------------------------
NARDL_BOUNDS_5PCT <- c(I0 = 3.79, I1 = 4.85)   # Pesaran et al. (2001), k = 2

restriction_matrix <- function(model, plus = character(), minus = character()) {
  nm <- names(coef(model))
  R <- matrix(0, nrow = 1, ncol = length(nm), dimnames = list(NULL, nm))
  R[1, plus] <- 1; R[1, minus] <- -1
  R
}
bounds_matrix <- function(model, dep) {
  nm <- names(coef(model))
  R <- matrix(0, nrow = 3, ncol = length(nm), dimnames = list(NULL, nm))
  R[1, sprintf("L(%s, 1)", dep)] <- 1
  R[2, "L(brent_pos, 1)"] <- 1
  R[3, "L(brent_neg, 1)"] <- 1
  R
}
wald_F <- function(model, R, robust = FALSE) {
  rhs <- rep(0, nrow(R))
  out <- if (robust) linearHypothesis(model, R, rhs = rhs, vcov. = vcovHAC)
         else        linearHypothesis(model, R, rhs = rhs)
  c(F = out$F[2], p = out$`Pr(>F)`[2])
}

wald_rows <- list()
for (f in names(nardl)) {
  m <- nardl[[f]]
  nm <- names(coef(m))
  pos_sr <- grep("^L\\(d\\(brent_pos\\)", nm, value = TRUE)
  neg_sr <- grep("^L\\(d\\(brent_neg\\)", nm, value = TRUE)
  R_lr <- restriction_matrix(m, plus = "L(brent_pos, 1)", minus = "L(brent_neg, 1)")
  R_sr <- restriction_matrix(m, plus = pos_sr, minus = neg_sr)   # sum of pos = sum of neg
  b    <- wald_F(m, bounds_matrix(m, f))
  wald_rows[[f]] <- data.frame(
    Fuel = fuel_names[f],
    Bounds_F = unname(b["F"]), Bounds_p = unname(b["p"]),
    LR_asym_F = unname(wald_F(m, R_lr)["F"]),
    LR_asym_p = unname(wald_F(m, R_lr)["p"]),
    LR_asym_p_HAC = unname(wald_F(m, R_lr, robust = TRUE)["p"]),
    SR_asym_F = unname(wald_F(m, R_sr)["F"]),
    SR_asym_p = unname(wald_F(m, R_sr)["p"]),
    SR_asym_p_HAC = unname(wald_F(m, R_sr, robust = TRUE)["p"]),
    row.names = NULL)
}
wald_table <- do.call(rbind, wald_rows)
wald_table$Bounds_conclusion <- ifelse(wald_table$Bounds_F > NARDL_BOUNDS_5PCT["I1"],
                                       "Cointegrated (F > I(1) bound, 5%)",
                                       ifelse(wald_table$Bounds_F < NARDL_BOUNDS_5PCT["I0"],
                                              "No cointegration", "Inconclusive"))
cat("5% critical bounds (k = 2): I(0) =", NARDL_BOUNDS_5PCT["I0"],
    ", I(1) =", NARDL_BOUNDS_5PCT["I1"], "\n")
print(wald_table)
save_table(wald_table, "table_nardl_bounds_and_wald",
           "NARDL bounds test and asymmetry tests (classical and HAC p-values)", digits = 4)

# --- 11.4 Long-run multipliers, dynamic multipliers, speed of adjustment ----
# Dynamic multipliers by forward simulation of the estimated NARDL after a
# permanent shock to Brent (positive: x+ rises; negative: x- falls).

dynamic_multiplier <- function(model, dep, shock = c("pos", "neg"),
                               size = DM_SHOCK, horizon = DM_HORIZON) {
  shock <- match.arg(shock)
  b <- coef(model); nm <- names(b)
  phi   <- unname(b[sprintf("L(%s, 1)", dep)])
  b_pos <- unname(b["L(brent_pos, 1)"]); b_neg <- unname(b["L(brent_neg, 1)"])
  rho   <- unname(b[grep(sprintf("^L\\(d\\(%s\\)", dep), nm)])
  g_pos <- unname(b[grep("^L\\(d\\(brent_pos\\)", nm)])
  g_neg <- unname(b[grep("^L\\(d\\(brent_neg\\)", nm)])
  const <- if ("(Intercept)" %in% nm) unname(b["(Intercept)"]) else 0
  nA <- length(rho); nB <- length(g_pos)

  burn  <- max(nA, nB) + 5
  total <- burn + horizon + 1
  t0    <- burn + 1                       # period of the shock

  simulate <- function(pos, neg) {
    y <- numeric(total)
    jB <- 0:(nB - 1)
    for (t in (burn + 1):total) {
      ecm <- phi * y[t - 1] + b_pos * pos[t - 1] + b_neg * neg[t - 1]
      ar  <- sum(rho * (y[t - seq_len(nA)] - y[t - seq_len(nA) - 1]))
      sr_pos <- sum(g_pos * (pos[t - jB] - pos[t - jB - 1]))
      sr_neg <- sum(g_neg * (neg[t - jB] - neg[t - jB - 1]))
      y[t] <- y[t - 1] + const + ecm + ar + sr_pos + sr_neg
    }
    y
  }

  pos_shock <- neg_shock <- numeric(total)
  if (shock == "pos") pos_shock[t0:total] <- size else neg_shock[t0:total] <- -size
  y_base  <- simulate(numeric(total), numeric(total))
  y_treat <- simulate(pos_shock, neg_shock)

  diff_path <- (y_treat - y_base)[t0 + 0:horizon]
  dm  <- if (shock == "neg") abs(diff_path) / abs(size) else diff_path / abs(size)
  lrm <- if (shock == "pos") -b_pos / phi else abs(-b_neg / phi)
  list(dm = dm, lrm = lrm, phi = phi)
}

dm <- lapply(names(nardl), function(f) {
  pos <- dynamic_multiplier(nardl[[f]], f, "pos")
  neg <- dynamic_multiplier(nardl[[f]], f, "neg")
  list(pos = pos$dm, neg = neg$dm, asym = pos$dm - neg$dm,
       lrm_pos = pos$lrm, lrm_neg = neg$lrm, phi = pos$phi)
})
names(dm) <- names(nardl)

# Dynamic multipliers at selected horizons (h = 0, 1, 4, 8, 16, 30, 50 weeks)
dm_weeks <- c(0, 1, 4, 8, 16, 30, 50)
dm_table <- do.call(rbind, lapply(names(dm), function(f) {
  d <- dm[[f]]
  rbind(
    data.frame(Fuel = fuel_names[f], Horizon = as.character(dm_weeks),
               Increase = d$pos[dm_weeks + 1], Decrease = d$neg[dm_weeks + 1],
               Asymmetry = d$asym[dm_weeks + 1], Asymmetry_pp = 100 * d$asym[dm_weeks + 1]),
    data.frame(Fuel = fuel_names[f], Horizon = "Long run",
               Increase = d$lrm_pos, Decrease = d$lrm_neg,
               Asymmetry = d$lrm_pos - d$lrm_neg,
               Asymmetry_pp = 100 * (d$lrm_pos - d$lrm_neg))
  )
}))
print(dm_table)
save_table(dm_table, "table_nardl_dynamic_multipliers",
           "NARDL dynamic multipliers (1% permanent Brent shock)", digits = 4)

# Long-run multiplier ratio (upward / downward)
lr_ratio <- data.frame(Fuel = fuel_names,
                       LRM_increase = sapply(dm, `[[`, "lrm_pos"),
                       LRM_decrease = sapply(dm, `[[`, "lrm_neg"), row.names = NULL)
lr_ratio$Ratio <- lr_ratio$LRM_increase / lr_ratio$LRM_decrease
print(lr_ratio)
save_table(lr_ratio, "table_nardl_long_run_ratio", "Long-run multipliers", digits = 4)

# Convergence check: gap from the long-run multiplier at h = 30 and 50
convergence <- do.call(rbind, lapply(names(dm), function(f) {
  d <- dm[[f]]
  data.frame(Fuel = fuel_names[f], Horizon = c(30, 50),
             Gap_increase_pct = 100 * (1 - d$pos[c(31, 51)] / d$lrm_pos),
             Gap_decrease_pct = 100 * (1 - d$neg[c(31, 51)] / d$lrm_neg))
}))
print(convergence)

# Speed of adjustment (phi) and half-life
adjust_table <- data.frame(
  Fuel = fuel_names,
  ECM_coefficient = sapply(dm, `[[`, "phi"),
  Speed = abs(sapply(dm, `[[`, "phi")),
  Half_life_weeks = log(0.5) / log(1 - abs(sapply(dm, `[[`, "phi"))),
  Pct_adjusted_4w = 100 * (1 - (1 - abs(sapply(dm, `[[`, "phi")))^4),
  Pct_adjusted_8w = 100 * (1 - (1 - abs(sapply(dm, `[[`, "phi")))^8),
  row.names = NULL
)
print(adjust_table)
save_table(adjust_table, "table_nardl_adjustment_speed",
           "Speed of adjustment to long-run equilibrium", digits = 4)

# Plot: dynamic multipliers
df_dm <- do.call(rbind, lapply(names(dm), function(f) {
  data.frame(weeks = rep(0:DM_HORIZON, 2),
             value = c(dm[[f]]$pos, dm[[f]]$neg),
             type = rep(c("Oil increase", "Oil decrease"), each = DM_HORIZON + 1),
             fuel = fuel_names[f])
}))
p_dm <- ggplot(df_dm, aes(x = weeks, y = value, color = type)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ fuel) +
  scale_color_manual(values = c("Oil increase" = "#E74C3C", "Oil decrease" = "#3498DB")) +
  labs(title = "NARDL dynamic multipliers", subtitle = "Permanent 1% Brent shock",
       x = "Weeks", y = "Cumulative response of retail price", color = "") +
  theme_minimal() + theme(legend.position = "bottom")
save_plot(p_dm, "fig_nardl_dynamic_multipliers", width = 11)

# --- 11.5 NARDL residual diagnostics -----------------------------------------
nardl_diag <- do.call(rbind, lapply(names(nardl), function(f) {
  m <- nardl[[f]]; e <- residuals(m)
  bg <- bgtest(m, order = 4); jb <- jarque.bera.test(as.numeric(e)); bp <- bptest(m)
  data.frame(Fuel = fuel_names[f], N = length(e),
             Adj_R_squared = summary(m)$adj.r.squared,
             BG4_stat = unname(bg$statistic), BG4_p = bg$p.value,
             JB_stat = unname(jb$statistic), JB_p = jb$p.value,
             BP_stat = unname(bp$statistic), BP_p = bp$p.value, row.names = NULL)
}))
print(nardl_diag)
save_table(nardl_diag, "table_nardl_diagnostics",
           "NARDL diagnostics (Breusch-Godfrey, Jarque-Bera, Breusch-Pagan)", digits = 4)

png(file.path(OUT_DIR, "figures", "fig_nardl_residuals.png"),
    width = 1400, height = 1000, res = 150)
par(mfrow = c(2, 4), mar = c(4, 4, 3, 1))
for (f in names(nardl)) {
  e <- as.numeric(residuals(nardl[[f]]))
  plot(e, type = "l", main = paste(fuel_names[f], "- residuals"), ylab = "", xlab = "Index")
  abline(h = 0, col = "red", lty = 2)
  qqnorm(e, main = paste(fuel_names[f], "- Q-Q plot")); qqline(e, col = "red")
  hist(e, breaks = 30, freq = FALSE, main = paste(fuel_names[f], "- distribution"), xlab = "")
  curve(dnorm(x, mean(e), sd(e)), add = TRUE, col = "red", lwd = 2)
  acf(e, main = paste(fuel_names[f], "- ACF"))
}
dev.off()


# ----------------------------------------------------------------------------
# 12. CHECK AGAINST THE VALUES REPORTED IN THE THESIS (optional)
# ----------------------------------------------------------------------------
# Reference values below were copied from the thesis text and tables. A
# mismatch does not necessarily mean an error: it tells you which number to
# look at (data version, lag setting, SE type, rounding).

if (RUN_THESIS_CHECK) {
  section("12. CHECK AGAINST THESIS VALUES")
  v5g <- vecm[[paste0("gasoline_", VECM_LAGS_MAIN)]]
  v5d <- vecm[[paste0("diesel_",   VECM_LAGS_MAIN)]]
  gr  <- function(label) granger_table$F_stat[granger_table$Direction == label]
  ur  <- function(row, col) unit_root_table[[col]][row]

  chk <- function(metric, thesis, got, tol) {
    data.frame(Metric = metric, Thesis = thesis, Reproduced = round(got, 4),
               AbsDiff = round(abs(thesis - got), 4), Match = abs(thesis - got) <= tol)
  }
  thesis_check <- rbind(
    chk("N observations (descriptive stats)", 993, nrow(dat), 0),
    chk("ADF statistic - gasoline", -2.886, ur(1, "ADF"), 0.002),
    chk("ADF statistic - diesel",   -2.474, ur(2, "ADF"), 0.002),
    chk("ADF statistic - Brent 3M", -3.064, ur(3, "ADF"), 0.002),
    chk("Long-run elasticity - gasoline", 0.77, eg$gasoline$beta, 0.006),
    chk("Long-run elasticity - diesel",   0.81, eg$diesel$beta,   0.006),
    chk("VECM(5) ECM coefficient - gasoline", -0.0275, v5g$alpha, 0.0005),
    chk("VECM(5) ECM coefficient - diesel",   -0.0256, v5d$alpha, 0.0005),
    chk("VECM(5) adjusted R2 - gasoline", 0.484, v5g$adj_r2, 0.002),
    chk("VECM(5) adjusted R2 - diesel",   0.510, v5d$adj_r2, 0.002),
    chk("VECM(5) sample size", 987, v5g$n, 0),
    chk("Granger F: Oil -> Gasoline", 18.199, gr("Oil -> Gasoline"), 0.05),
    chk("Granger F: Oil -> Diesel",   19.901, gr("Oil -> Diesel"),   0.05),
    chk("NARDL bounds F - gasoline", 13.220, wald_table$Bounds_F[1], 0.01),
    chk("NARDL bounds F - diesel",    6.774, wald_table$Bounds_F[2], 0.01),
    chk("NARDL adjusted R2 - gasoline", 0.5057, nardl_diag$Adj_R_squared[1], 0.0005),
    chk("NARDL adjusted R2 - diesel",   0.5312, nardl_diag$Adj_R_squared[2], 0.0005),
    chk("Asymmetry at 4 weeks - gasoline", 0.2631, dm$gasoline$asym[5], 0.0005),
    chk("Asymmetry at 4 weeks - diesel",   0.2729, dm$diesel$asym[5],   0.0005),
    chk("Long-run multiplier (increase) - gasoline", 0.7414, dm$gasoline$lrm_pos, 0.0005),
    chk("Long-run multiplier (increase) - diesel",   0.9031, dm$diesel$lrm_pos,   0.0005),
    chk("NARDL half-life (weeks) - gasoline", 10.60, adjust_table$Half_life_weeks[1], 0.05),
    chk("NARDL half-life (weeks) - diesel",   19.93, adjust_table$Half_life_weeks[2], 0.05)
  )
  print(thesis_check, row.names = FALSE)
  save_table(thesis_check, "table_check_vs_thesis", "Reproduced values vs thesis", digits = 4)
  cat("\nChecks matching:", sum(thesis_check$Match), "of", nrow(thesis_check), "\n")
}

# Record the software versions used
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "session_info.txt"))
section("DONE - tables in output/tables, figures in output/figures")
