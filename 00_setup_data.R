# ==============================================================================
# 08 - GRAPHICS
# ==============================================================================
# source("R/00_setup_data.R") must be run first.
# ==============================================================================

library(scales)

# ------------------------------------------------------------------------
# 1. Log price levels over time
# ------------------------------------------------------------------------
df_levels <- data.frame(
  Date     = time(log_gasoline),
  Gasoline = as.numeric(log_gasoline),
  Diesel   = as.numeric(log_diesel),
  Brent    = as.numeric(log_brent_3M)
)

df_levels_long <- pivot_longer(df_levels, cols = c("Gasoline", "Diesel", "Brent"),
                                names_to = "Series", values_to = "Value")

ggplot(df_levels_long, aes(x = Date, y = Value, color = Series)) +
  geom_line(linewidth = 1) +
  labs(title = "Weekly Log Price Series (2005-2024)", x = "Year", y = "Log(price)") +
  theme_minimal()

# ------------------------------------------------------------------------
# 2. Highlight major oil-market crisis periods
# ------------------------------------------------------------------------
crisis_periods <- data.frame(
  label = c("2008 Financial Crisis", "COVID-19", "Russia-Ukraine War"),
  start = as.Date(c("2008-06-01", "2020-03-01", "2022-02-24")),
  end   = as.Date(c("2009-06-30", "2020-10-31", "2024-01-01"))
)

df_dates <- data.frame(
  Date = as.Date(dati_clean$Date),
  Gasoline = as.numeric(ts_gasoline),
  Diesel   = as.numeric(ts_diesel),
  Brent_3M = as.numeric(ts_brent_3M)
)
df_dates <- df_dates[!is.na(df_dates$Date), ]

df_dates_long <- pivot_longer(df_dates, cols = c("Gasoline", "Diesel", "Brent_3M"),
                               names_to = "Series", values_to = "Price")

ggplot(df_dates_long, aes(x = Date, y = Price, color = Series)) +
  geom_rect(data = crisis_periods, inherit.aes = FALSE,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = "grey85", alpha = 0.5) +
  geom_line(linewidth = 0.8) +
  labs(title = "Retail Fuel and Oil Prices with Major Crisis Periods Highlighted",
       x = "Year", y = "Price (EUR)") +
  theme_minimal()

# ------------------------------------------------------------------------
# 3. Log-differences (weekly % changes) - useful to visually inspect
#    volatility clustering ahead of the Markov-switching analysis
# ------------------------------------------------------------------------
df_changes <- data.frame(
  Date = df_dates$Date[-1],
  D_Log_Gasoline = as.numeric(d_log_gasoline),
  D_Log_Diesel   = as.numeric(d_log_diesel),
  D_Log_Brent_3M = as.numeric(d_log_brent_3M)
)

df_changes_long <- pivot_longer(df_changes,
                                 cols = c("D_Log_Gasoline", "D_Log_Diesel", "D_Log_Brent_3M"),
                                 names_to = "Series", values_to = "Change")

ggplot(df_changes_long, aes(x = Date, y = Change)) +
  geom_line(color = "steelblue") +
  facet_wrap(~Series, ncol = 1, scales = "free_y") +
  labs(title = "Weekly Log-Differences (Volatility Clustering)", x = "Year", y = "% change") +
  theme_minimal()
