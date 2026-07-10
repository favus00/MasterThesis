# ==============================================================================
# RUN ALL - executes the full analysis pipeline in order
# ==============================================================================
# Run from the project root: Rscript R/run_all.R
# Requires Dati_completo.xlsx in the project root (see README).
# ==============================================================================

if (!dir.exists("output")) dir.create("output")

scripts <- c(
  "R/00_setup_data.R",
  "R/01_descriptive_stats.R",
  "R/02_stationarity_cointegration.R",
  "R/03_vecm.R",
  "R/04_asymmetric_models.R",
  "R/05_markov_switching.R",
  "R/06_granger_causality.R",
  "R/07_forecasting_irf.R",
  "R/08_graphics.R"
)

for (s in scripts) {
  cat("\n\n==========", s, "==========\n")
  source(s)
}
