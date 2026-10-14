# Asymmetric price transmission in Italian fuel markets (2005–2024)

MSc thesis in Economics and Econometrics, University of Bologna.
Weekly data, cointegration, VECM and NARDL in R.

**Main result.** Italian retail gasoline and diesel prices follow crude oil increases more than twice as fast as decreases ("rockets and feathers"). At the four-week horizon the gap between the two responses is **26 percentage points for gasoline and 27 for diesel**, and it disappears by week eight. A small but statistically significant asymmetry remains in the long run.

## Key results

| | Gasoline | Diesel |
|---|---|---|
| Long-run elasticity of retail price to Brent (Engle-Granger) | 0.77 | 0.81 |
| NARDL(6,6,6) bounds test, F-statistic (5% upper bound 4.85) | 13.22 | 6.77 |
| Pass-through gap, oil increase minus decrease, at 4 weeks | 26.3 pp | 27.3 pp |
| Long-run multiplier after an oil increase | 0.741 | 0.903 |
| Half-life of deviations from equilibrium (NARDL) | 10.6 weeks | 19.9 weeks |

Oil Granger-causes both retail prices (F = 18.2 for gasoline, 19.9 for diesel; 5 lags, HAC).

## What the code does

`thesis_pipeline.R` runs the whole analysis, in the order of the thesis chapters:

1. Data preparation: Brent converted to EUR with maturity-matched exchange rates
2. Descriptive statistics and figures
3. Unit root tests (ADF, Phillips-Perron, KPSS, Zivot-Andrews)
4. VAR lag selection
5. Engle-Granger cointegration
6. VECM (5 lags main specification, 3 lags as robustness)
7. Residual diagnostics
8. Granger causality
9. Impulse response functions
10. Futures versus spot oil prices
11. NARDL(6,6,6) with bounds test, asymmetry tests, dynamic multipliers
12. Automatic check of the results against the values reported in the thesis

## How to reproduce

1. Put the data file `Dati_completo.xlsx` in the same folder as the script (see *Data*).
2. In R:

```r
setwd("path/to/MasterThesis")
install.packages(c("readxl", "urca", "tseries", "vars", "lmtest", "sandwich",
                   "car", "dynlm", "xtable", "ggplot2", "moments"))
source("thesis_pipeline.R")
```

Tables (`.csv` and `.tex`) are written to `output/tables` and figures to `output/figures`. The script was developed with R 4.5; `output/session_info.txt` records the package versions of each run. Bootstrap confidence bands use a fixed seed. With `IRF_RUNS <- 50` a first run is much faster than the default 500 replications.

### Reproducibility check

Section 12 of the script compares 23 reference values with those in the thesis. On the original data, **21 match to the fourth decimal** (VECM, Granger, bounds test, NARDL, multipliers, half-lives). The two remaining differences come from a labelling error in the thesis, described under *Errata*. The comparison is saved in `output/tables/table_check_vs_thesis.csv`.

## Data availability and reproducibility

The analysis is fully scripted: `thesis_pipeline.R` goes from the raw input
file to every table and figure. The input file (`Dati_completo.xlsx`) is not
included, because the Brent futures, Brent spot and exchange-rate series come
from LSEG and cannot be redistributed. The retail fuel prices come from the
EC Weekly Oil Bulletin, which is public.

To reproduce the results you need a file with the structure described below
(for example, built from an LSEG subscription). On the original data, the
script reproduces the thesis results; the comparison is saved in
`output/tables/table_check_vs_thesis.csv` (21 of 23 reference values match, the
other two are explained under Errata).

| Column | Content | Source |
|---|---|---|
| `Date` | Week | |
| `Gasoline`, `Diesel` | Italian pre-tax retail prices, EUR per 1,000 litres | EC Weekly Oil Bulletin (public) |
| `ICE - Brent  Continuous 3rd future Settlment Price` | Brent 3rd continuous future, USD/bbl | LSEG |
| `EXCHANGE RATE 3M FWD` | USD per EUR, 3-month forward | LSEG |
| `Brent Spot USD/BBL` | Brent spot, USD/bbl | LSEG |
| `USD TO EUR (LSEG DS) - EXCHANGE RATE` | USD per EUR, spot | LSEG |

(The double space in the third column name is as in the source file.) Weeks with missing values are dropped, not interpolated; the final sample has 993 weekly observations.

## Repository structure

```
thesis_pipeline.R      # full analysis, single script
output/                # tables and figures produced by the script
exploratory/           # additional analyses NOT used for the thesis results
Asymmetric price transmission ... .pdf   # the thesis
```

`exploratory/` holds an asymmetric error-correction model and a two-regime Markov-switching model that were explored alongside the thesis. They are kept for reference: the thesis only mentions regime switching as a direction for future work.

## Errata

Comparing the code output with the thesis text turned up three points in the unit root section (Table 3 and its discussion):

1. **Gasoline and Brent rows swapped.** In Table 3 all four statistics of gasoline and Brent are interchanged (diesel is correct). The correct values are:

   | | ADF | PP | KPSS | Zivot-Andrews |
   |---|---|---|---|---|
   | Gasoline | −3.063 | −2.816 | 2.342 | −3.339 |
   | Diesel | −2.474 | −2.315 | 1.604 | −3.131 |
   | Brent (3M) | −2.886 | −3.111 | 0.883 | −4.333 |

   The statements in the text that quote these figures (ADF and PP "strongest evidence" for gasoline, KPSS values, the Zivot-Andrews range) follow the swapped labels.
2. **KPSS critical value.** The value shown as the 1% critical value (0.574) is the 2.5% one; the 1% value is 0.739. All three KPSS statistics still exceed it, so the conclusion is unchanged.
3. **Wording on rejection.** The text states that neither ADF nor PP rejects the unit root at conventional levels. With the correct statistics, the ADF test rejects at the 5% level (not at 1%) for gasoline and Brent, and the PP test rejects at 5% for Brent only; diesel is not rejected by either. KPSS and Zivot-Andrews, together with stationarity in first differences, still support treating the series as I(1).

Impulse responses are computed for a one-standard-deviation orthogonalised shock to Brent (Cholesky ordering, Brent first), not a 1% shock.

## References

Shin, Y., Yu, B. and Greenwood-Nimmo, M. (2014). Modelling asymmetric cointegration and dynamic multipliers in a nonlinear ARDL framework.
Pesaran, M.H., Shin, Y. and Smith, R.J. (2001). Bounds testing approaches to the analysis of level relationships.
Engle, R.F. and Granger, C.W.J. (1987). Co-integration and error correction.

## License

MIT (code). The thesis PDF and the LSEG data are not covered by this license.
