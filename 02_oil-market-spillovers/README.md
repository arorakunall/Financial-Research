# 02 · Liquidity and Volatility Spillovers in the Global Oil Market

**Status:** Completed — Time Series Analytics Project
**Submitted to:** Dr. Sahaj Wadhwa, Delhi School of Economics, University of Delhi
**Period:** April 2026

---

## Overview

This study examines **price return, realized volatility, and trading volume spillovers** between WTI and Brent crude oil futures over 10 Indian financial years (April 2016 – March 2026), covering 2,457 daily trading observations. The analysis spans three major global crises — the COVID-19 pandemic (2020), the Russia-Ukraine war (2022), and the Iran-related Middle East conflict (2023).

The core question: Who leads and who follows between WTI and Brent — and does that leadership hold across returns, volatility, and liquidity simultaneously?

---

## Key Findings

- **Price Leadership:** Brent explains **86.74%** of WTI return variance while WTI explains only 0.10% of Brent's — establishing Brent as the dominant, unidirectional price transmitter.
- **Volatility Dominance:** Brent_Vol Granger-causes WTI_Vol (chi-sq = 150.28, p < 0.0001) and explains **59.86%** of WTI_Vol's forecast error variance. The reverse channel is weak (3.25%).
- **Liquidity Transmission:** WTI volume Granger-causes Brent volume (p = 0.0099), demonstrating cross-market liquidity spillovers. Volume series remain 95–98% own-driven in FEVD — disturbances are transitory.
- **Total Connectedness Index:** Returns 43.42% · Volatility 31.56% · Liquidity 1.37% — revealing a clear integration hierarchy across market dimensions.
- **Crisis Robustness:** COVID-19, Russia-Ukraine, and Iran event dummies absorbed endogenously within VAR dynamics (R² = 97.3–97.8% for volatility equations). IRF confirms shocks dissipate within 4–5 trading days.

---

## Data

| File | Description |
|---|---|
| `oil_processed.csv` | Cleaned daily dataset — WTI and Brent log returns, 21-day rolling realized volatility (annualized), log-volume changes, and crisis dummies |

**Source:** Yahoo Finance (via `yfinance`)
**Period:** April 2016 – March 2026 · **Observations:** 2,457 · **Frequency:** Daily

---

## Methodology

| Step | Tool | Detail |
|---|---|---|
| Data collection & prep | Python (`notebooks/`) | Log returns, rolling volatility, log-volume changes, crisis dummies |
| Stationarity testing | EViews | ADF tests — all 6 series confirmed I(0) |
| Lag selection | EViews | SIC criterion → lag 4 selected |
| VAR estimation | EViews | VAR(4), 6 variables, with COVID/Russia/Iran dummies |
| VAR stability | EViews | Inverse AR roots — all within unit circle |
| Granger causality | EViews | Pairwise Block Exogeneity Wald Tests |
| Connectedness | EViews | Diebold-Yilmaz (2012) via 10-day Cholesky FEVD |
| Shock dynamics | EViews | Impulse Response Functions (10-day horizon) |
| Visualisation | Python (`notebooks/`) | Price, return, and volatility plots with crisis markers |

---

## Repository Structure

```
02_oil-market-spillovers/
│
├── README.md
│
├── data/
│   └── oil_processed.csv           ← cleaned dataset (all 6 variables + dummies)
│
├── notebooks/
│   └── Time_Series_Project.ipynb   ← data prep, feature construction, visualisation
│
└── outputs/
    ├── figures/
    │   ├── oil_market_plots.png                    ← prices, returns, volatility (2016–2026)
    │   ├── Impulse_Response_Function.pdf           ← individual IRF panels
    │   ├── Impulse_Response_Function_Stacked.pdf   ← stacked IRF grid
    │   ├── Variance_Decomposition.pdf              ← FEVD stacked bar charts
    │   ├── Variance_decomposition_Stack.pdf        ← alternative FEVD layout
    │   └── AR_Graph___Lag_1_4.pdf                  ← VAR stability — inverse roots
    │
    ├── tables/
    │   ├── Descriptive_Stats.pdf
    │   ├── Lag_Choosing_Criteria.pdf
    │   ├── VAR_Equation.pdf                        ← full VAR(4) coefficient table
    │   ├── Granger_Causality.pdf                   ← Block Exogeneity Wald Test results
    │   ├── variance_decomposition_table.pdf        ← numerical FEVD table (10-day)
    │   ├── BRENT_RETURN_ADF_TEST.pdf
    │   ├── BRENT_VOL_ADF_TEST.pdf
    │   ├── BRENT_VOLUME_CHANGE_ADF_TEST.pdf
    │   ├── WTI_RETURN_ADF_TEST.pdf
    │   ├── WTI_VOL_ADF_TEST.pdf
    │   └── WTI_VOLUME_CHANGE_ADF_TEST.pdf
    │
    └── Oil_Spillovers.pdf                       ← full report
```

---

## Dependencies

```
pandas
numpy
matplotlib
yfinance
openpyxl
jupyter
```

Install all at once:
```bash
pip install pandas numpy matplotlib yfinance openpyxl jupyter
```

> **Note:** VAR estimation, ADF tests, Granger causality, FEVD, and IRF were conducted in **EViews**. The Jupyter notebook handles data collection, preprocessing, and visualisation only.

---

## Reference

Diebold, F. X., & Yilmaz, K. (2012). Better to give than to receive: Predictive directional measurement of volatility spillovers. *International Journal of Forecasting, 28*(1), 57–66.
