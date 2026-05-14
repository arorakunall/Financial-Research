# 01 · Directional Prediction of Indian Thematic Equity Indices Using ML Ensemble Classifiers

**Status:** Under Departmental Review, Department of Commerce, Delhi School of Economics, University of Delhi
**Supervisor:** Dr. Vibhuti Vashisth
**Period:** September 2025 – April 2026
**Preprint:** [SSRN — Abstract ID 6748042](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=6748042)

---

## Overview

This study evaluates the directional price-movement predictability of five Nifty thematic equity indices — **Housing, Consumption, PSE, Transportation & Logistics, and ESG** — using Random Forest (RF) and Gradient Boosting (GB) ensemble classifiers, benchmarked against Logistic Regression.

The core question: Can machine learning models predict the direction of thematic index returns? And does predictability improve at longer forecast horizons?

---

## Key Findings

- **Horizon Effect:** 1-day GB accuracy approximates the random-walk baseline (50–62%). At the 20-day horizon, accuracy rises to **90.8–95.1%** across all five indices — an improvement of 29–45 percentage points — with ROC-AUC of 93.3–99.0%.
- **Universal Super-Feature:** MA_200 ranks #1 and Price_vs_MA200 ranks #2 or #3 across all five indices via Mutual Information ranking and GB feature importance — a finding that extends Sadorsky (2021) from clean energy stocks to Indian thematic markets.
- **Backtest:** Out-of-sample simulation yields positive excess returns for Housing (+8.9 pp, Sharpe 1.26), PSE (+11.7 pp, Sharpe 0.91), and Consumption (+1.0 pp). 5-day holding periods unlock Transportation alpha (Sharpe 1.72).
- **Position Sizing:** At GB precision of 90.3–95.8%, the Kelly Criterion supports Half-Kelly fractions of 42–47% of risk capital per signal — consistent with institutional discretionary allocation norms.

---

## Data

Five daily OHLCV price series sourced from NSE India, stored as individual Excel files:

| File | Index | Observations |
|---|---|---|
| `Nifty-Housing-PR.xlsx` | Nifty Housing | 939 |
| `Nifty-INDIA-Consumption.xlsx` | Nifty India Consumption | 939 |
| `Nifty-PSE.xlsx` | Nifty PSE | 939 |
| `Nifty-Trans-Log-PR.xlsx` | Nifty Transportation & Logistics | 919 |
| `Nifty-100-ESG.xlsx` | Nifty100 ESG | 935 |

**Period:** April 2021 – October 2025 · **Frequency:** Daily

---

## Feature Engineering

41 technical indicators across 6 categories, with a two-stage selection pipeline retaining 18–20 features per index:

| Category | Features |
|---|---|
| Trend | MA-20, MA-50, MA-200, Price_vs_MA20/50/200 |
| Bollinger Bands | BollW, %B |
| Momentum | MACD, RSI-14, Stochastic %K/%D, Williams %R, CCI, ROC (3 scales) |
| Volume | OBV, OBV_ROC5, Vol_Ratio |
| Volatility | Rolling std (5/10/20-day), ATR14 |

**Selection pipeline:** Stage 1 — Mutual Information ranking · Stage 2 — GB feature importance

---

## Methodology

```
Raw OHLCV (Excel)
      │
      ▼
Feature Engineering (41 indicators)
      │
      ▼
Two-Stage Feature Selection
      │
      ▼
Model Training
Random Forest · Gradient Boosting · Logistic Regression
      │
      ▼
Walk-Forward Expanding-Window Validation
      │
      ▼
Evaluation: Accuracy · ROC-AUC · Precision · Sharpe
      │
      ▼
SHAP Interpretability · Kelly Criterion Position Sizing
```

**Forecast horizons:** 1, 5, 10, 15, 20 trading days
**Validation:** Chronological out-of-sample split; walk-forward expanding window
**Robustness:** 140 parameter configurations tested

---

## Repository Structure

```
01_thematic-equity-ml-forecasting/
│
├── README.md
│
├── data/
│   ├── Nifty-Housing-PR.xlsx
│   ├── Nifty-INDIA-Consumption.xlsx
│   ├── Nifty-PSE.xlsx
│   ├── Nifty-Trans-Log-PR.xlsx
│   └── Nifty-100-ESG.xlsx
│
├── notebooks/
│   └── Indian_Thematic_ML_Final.ipynb
│
└── outputs/
    ├── figures/              ← SHAP plots, ROC curves, accuracy charts, backtest graphs
    └── report/               ← full paper PDF
```

---

## Dependencies

```
pandas
numpy
scikit-learn
xgboost
shap
matplotlib
seaborn
openpyxl
jupyter
```

Install all at once:
```bash
pip install pandas numpy scikit-learn xgboost shap matplotlib seaborn openpyxl jupyter
```

---

## Citation

> Arora, K. (2026). *Directional Prediction of Indian Thematic Equity Indices Using Machine Learning Ensemble Classifiers: A Multi-Horizon Analysis.* Under Departmental Review, Department of Commerce, Delhi School of Economics, University of Delhi. Available at SSRN: https://papers.ssrn.com/sol3/papers.cfm?abstract_id=6748042
