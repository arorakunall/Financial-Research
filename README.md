# Financial-Research

Empirical research in quantitative finance — ML ensemble forecasting, time-series econometrics, and financial risk modelling.

**Python · R · EViews · scikit-learn · XGBoost · SHAP · NARDL · VAR**

---

## About

This repository collects my research in empirical finance, spanning machine-learning-based directional forecasting, econometric spillover and connectedness modelling, and nonlinear commodity–equity transmission. Each project is self-contained with its own code, outputs, and write-up.

I am an MBA student in Business Analytics at the Delhi School of Economics, University of Delhi (FRM Part I), with research interests in financial risk modelling, statistical learning under distribution shift, market-volatility dynamics, and risk-aware portfolio construction.

---

## Research at a glance

| # | Project | Core method | Asset class | Headline result | Output |
| --- | --- | --- | --- | --- | --- |
| 01 | Thematic Equity ML Forecasting | RF / Gradient Boosting + SHAP | Indian thematic indices | Directional accuracy rises to **90–95% at the 20-day horizon** | [SSRN preprint](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=6748042) |
| 02 | Oil-Market Spillovers | VAR(4) + Diebold–Yilmaz connectedness | WTI / Brent crude | Return spillover TCI **43.4%**; volatility **31.6%**; liquidity near-segmented | Project report |
| 03 | Mining-Equity Asymmetry | NARDL + Markov-switching | Australian iron-ore equities | Long-run asymmetry (**L⁺ 0.72 vs L⁻ 0.52**) is regime-invariant; short-run intensity triples in stress | [SSRN preprint](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=6873399) |

---

## Projects

### [01 · Thematic Equity ML Forecasting](./01_thematic-equity-ml-forecasting)

**Directional Prediction of Indian Thematic Equity Indices Using Machine Learning Ensemble Classifiers: A Multi-Horizon Analysis**

Evaluates the directional price-movement predictability of five Nifty thematic indices — Housing, Consumption, PSE, Transportation & Logistics, and ESG — using Random Forest and Gradient Boosting ensemble classifiers across five forecast horizons (1–20 trading days).

The pipeline engineers 41 technical and statistical indicators, applies two-stage feature selection, and trains ensemble classifiers benchmarked against logistic regression. The headline finding is a **"horizon effect": directional accuracy climbs steadily with the forecast horizon, reaching 90–95% at 20 days**, with the 200-day moving average emerging as the single most important feature across every index (via SHAP). A random 80/20 split (consistent with Sadorsky, 2021) is the primary design, with walk-forward validation and Kelly-criterion position sizing as robustness layers.

📄 **Preprint:** [SSRN — Abstract ID 6748042](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=6748042)

`Python` `scikit-learn` `XGBoost` `SHAP` `Random Forest` `Gradient Boosting` `Walk-Forward Validation`

---

### [02 · Oil-Market Liquidity & Volatility Spillovers](./02_oil-market-spillovers)

**Liquidity and Volatility Spillovers in the Global Oil Market: A Time-Series Analytics Study**

Examines price-return, realized-volatility, and trading-volume spillovers between WTI and Brent crude-oil futures over ten Indian financial years (April 2016 – March 2026) using a six-variable VAR(4) system and the Diebold–Yilmaz (2012) connectedness framework.

The study estimates three separate bivariate total connectedness indices — **returns (43.4%), volatility (31.6%), and liquidity (1.4%)** — showing that the two benchmarks are tightly integrated in price and risk but almost segmented in liquidity. Connectedness is tracked through three major global crises (COVID-19, the Russia–Ukraine war, and the 2023 Iran-related Middle East conflict), with Granger-causality, Cholesky FEVD, and impulse-response analysis characterising the direction and persistence of transmission. Completed as a Time Series Analytics coursework project.

`Python` `EViews` `VAR` `Granger Causality` `Cholesky FEVD` `IRF`

---

### [03 · Asymmetric Iron-Ore Shock Transmission in Australian Mining Equities](./03_mining-asymmetry-nardl)

**A Nonlinear ARDL and Regime-Switching Analysis**

Tests whether iron-ore price shocks transmit asymmetrically into the equity returns of BHP, Rio Tinto, and Fortescue over January 2011 – May 2026 (N = 3,013 trading days), and whether any asymmetry is structural or crisis-driven.

A NARDL framework decomposes cumulative iron-ore changes into positive and negative partial sums; a two-state Markov-switching model identifies calm and stress regimes, which are embedded back into the NARDL via filtered-probability interaction terms. The central result is a robust **long-run asymmetry — L⁺ = 0.720 vs L⁻ = 0.523** (Wald F = 13.90, p < 0.001) — that is **regime-invariant** and survives a two-stage block bootstrap propagating first-stage Markov uncertainty, while **short-run transmission intensity roughly triples in the stress regime**. The asymmetry concentrates by business model: the pure-play producer (Fortescue) responds three-to-four times as strongly as the diversified majors. Cointegration evidence is framed cautiously, resting on a significant error-correction term and the bounds t-test rather than a marginal F-statistic.

📄 **Preprint:** [SSRN — Abstract ID 6873399](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=6873399)

`R` `NARDL` `Markov-Switching` `Bounds Testing` `Block Bootstrap` `HAC Inference`

---

## Repository structure

```
financial-research/
│
├── README.md
│
├── 01_thematic-equity-ml-forecasting/
│   ├── README.md
│   ├── data/                  ← 5 index Excel files
│   ├── notebooks/             ← Jupyter notebook
│   └── outputs/               ← figures and report PDF
│
├── 02_oil-market-spillovers/
│   ├── README.md
│   ├── data/                  ← processed CSV
│   ├── notebooks/             ← Jupyter notebook
│   └── outputs/
│       ├── figures/           ← price/return/volatility plots
│       ├── tables/            ← ADF, VAR, Granger, FEVD, IRF (PDFs)
│       └── TimeSeriesProject.pdf
│
└── 03_mining-asymmetry-nardl/
    ├── README.md
    ├── code/                  ← R pipeline (self-installing deps)
    └── outputs/
        ├── mining_asymmetry_SSRN.pdf
        ├── tables/            ← descriptives, ADF, regime probs, lag grid, bootstrap CSVs
        └── figures/           ← partial sums, multipliers, regime probability, CUSUM, etc.
```

---

## Skills & tools

| Category | Tools |
| --- | --- |
| Programming | Python, R, EViews |
| ML libraries | scikit-learn, XGBoost, SHAP, pandas, NumPy |
| Econometric methods | NARDL, VAR, Diebold–Yilmaz connectedness, Markov-switching, bounds testing, Granger causality, Cholesky FEVD, IRF, ADF, HAC / Newey–West, block bootstrap |
| ML methods | Random Forest, Gradient Boosting, Logistic Regression, walk-forward validation, SMOTE |
| Financial platforms | Bloomberg Terminal, TradingView |

---

## Contact

**Kunal Arora**

MBA in Business Analytics, Delhi School of Economics, University of Delhi · FRM Part I

[LinkedIn](https://www.linkedin.com/in/kunalar/) · <workkunal1@gmail.com>
