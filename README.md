# Financial-Research

Empirical research in quantitative finance — ML ensemble forecasting, time series econometrics, and financial risk modelling.

**Python · EViews · scikit-learn · XGBoost · SHAP**

---

## About

This repository collects my research in empirical finance, spanning machine learning-based directional forecasting and econometric spillover modelling. Each project is self-contained with its own data, code, and outputs.

I am an MBA student in Business Analytics at the Delhi School of Economics, University of Delhi, with research interests in financial risk modelling, statistical learning under distribution shift, market volatility dynamics, and risk-aware portfolio construction.

---

## Projects

### [01 · Thematic Equity ML Forecasting](./01_thematic-equity-ml-forecasting/)

**Directional Prediction of Indian Thematic Equity Indices Using Machine Learning Ensemble Classifiers: A Multi-Horizon Analysis**

Evaluates the directional price-movement predictability of five Nifty thematic indices — Housing, Consumption, PSE, Transportation & Logistics, and ESG — using Random Forest and Gradient Boosting ensemble classifiers across five forecast horizons (1–20 trading days). Features a 41-indicator feature engineering pipeline, two-stage feature selection, SHAP-based interpretability, walk-forward validation, and Kelly Criterion position sizing.

📄 **Preprint:** [SSRN — Abstract ID 6748042](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=6748042)

`Python` `scikit-learn` `XGBoost` `SHAP` `Random Forest` `Gradient Boosting`

---

### [02 · Oil Market Liquidity & Volatility Spillovers](./02_oil-market-spillovers/)

**Liquidity and Volatility Spillovers in the Global Oil Market: A Time Series Analytics Study**

Examines price return, realized volatility, and trading volume spillovers between WTI and Brent crude oil futures over 10 Indian financial years (April 2016 – March 2026) using a six-variable VAR(4) system and the Diebold-Yilmaz (2012) connectedness framework. Covers three major global crises: COVID-19, the Russia-Ukraine war, and the 2023 Iran-related Middle East conflict.

`Python` `EViews` `VAR` `Granger Causality` `Cholesky FEVD` `IRF`

---

## Repository Structure

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
└── 02_oil-market-spillovers/
    ├── README.md
    ├── data/                  ← processed CSV
    ├── notebooks/             ← Jupyter notebook
    └── outputs/
        ├── figures/           ← price/return/volatility plots
        ├── tables/            ← ADF tests, VAR output, Granger, FEVD, IRF (PDFs)
        └── TimeSeriesProject.pdf
```

---

## Skills & Tools

| Category | Tools |
|---|---|
| Programming | Python, EViews |
| ML Libraries | scikit-learn, XGBoost, SHAP, pandas, NumPy |
| Econometric Methods | VAR, Granger Causality, Cholesky FEVD, IRF, ADF |
| ML Methods | Random Forest, Gradient Boosting, Logistic Regression, Walk-Forward Validation, SMOTE |
| Financial Platforms | Bloomberg Terminal, TradingView |

---

## Contact

**Kunal Arora**

MBA in Business Analytics, Delhi School of Economics, University of Delhi · FRM Part I 

[LinkedIn](https://www.linkedin.com/in/kunalar/) · workkunal1@gmail.com
