# 03 · Asymmetric Iron-Ore Shock Transmission in Australian Mining Equities

**A Nonlinear ARDL and Regime-Switching Analysis**

A study of whether iron-ore price shocks transmit *asymmetrically* into the equity returns of Australia's three largest listed producers — BHP, Rio Tinto, and Fortescue — over **January 2011 – May 2026 (N = 3,013 trading days)**, and whether that asymmetry is a stable structural feature of the relationship or a creature of crisis periods.

📄 **Preprint:** [SSRN — Abstract ID 6873399](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=6873399) <!-- verify this is the live link before publishing -->

`R` `NARDL` `Markov-Switching` `Bounds Testing` `Block Bootstrap` `HAC Inference`

---

## Research questions

1. **Long-run asymmetry** — Does iron-ore transmission differ between sustained price *increases* and *decreases*?
2. **Exposure gradient** — Does transmission scale with a firm's iron-ore concentration, so a pure-play producer responds more than diversified majors?
3. **Structural vs. state-dependent** — Is the asymmetry stable across calm and stressed markets, or does it change with the regime?

## Data

| Series | Source | Role |
| --- | --- | --- |
| BHP, Rio Tinto, Fortescue (adj. close) | Yahoo Finance | Equally-weighted portfolio — dependent variable |
| Iron ore 62% Fe CFR China (USD/t) | Investing.com CSV | Primary commodity regressor |
| Copper futures | Yahoo Finance | Secondary-commodity robustness check |
| AUD/USD, VIX, S&P 500 | Yahoo Finance | Controls (S&P 500 preferred over ASX 200 to avoid mechanical overlap) |

Equity and macro series are pulled programmatically via `quantmod`; the iron-ore 62% Fe CFR China series (Investing.com) is included as a CSV in `data/`.

## Methodology

- **NARDL** (Shin, Yu & Greenwood-Nimmo, 2014) — decomposes cumulative iron-ore log changes into positive (x⁺) and negative (x⁻) partial sums within a single error-correction structure; long-run multipliers L⁺ = −θ⁺/ρ, L⁻ = −θ⁻/ρ.
- **Bounds testing** (Pesaran, Shin & Smith, 2001) — F- and t-bounds run through `dynamac::pssbounds` with correct sample-specific critical values.
- **Markov-switching model** (Hamilton, 1989) — two-state (calm / stress) regime identification with switching intercept, AR coefficient, and variance.
- **Regime-conditional NARDL** — a stress dummy built from **filtered** Markov probabilities (no look-ahead bias) interacted with long- and short-run terms on the full continuous sample.
- **Inference** — Newey–West HAC standard errors throughout; a **two-stage moving-block bootstrap** that refits the entire pipeline (resample → refit Markov → rebuild dummy → re-estimate NARDL) to propagate first-stage regime-estimation uncertainty.
- **Diagnostics** — ADF unit-root tests, ARCH-LM, Breusch–Godfrey, Ramsey RESET, and an OLS-CUSUM parameter-stability test.

## Key findings

- **Robust long-run asymmetry.** A sustained iron-ore rise lifts the portfolio's long-run equilibrium value more than an equal fall lowers it — **L⁺ = 0.720 vs L⁻ = 0.523**, a 0.197 gap (Wald F = 13.90, p < 0.001), stable across market controls and the full (p, q) lag grid.
- **Concentration drives it.** Fortescue (pure-play) shows multipliers **three-to-four times** the diversified majors (L⁺ = 1.314 vs 0.389 BHP, 0.334 Rio Tinto) — the firms are distinct, not interchangeable, iron-ore instruments.
- **A two-part regime structure.** The long-run asymmetry is **regime-invariant** (calm gap 0.191 vs stress gap 0.198; difference test F = 0.24, p = 0.62) and survives the two-stage bootstrap — its regime-difference CI [−0.0012, 0.0010] is a well-identified zero. Meanwhile **short-run transmission intensity roughly triples in stress** in both directions.
- **Honest cointegration framing.** The bounds F (3.66) is inconclusive at 5%; the long-run reading rests on a significantly negative error-correction term (ρ = −0.016, t = −4.17; ~43-day half-life) and the bounds t-test, rather than a single marginal statistic.

**Implication:** hedges built on symmetric, constant exposure are systematically miscalibrated; hedge *levels* can rest on stable long-run parameters, while rebalancing *intensity* should rise in stressed markets.

## Repository contents

```
03_mining-asymmetry-nardl/
├── README.md
├── code/
│   └── mining_asymmetry_final_8.R      ← full pipeline (~1,100 lines, self-installing deps)
├── data/
│   └── Iron_ore_fines_62__Fe_CFR_Futures_Historical_Data.csv   ← iron ore 62% Fe CFR China (Investing.com), 2011–2026
└── outputs/
    ├── mining_asymmetry_SSRN.pdf       ← the working paper
    ├── tables/
    │   ├── 01_descriptive_statistics.csv     ← Table 1 (daily log-return moments)
    │   ├── 02_adf_unit_root_tests.csv        ← Table 2 (ADF, level & first diff.)
    │   ├── 03_regime_classification.csv      ← daily smoothed/filtered Markov probs.
    │   ├── 04_lag_robustness.csv             ← Table 7 (L⁺, L⁻, gap, Wald p, RESET p over p×q grid)
    │   ├── 05_bootstrap_classification.csv   ← classification-only bootstrap SEs/CIs
    │   └── 06_bootstrap_full_twostage.csv    ← Table 6 (two-stage bootstrap, propagates Markov uncertainty)
    └── figures/                              ← exported plots (iron-ore price, partial sums,
                                                long-run multipliers, firm-level, regime probability,
                                                regime-conditional panel, OLS-CUSUM)
```

## Reproducing

The script is self-contained — it checks for and installs any missing packages on first run.

```r
# from the project root, in R / RStudio
source("code/mining_asymmetry_final_8.R")
```

**Dependencies** (auto-installed): `quantmod`, `xts`, `zoo`, `tseries`, `FinTS`, `dynlm`, `car`, `lmtest`, `sandwich`, `MSwM`, `PerformanceAnalytics`, `dynamac`, `strucchange`, `ggplot2`, `dplyr`, `tidyr`, `scales`, `patchwork`, `lubridate`.

**Note on data:** equities, copper, FX, VIX and the S&P 500 download automatically via `quantmod`. The iron-ore 62% Fe CFR China series is provided in `data/`; point the script's iron-ore read path at that file.

## Citation

> Arora, K. (2026). *Asymmetric Iron-Ore Shock Transmission in Australian Mining Equities: A Nonlinear ARDL and Regime-Switching Analysis.* SSRN Working Paper.

---

*This paper has not been peer reviewed. Comments and suggestions are welcome.*
