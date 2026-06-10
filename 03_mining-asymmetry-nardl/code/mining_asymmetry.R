# =============================================================================
# ASYMMETRIC COMMODITY SHOCK TRANSMISSION IN AUSTRALIAN MINING EQUITIES
# A Nonlinear ARDL and Regime-Switching Analysis  —  VERSION 8 (re-run build)
#
# Author : Kunal Arora
#
# WHAT CHANGED FROM v7  (all marked inline with  # === v8 CHANGE ===):
#   (1) SPEC CONSISTENCY: copper, market-control robustness, and firm-level
#       NARDLs now use the SAME BIC-selected (p,q) lags as the main model,
#       via one shared formula builder. (Fixes "same specification" claim.)
#   (2) BOUNDS TEST: now reports BOTH the F-bounds and t-bounds tests and runs
#       them through dynamac::pssbounds (correct PSS/Narayan CVs, auto-selected
#       for the sample size). Replaces the lone, marginal HAC-F comparison.
#   (3) LAG-ROBUSTNESS: full p x q grid re-estimated; L+, L-, gap, Wald p and
#       RESET p reported for every spec (supports the RESET defence).
#   (4) GENERATED-REGRESSOR BOOTSTRAP: regime-conditional results now carry
#       bootstrap SEs that propagate first-stage (Markov) uncertainty.
#   (5) SAMPLE HARMONISATION: the Markov model is now fit on the SAME merged
#       sample as the NARDL, so N is consistent across all tables.
#   (6) REGIME LABELS: outputs use "Stress (high-vol)" / "Calm (low-vol)"
#       instead of "Regime 1/2" to remove the numbering ambiguity.
#   (7) CUSUM is run on the EXACT main NARDL model (incl. q lagged differences).
#
# Data    : Yahoo Finance (equities, copper, FX, VIX, S&P 500)
#           + Investing.com CSV (iron ore 62% Fe CFR China)
# =============================================================================


# =============================================================================
# SECTION 0 — PACKAGE SETUP
# =============================================================================
required_packages <- c(
  # --- Data ---
  "quantmod", "xts", "zoo",
  # --- Analysis ---
  "tseries", "FinTS", "dynlm", "car", "lmtest", "sandwich", "MSwM",
  "PerformanceAnalytics",
  "dynamac",      # === v8 CHANGE (2): pssbounds() — correct bounds-test CVs
  "strucchange",  # === v8 CHANGE (7): efp()/sctest() promoted to setup
  # --- Visualisation ---
  "ggplot2", "dplyr", "tidyr", "scales", "patchwork", "lubridate"
)

new_pkg <- required_packages[!(required_packages %in%
                                 installed.packages()[, "Package"])]
if (length(new_pkg) > 0) {
  message("Installing missing packages: ", paste(new_pkg, collapse = ", "))
  install.packages(new_pkg, dependencies = TRUE)
}
invisible(lapply(required_packages, library, character.only = TRUE))
message("All packages loaded.\n")


# =============================================================================
# SECTION 1 — DATA DOWNLOAD
# =============================================================================
start_date <- "2011-01-01"   # post-spot-pricing era (annual benchmarks ended 2010)
end_date   <- Sys.Date()

message("Downloading data from Yahoo Finance...")
message("Date range: ", start_date, " to ", format(end_date, "%Y-%m-%d"), "\n")

# 1.1  Australian Mining Equities
bhp_raw <- getSymbols("BHP.AX", src = "yahoo", from = start_date, to = end_date, auto.assign = FALSE)
rio_raw <- getSymbols("RIO.AX", src = "yahoo", from = start_date, to = end_date, auto.assign = FALSE)
fmg_raw <- getSymbols("FMG.AX", src = "yahoo", from = start_date, to = end_date, auto.assign = FALSE)
message("Equity data downloaded: BHP.AX, RIO.AX, FMG.AX")

# 1.2  Iron Ore (Investing.com CSV — clean weekday-only series)
io_csv_file <- "Iron_ore_fines_62__Fe_CFR_Futures_Historical_Data.csv"
if (!file.exists(io_csv_file)) {
  stop(paste0("Iron ore CSV not found in working directory: ", getwd(), "\n",
              "Download from: https://www.investing.com/commodities/",
              "iron-ore-62-cfr-futures-historical-data\nSave as: ", io_csv_file))
}
io_csv_raw <- read.csv(io_csv_file, stringsAsFactors = FALSE)
io_csv_raw$Date  <- as.Date(io_csv_raw$Date, format = "%m/%d/%Y")
io_csv_raw$Price <- as.numeric(gsub(",", "", io_csv_raw$Price))
io_csv_raw       <- io_csv_raw[order(io_csv_raw$Date), ]
io_csv_raw       <- io_csv_raw[io_csv_raw$Date >= as.Date(start_date) &
                                 io_csv_raw$Date <= end_date, ]
iron_ore_xts <- xts(io_csv_raw$Price, order.by = io_csv_raw$Date)
colnames(iron_ore_xts) <- "IronOre"
message(sprintf("Iron ore CSV loaded: %d obs | %s to %s | NAs: %d",
                nrow(io_csv_raw), format(min(io_csv_raw$Date), "%Y-%m-%d"),
                format(max(io_csv_raw$Date), "%Y-%m-%d"), sum(is.na(iron_ore_xts))))

copper_raw <- getSymbols("HG=F",     src = "yahoo", from = start_date, to = end_date, auto.assign = FALSE)
audusd_raw <- getSymbols("AUDUSD=X", src = "yahoo", from = start_date, to = end_date, auto.assign = FALSE)
vix_raw    <- getSymbols("^VIX",     src = "yahoo", from = start_date, to = end_date, auto.assign = FALSE)
spx_raw    <- getSymbols("^GSPC",    src = "yahoo", from = start_date, to = end_date, auto.assign = FALSE)
message("Copper, AUD/USD, VIX, S&P 500 downloaded.\n")


# =============================================================================
# SECTION 2 — DATA CONSTRUCTION
# =============================================================================
# 2.1  Adjusted prices (equities) / closing prices (rest)
bhp_p  <- Ad(bhp_raw); rio_p <- Ad(rio_raw); fmg_p <- Ad(fmg_raw)
io_p   <- iron_ore_xts
cu_p   <- Cl(copper_raw); fx_p <- Cl(audusd_raw); vix_p <- Cl(vix_raw); spx_p <- Cl(spx_raw)
colnames(bhp_p) <- "BHP"; colnames(rio_p) <- "RIO"; colnames(fmg_p) <- "FMG"
colnames(io_p)  <- "IronOre"; colnames(cu_p) <- "Copper"; colnames(fx_p) <- "AUDUSD"
colnames(vix_p) <- "VIX"; colnames(spx_p) <- "SP500"

# Interior NA interpolation (disclose in paper: creates a few near-zero returns)
bhp_p <- na.approx(bhp_p, na.rm = FALSE); rio_p <- na.approx(rio_p, na.rm = FALSE)
fmg_p <- na.approx(fmg_p, na.rm = FALSE); io_p  <- na.approx(io_p,  na.rm = FALSE)
cu_p  <- na.approx(cu_p,  na.rm = FALSE); fx_p  <- na.approx(fx_p,  na.rm = FALSE)
vix_p <- na.approx(vix_p, na.rm = FALSE); spx_p <- na.approx(spx_p, na.rm = FALSE)
message("Prices extracted; interior NAs interpolated.")

# 2.2  Log prices
ln_bhp <- log(bhp_p); ln_rio <- log(rio_p); ln_fmg <- log(fmg_p)
ln_io  <- log(io_p);  ln_cu  <- log(cu_p);  ln_fx  <- log(fx_p)
ln_vix <- log(vix_p); ln_spx <- log(spx_p)

# 2.3  Equally weighted portfolio (avg of log prices => avg of log returns)
equity_ln <- na.omit(merge(ln_bhp, ln_rio, ln_fmg, join = "inner"))
colnames(equity_ln) <- c("lnBHP", "lnRIO", "lnFMG")
portfolio_ln <- xts(rowMeans(equity_ln), order.by = index(equity_ln))
colnames(portfolio_ln) <- "lnPortfolio"
portfolio_ret <- na.omit(diff(portfolio_ln))          # full equity-sample returns (visual only)
colnames(portfolio_ret) <- "PortfolioRet"
message("Portfolio constructed (equally weighted BHP/RIO/FMG).")

# 2.4  Merge everything on common dates
all_ln <- na.omit(merge(portfolio_ln, ln_io, ln_cu, ln_fx, ln_vix, ln_spx, join = "inner"))
colnames(all_ln) <- c("lnPortfolio", "lnIronOre", "lnCopper", "lnAUDUSD", "lnVIX", "lnSP500")
all_ret <- na.omit(diff(all_ln))
colnames(all_ret) <- c("Portfolio", "IronOre", "Copper", "AUDUSD", "VIX", "SP500")

# === v8 CHANGE (5): portfolio returns ON THE MERGED SAMPLE, used for the
# Markov model so the regime sample == the NARDL sample (one consistent N).
portfolio_ret_merged <- na.omit(diff(all_ln$lnPortfolio))
colnames(portfolio_ret_merged) <- "PortfolioRet"

message(sprintf("Final merged dataset: %d obs | %s to %s\n",
                nrow(all_ln), format(start(all_ln), "%Y-%m-%d"),
                format(end(all_ln), "%Y-%m-%d")))


# =============================================================================
# SECTION 3 — PRELIMINARY ANALYSIS
# =============================================================================
# 3.1  Descriptive statistics
message("=== DESCRIPTIVE STATISTICS (Log Returns) ===")
ret_df <- as.data.frame(all_ret)
desc_stats <- data.frame(
  Series      = names(ret_df),
  N           = sapply(ret_df, function(x) sum(!is.na(x))),
  Mean_pct    = sapply(ret_df, mean, na.rm = TRUE) * 100,
  StdDev_pct  = sapply(ret_df, sd,   na.rm = TRUE) * 100,
  Min_pct     = sapply(ret_df, min,  na.rm = TRUE) * 100,
  Max_pct     = sapply(ret_df, max,  na.rm = TRUE) * 100,
  Skewness    = sapply(ret_df, skewness, na.rm = TRUE),
  Ex_Kurtosis = sapply(ret_df, kurtosis, na.rm = TRUE),
  row.names = NULL)
print(desc_stats, digits = 4)

# 3.2  ADF unit root tests
message("\n=== ADF UNIT ROOT TESTS ===")
adf_df <- as.data.frame(all_ln); adf_vars <- names(adf_df)
adf_results <- do.call(rbind, lapply(adf_vars, function(v) {
  lev <- tseries::adf.test(na.omit(adf_df[[v]]))
  dif <- tseries::adf.test(na.omit(diff(adf_df[[v]])))
  order <- if (lev$p.value < 0.05) "I(0)" else
    if (dif$p.value < 0.05) "I(1)" else "I(2) CHECK"
  data.frame(Variable = v, ADF_Level = round(lev$statistic, 3),
             p_Level = round(lev$p.value, 3), ADF_Diff = round(dif$statistic, 3),
             p_Diff = round(dif$p.value, 3), Order = order, stringsAsFactors = FALSE)
}))
print(adf_results, row.names = FALSE)
if (any(grepl("I\\(2\\)", adf_results$Order)))
  warning("One or more variables appear I(2). Review before NARDL.") else
    message("All variables I(0) or I(1). NARDL valid.")

# 3.3  ARCH-LM test (justifies the regime model)  — on the MERGED-sample returns
message("\n=== ARCH-LM TEST (Portfolio Returns, merged sample) ===")
port_vec_merged <- as.numeric(portfolio_ret_merged)
arch_res <- FinTS::ArchTest(port_vec_merged, lags = 10)
message(sprintf("  Chi-sq(10) = %.4f | p = %.6f", arch_res$statistic, arch_res$p.value))
message(if (arch_res$p.value < 0.05) "  -> ARCH effects confirmed; regime model justified."
        else "  -> No significant ARCH effects.")

# 3.4  Correlations
message("\n=== PAIRWISE CORRELATIONS (Log Returns) ===")
print(round(cor(as.matrix(all_ret), use = "complete.obs"), 3))


# =============================================================================
# SECTION 4 — EXPLORATORY VISUALISATIONS
# =============================================================================
theme_paper <- theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "grey45"),
        axis.title = element_text(size = 9),
        panel.grid.minor = element_blank(), legend.position = "bottom")
theme_set(theme_paper)

events <- data.frame(
  date = as.Date(c("2016-01-12","2019-01-25","2020-03-23","2021-05-10","2021-09-20",
                   "2022-02-24","2023-07-01","2025-04-02","2025-10-01")),
  label = c("Iron Ore\nTrough","Vale Dam\nCollapse","COVID\nCrash","Iron Ore\nPeak",
            "Evergrande\nCrisis","Russia-\nUkraine","Safeguard\nMech.","US Tariff\nShock",
            "China-BHP\nStandoff"),
  y_pos = c(0.90,0.75,0.90,0.90,0.75,0.90,0.75,0.90,0.75))

# 4.1  Portfolio returns
ret_plot_df <- data.frame(Date = index(portfolio_ret), Return = as.numeric(portfolio_ret)*100)
p_returns <- ggplot(ret_plot_df, aes(Date, Return)) +
  geom_line(color = "steelblue", linewidth = 0.35, alpha = 0.85) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_vline(data = events, aes(xintercept = date), linetype = "dashed",
             color = "grey35", linewidth = 0.55) +
  geom_text(data = events, aes(x = date + 30,
            y = max(ret_plot_df$Return, na.rm = TRUE) * y_pos, label = label),
            size = 2.5, hjust = 0, color = "grey25") +
  labs(title = "Daily Log Returns — BHP / RIO / FMG Equally Weighted Portfolio",
       subtitle = "ASX-listed Australian iron ore miners (2011 – present)",
       x = NULL, y = "Log Return (%)") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y")
print(p_returns)

# 4.2  Iron ore price
io_plot_df <- data.frame(Date = index(io_p), Price = as.numeric(io_p))
p_ironore <- ggplot(io_plot_df, aes(Date, Price)) +
  geom_line(color = "firebrick3", linewidth = 0.5) +
  geom_vline(data = events, aes(xintercept = date), linetype = "dashed",
             color = "grey35", linewidth = 0.55) +
  geom_text(data = events, aes(x = date + 30,
            y = max(io_plot_df$Price, na.rm = TRUE) * y_pos, label = label),
            size = 2.5, hjust = 0, color = "grey25") +
  labs(title = "Iron Ore Price — IOF 62% Fe CFR China",
       subtitle = "USD per dry metric tonne (2011 – present)", x = NULL, y = "Price (USD / t)") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y")
print(p_ironore)

# 4.3  Rolling correlation (motivates the NARDL)
io_ret_xts <- na.omit(diff(io_p))
aligned_rets <- na.omit(merge(portfolio_ret, io_ret_xts, join = "inner"))
colnames(aligned_rets) <- c("Portfolio", "IronOre")
roll_cor <- rollapply(aligned_rets, width = 60,
                      FUN = function(x) cor(x[,1], x[,2], use = "complete.obs"),
                      by.column = FALSE, align = "right")
roll_cor_df <- data.frame(Date = index(roll_cor), Correlation = as.numeric(roll_cor))
p_rollcor <- ggplot(roll_cor_df, aes(Date, Correlation)) +
  geom_line(color = "darkorange", linewidth = 0.5) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_hline(yintercept = mean(roll_cor_df$Correlation, na.rm = TRUE),
             linetype = "dotted", color = "grey50") +
  geom_vline(data = events, aes(xintercept = date), linetype = "dashed",
             color = "grey35", linewidth = 0.5) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(title = "60-Day Rolling Correlation: Portfolio vs Iron Ore Returns",
       subtitle = "Dotted line = full-sample average", x = NULL, y = "Pearson Correlation") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y")
print(p_rollcor)


# =============================================================================
# SECTION 5 — PARTIAL SUM DECOMPOSITION
# =============================================================================
message("=== PARTIAL SUM DECOMPOSITION ===")
ln_df <- as.data.frame(all_ln); ln_df$Date <- as.Date(rownames(ln_df)); rownames(ln_df) <- NULL
d_io <- c(NA, diff(ln_df$lnIronOre)); d_cu <- c(NA, diff(ln_df$lnCopper))
ln_df$POS_IO <- cumsum(ifelse(is.na(d_io), 0, pmax(d_io, 0)))
ln_df$NEG_IO <- cumsum(ifelse(is.na(d_io), 0, pmin(d_io, 0)))
ln_df$POS_CU <- cumsum(ifelse(is.na(d_cu), 0, pmax(d_cu, 0)))
ln_df$NEG_CU <- cumsum(ifelse(is.na(d_cu), 0, pmin(d_cu, 0)))
nardl_data <- ln_df[!is.na(d_io), ]
nardl_data <- nardl_data[complete.cases(nardl_data), ]
message(sprintf("NARDL dataset: %d obs | %s to %s\n",
                nrow(nardl_data), min(nardl_data$Date), max(nardl_data$Date)))

decomp_long <- tidyr::pivot_longer(
  data.frame(Date = nardl_data$Date, Positive = nardl_data$POS_IO, Negative = nardl_data$NEG_IO),
  -Date, names_to = "Type", values_to = "Value")
p_decomp <- ggplot(decomp_long, aes(Date, Value, color = Type)) +
  geom_line(linewidth = 0.5) +
  scale_color_manual(values = c("Positive" = "darkgreen", "Negative" = "firebrick"),
                     labels = c("Positive" = "Cumulative positive (x+)",
                                "Negative" = "Cumulative negative (x-)"), name = "") +
  geom_vline(data = events, aes(xintercept = date), linetype = "dashed",
             color = "grey40", linewidth = 0.5, inherit.aes = FALSE) +
  labs(title = "Partial Sum Decomposition of Iron Ore Log Returns",
       subtitle = "x+ accumulates price rises | x- accumulates price falls",
       x = NULL, y = "Cumulative partial sum") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y")
print(p_decomp)


# =============================================================================
# === v8 CHANGE (1): SHARED NARDL FORMULA BUILDER
# -----------------------------------------------------------------------------
# One builder used by the main model, copper, the market-control robustness,
# and the firm-level models — so every specification carries the SAME (p,q)
# lag structure. This is what makes the "same specification" claim true.
#   dep      : name of the dependent log-price level (e.g. "lnPortfolio")
#   pos,neg  : names of the positive / negative partial-sum regressors
#   controls : character vector of long-run control levels (may be empty)
#   p        : number of lagged Δy (own short-run) terms
#   q        : number of lagged Δx+ / Δx- (commodity short-run) terms
# =============================================================================
build_nardl_formula <- function(dep, pos, neg, controls, p, q) {
  ar    <- paste(sprintf("L(d(%s), %d)", dep, 1:p), collapse = " + ")
  ctrlL <- if (length(controls) > 0)
    paste(sprintf("L(%s, 1)", controls), collapse = " + ") else NULL
  qpos  <- if (q > 0) paste(sprintf("L(d(%s), %d)", pos, 1:q), collapse = " + ") else NULL
  qneg  <- if (q > 0) paste(sprintf("L(d(%s), %d)", neg, 1:q), collapse = " + ") else NULL
  rhs <- c(sprintf("L(%s, 1)", dep), sprintf("L(%s, 1)", pos), sprintf("L(%s, 1)", neg),
           ctrlL, ar, sprintf("d(%s)", pos), sprintf("d(%s)", neg), qpos, qneg)
  rhs <- rhs[!is.null(rhs) & nzchar(rhs)]
  as.formula(paste0("d(", dep, ") ~ ", paste(rhs, collapse = " + ")))
}


# =============================================================================
# SECTION 6 — NARDL ESTIMATION: IRON ORE (MAIN SPECIFICATION)
# =============================================================================
message("=== NARDL ESTIMATION: IRON ORE (MAIN SPECIFICATION) ===")

# --- 6.0  BIC lag selection over p in {2..5} x q in {0,1} -------------------
nardl_ts_full <- ts(
  nardl_data[, c("lnPortfolio","POS_IO","NEG_IO","lnAUDUSD","lnVIX","lnSP500")],
  start = c(as.integer(format(min(nardl_data$Date), "%Y")),
            as.integer(format(min(nardl_data$Date), "%j"))), frequency = 252)

ctrl_main <- c("lnAUDUSD", "lnVIX", "lnSP500")
lag_grid  <- expand.grid(p = 2:5, q = 0:1)

lag_results <- do.call(rbind, lapply(seq_len(nrow(lag_grid)), function(i) {
  f <- build_nardl_formula("lnPortfolio","POS_IO","NEG_IO", ctrl_main,
                           lag_grid$p[i], lag_grid$q[i])
  m <- tryCatch(dynlm(f, data = nardl_ts_full), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  data.frame(p = lag_grid$p[i], q = lag_grid$q[i],
             K = length(coef(m)), N = length(residuals(m)),
             AIC = round(AIC(m), 2), BIC = round(BIC(m), 2))
}))
print(lag_results, row.names = FALSE)
best_bic <- lag_results[which.min(lag_results$BIC), ]
p_sel <- best_bic$p; q_sel <- best_bic$q
message(sprintf("BIC-selected specification: p = %d, q = %d  (used everywhere)\n", p_sel, q_sel))

# --- 6.1  Main NARDL --------------------------------------------------------
formula_main <- build_nardl_formula("lnPortfolio","POS_IO","NEG_IO", ctrl_main, p_sel, q_sel)
nardl_iron   <- dynlm(formula_main, data = nardl_ts_full)
message("--- Main NARDL (Iron Ore) ---"); print(summary(nardl_iron))
nardl_iron_hac <- coeftest(nardl_iron, vcov = NeweyWest(nardl_iron, lag = 4))
message("\n--- HAC-Robust Coefficients ---"); print(nardl_iron_hac)

# --- 6.2  Bounds test: F-bounds AND t-bounds via dynamac::pssbounds --------
# === v8 CHANGE (2) ===
# PSS (2001) give two cointegration tests: the joint F on all level terms and
# the t on the lagged dependent (rho). We report BOTH and let pssbounds() apply
# the correct critical-value table for the sample size (asymptotic PSS at
# large N; Narayan small-sample tables only where they apply). Note: the bounds
# F/t are computed WITHOUT HAC, as the PSS/Narayan CVs are derived for the
# standard statistics; the HAC versions are reported alongside for transparency.
message("\n--- Bounds Test for Cointegration (F-bounds + t-bounds) ---")
bounds_hyp <- c("L(lnPortfolio, 1) = 0", "L(POS_IO, 1) = 0", "L(NEG_IO, 1) = 0",
                "L(lnAUDUSD, 1) = 0", "L(lnVIX, 1) = 0", "L(lnSP500, 1) = 0")

# Standard (non-HAC) F for the bounds comparison
bounds_F_std <- linearHypothesis(nardl_iron, bounds_hyp)
f_stat_std   <- bounds_F_std$F[2]
# HAC F (reported for transparency only)
bounds_F_hac <- linearHypothesis(nardl_iron, bounds_hyp, vcov = NeweyWest(nardl_iron, lag = 4))
f_stat_hac   <- bounds_F_hac$F[2]

# t-statistic on rho (lagged dependent): standard and HAC
sm        <- summary(nardl_iron)$coefficients
t_rho_std <- sm["L(lnPortfolio, 1)", "t value"]
t_rho_hac <- nardl_iron_hac["L(lnPortfolio, 1)", "z value"]
N_bounds  <- length(residuals(nardl_iron))

message(sprintf("  F-bounds (standard): %.4f   |  F-bounds (HAC, info only): %.4f",
                f_stat_std, f_stat_hac))
message(sprintf("  t-bounds rho (standard): %.4f  |  t (HAC, info only): %.4f",
                t_rho_std, t_rho_hac))

# k = 5 long-run forcing variables: POS_IO, NEG_IO, AUDUSD, VIX, SP500
# Case III = unrestricted intercept, no trend
message("\n--- dynamac::pssbounds (Case III, k = 5) ---")
pss <- tryCatch(
  dynamac::pssbounds(obs = N_bounds, fstat = f_stat_std, tstat = t_rho_std,
                     case = 3, k = 5),
  error = function(e) { message("pssbounds error: ", e$message); NULL })
# pssbounds() prints the F-test and t-test verdicts against the correct bounds.

# --- 6.3  Wald test for long-run asymmetry (core hypothesis) ---------------
message("\n--- Wald Test: Long-Run Asymmetry (H0: theta+ = theta-) ---")
asym_test_iron <- linearHypothesis(nardl_iron, "L(POS_IO, 1) = L(NEG_IO, 1)",
                                   vcov = NeweyWest(nardl_iron, lag = 4))
print(asym_test_iron)

# --- 6.4  Long-run multipliers ---------------------------------------------
coefs <- coef(nardl_iron)
rho   <- coefs["L(lnPortfolio, 1)"]
L_pos <- -coefs["L(POS_IO, 1)"] / rho
L_neg <- -coefs["L(NEG_IO, 1)"] / rho
message(sprintf("\nrho = %.4f | L+ = %.4f | L- = %.4f | (L+ - L-) = %.4f",
                rho, L_pos, L_neg, L_pos - L_neg))

# --- 6.5  Residual diagnostics on the EXACT main model ---------------------
message("\n--- Residual Diagnostics (main NARDL) ---")
bg  <- bgtest(nardl_iron, order = 4)
rst <- resettest(nardl_iron, power = 2:3)
message(sprintf("  Breusch-Godfrey LM(4): LM = %.3f, p = %.4f", bg$statistic, bg$p.value))
message(sprintf("  Ramsey RESET (powers 2:3): F = %.3f, p = %.4f", rst$statistic, rst$p.value))


# =============================================================================
# === v8 CHANGE (3): LAG-ROBUSTNESS TABLE (supports the RESET defence)
# -----------------------------------------------------------------------------
# Re-estimate the FULL p x q grid and record L+, L-, gap, asymmetry-Wald p, and
# RESET p for each. If L+/L-/gap/Wald are stable across specs while RESET keeps
# rejecting, then whatever RESET reacts to does NOT contaminate the reported
# multipliers — exactly the argument we make in the paper.
# =============================================================================
message("\n=== LAG-ROBUSTNESS GRID (RESET defence) ===")
lag_robust <- do.call(rbind, lapply(seq_len(nrow(lag_grid)), function(i) {
  pp <- lag_grid$p[i]; qq <- lag_grid$q[i]
  f  <- build_nardl_formula("lnPortfolio","POS_IO","NEG_IO", ctrl_main, pp, qq)
  m  <- tryCatch(dynlm(f, data = nardl_ts_full), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  V  <- NeweyWest(m, lag = 4); b <- coef(m); r <- b["L(lnPortfolio, 1)"]
  w  <- linearHypothesis(m, "L(POS_IO, 1) = L(NEG_IO, 1)", vcov = V)
  rs <- tryCatch(resettest(m, power = 2:3)$p.value, error = function(e) NA)
  data.frame(p = pp, q = qq,
             Lpos = round(-b["L(POS_IO, 1)"]/r, 4),
             Lneg = round(-b["L(NEG_IO, 1)"]/r, 4),
             gap  = round((-b["L(POS_IO, 1)"]/r) - (-b["L(NEG_IO, 1)"]/r), 4),
             wald_p  = round(w$Pr[2], 4),
             reset_p = round(rs, 4),
             BIC = round(BIC(m), 1), row.names = NULL)
}))
print(lag_robust, row.names = FALSE)
message("Read: L+/L-/gap/Wald p stable across the grid => inference is invariant")
message("to the short-run dynamics RESET might be reacting to.\n")


# =============================================================================
# SECTION 7 — NARDL ESTIMATION: COPPER (ROBUSTNESS)   [now BIC-consistent]
# =============================================================================
# === v8 CHANGE (1): copper uses the SAME (p_sel, q_sel) as the main model ===
message("=== NARDL ESTIMATION: COPPER (ROBUSTNESS, BIC-consistent) ===")
nardl_ts_cu <- ts(
  nardl_data[, c("lnPortfolio","POS_CU","NEG_CU","lnAUDUSD","lnVIX","lnSP500")],
  start = c(as.integer(format(min(nardl_data$Date), "%Y")),
            as.integer(format(min(nardl_data$Date), "%j"))), frequency = 252)
formula_cu   <- build_nardl_formula("lnPortfolio","POS_CU","NEG_CU", ctrl_main, p_sel, q_sel)
nardl_copper <- dynlm(formula_cu, data = nardl_ts_cu)
message("--- HAC-Robust Coefficients (Copper) ---")
print(coeftest(nardl_copper, vcov = NeweyWest(nardl_copper, lag = 4)))
asym_test_cu <- linearHypothesis(nardl_copper, "L(POS_CU, 1) = L(NEG_CU, 1)",
                                 vcov = NeweyWest(nardl_copper, lag = 4))
message("--- Wald: Copper asymmetry ---"); print(asym_test_cu)
coefs_cu <- coef(nardl_copper); rho_cu <- coefs_cu["L(lnPortfolio, 1)"]
L_pos_cu <- -coefs_cu["L(POS_CU, 1)"] / rho_cu
L_neg_cu <- -coefs_cu["L(NEG_CU, 1)"] / rho_cu
message(sprintf("Copper L+ = %.4f | L- = %.4f  (corroborating evidence; check whether",
                L_pos_cu, L_neg_cu))
message("individual long-run coefficients are significant before claiming cointegration.)\n")


# =============================================================================
# SECTION 7B — MARKET-CONTROL ROBUSTNESS (THREE SPECS)  [now BIC-consistent]
# =============================================================================
message("=== SECTION 7B: MARKET-CONTROL ROBUSTNESS (BIC-consistent) ===")

# === v8 CHANGE (1): control-robustness models now use the shared builder ===
run_nardl_control <- function(control_vars, control_label, data_df, p, q) {
  message(sprintf("\n--- %s ---", control_label))
  core <- c("lnPortfolio", "POS_IO", "NEG_IO", "lnAUDUSD", "lnVIX")
  ts_data <- ts(data_df[, c(core, control_vars)],
                start = c(as.integer(format(min(data_df$Date), "%Y")),
                          as.integer(format(min(data_df$Date), "%j"))), frequency = 252)
  controls <- c("lnAUDUSD", "lnVIX", control_vars)
  f <- build_nardl_formula("lnPortfolio","POS_IO","NEG_IO", controls, p, q)
  mod <- tryCatch(dynlm(f, data = ts_data),
                  error = function(e) { message("  Error: ", e$message); NULL })
  if (is.null(mod)) return(NULL)
  V <- NeweyWest(mod, lag = 4); b <- coef(mod); r <- b["L(lnPortfolio, 1)"]
  asym <- tryCatch(linearHypothesis(mod, "L(POS_IO, 1) = L(NEG_IO, 1)", vcov = V),
                   error = function(e) NULL)
  Lpos <- -b["L(POS_IO, 1)"]/r; Lneg <- -b["L(NEG_IO, 1)"]/r
  p_val <- if (!is.null(asym)) asym$Pr[2] else NA
  # bounds F (standard) for this spec
  bh <- c("L(lnPortfolio, 1) = 0","L(POS_IO, 1) = 0","L(NEG_IO, 1) = 0",
          "L(lnAUDUSD, 1) = 0","L(lnVIX, 1) = 0")
  if (length(control_vars) > 0) bh <- c(bh, sprintf("L(%s, 1) = 0", control_vars))
  fb <- tryCatch(linearHypothesis(mod, bh)$F[2], error = function(e) NA)
  message(sprintf("  rho = %.4f | L+ = %.4f | L- = %.4f | Wald p = %.4f | bounds F = %.3f",
                  r, Lpos, Lneg, p_val, fb))
  list(label = control_label, Lpos = Lpos, Lneg = Lneg, asym_p = p_val,
       rho = r, boundsF = fb)
}

# ASX200 (for the endogeneity-contaminated comparison only)
asx_raw <- tryCatch(getSymbols("^AXJO", src = "yahoo", from = start_date,
                               to = end_date, auto.assign = FALSE),
                    error = function(e) { message("ASX200 unavailable: ", e$message); NULL })
if (!is.null(asx_raw)) {
  ln_asx <- log(na.approx(Cl(asx_raw), na.rm = FALSE))
  asx_df <- data.frame(Date = as.Date(index(ln_asx)), lnASX200 = as.numeric(ln_asx))
  nardl_data_m3 <- merge(nardl_data, asx_df, by = "Date", all.x = TRUE)
  nardl_data_m3$lnASX200[is.na(nardl_data_m3$lnASX200)] <-
    na.approx(nardl_data_m3$lnASX200, na.rm = FALSE)
  nardl_data_m3 <- na.omit(nardl_data_m3)
} else nardl_data_m3 <- NULL

m1 <- run_nardl_control(character(0), "Model 1: No market index (AUD/USD + VIX)",
                        nardl_data, p_sel, q_sel)
m2 <- run_nardl_control("lnSP500", "Model 2: S&P 500 (main control)",
                        nardl_data, p_sel, q_sel)
m3 <- if (!is.null(nardl_data_m3))
  run_nardl_control("lnASX200", "Model 3: ASX200 (endogeneity comparison only)",
                    nardl_data_m3, p_sel, q_sel) else NULL

message("\n=== ROBUSTNESS SUMMARY (asymmetry across market controls) ===")
message(sprintf("%-46s %-8s %-8s %-10s %-8s", "Specification", "L+", "L-", "L+ - L-", "Wald p"))
message(strrep("-", 84))
message(sprintf("%-46s %-8.4f %-8.4f %-10.4f %-8s", "Main (S&P 500, Section 6)",
                L_pos, L_neg, L_pos - L_neg, sprintf("%.4f", asym_test_iron$Pr[2])))
for (m in list(m1, m2, m3)) if (!is.null(m))
  message(sprintf("%-46s %-8.4f %-8.4f %-10.4f %-8.4f",
                  m$label, m$Lpos, m$Lneg, m$Lpos - m$Lneg, m$asym_p))
message(strrep("-", 84), "\n")


# =============================================================================
# SECTION 8 — MARKOV-SWITCHING REGIME IDENTIFICATION
# =============================================================================
# === v8 CHANGE (5): fit on portfolio_ret_merged so the regime sample equals
# the NARDL sample (one consistent N across Tables 1, regime table, and Sec 9).
# === v8 CHANGE (6): outputs labelled Stress/Calm, not Regime 1/2.
message("\n=== MARKOV-SWITCHING REGIME IDENTIFICATION (merged sample) ===")

port_vec <- as.numeric(portfolio_ret_merged)   # merged-sample returns
n        <- length(port_vec)
ms_df    <- data.frame(ret = port_vec[-1], lag1 = port_vec[-n])
lm_base  <- lm(ret ~ lag1, data = ms_df)

message("Fitting 2-state Markov-Switching model (30-60s)...")
ms_model <- tryCatch(
  MSwM::msmFit(lm_base, k = 2, sw = c(TRUE, TRUE, TRUE), p = 0,
               control = list(parallel = FALSE, maxiter = 1000, tol = 1e-8)),
  error = function(e) {
    message("Full switching failed (", e$message, "); trying intercept+variance only...")
    tryCatch(MSwM::msmFit(lm_base, k = 2, sw = c(TRUE, FALSE, TRUE), p = 0,
                          control = list(parallel = FALSE, maxiter = 1000, tol = 1e-6)),
             error = function(e2) { message("Both failed: ", e2$message); NULL }) })
if (is.null(ms_model)) stop("MS model could not be fitted.")
message("\n--- Markov-Switching Results ---"); summary(ms_model)

# 8.1  Extract / characterise regimes
smoothed_probs <- ms_model@Fit@smoProb
filtered_probs <- ms_model@Fit@filtProb
regime_dates   <- index(portfolio_ret_merged)[2:(nrow(smoothed_probs) + 1)]
n_dates <- min(length(regime_dates), nrow(smoothed_probs), nrow(filtered_probs))
regime_dates   <- regime_dates[1:n_dates]
smoothed_probs <- smoothed_probs[1:n_dates, , drop = FALSE]
filtered_probs <- filtered_probs[1:n_dates, , drop = FALSE]
colnames(smoothed_probs) <- c("P_Regime1", "P_Regime2")
colnames(filtered_probs) <- c("F_Regime1", "F_Regime2")

regime_df <- data.frame(
  Date = regime_dates,
  P_Regime1 = smoothed_probs[,1], P_Regime2 = smoothed_probs[,2],
  F_Regime1 = filtered_probs[,1], F_Regime2 = filtered_probs[,2],
  Regime = ifelse(smoothed_probs[,2] > 0.5, 2, 1))

r1_vol <- sd(ms_df$ret[regime_df$Regime == 1], na.rm = TRUE)
r2_vol <- sd(ms_df$ret[regime_df$Regime == 2], na.rm = TRUE)
r1_mean <- mean(ms_df$ret[regime_df$Regime == 1], na.rm = TRUE)
r2_mean <- mean(ms_df$ret[regime_df$Regime == 2], na.rm = TRUE)
high_vol <- if (r2_vol > r1_vol) 2 else 1
low_vol  <- if (r2_vol > r1_vol) 1 else 2

# === v8 CHANGE (6): Stress/Calm naming for all reporting ===
stress_vol  <- if (high_vol == 2) r2_vol  else r1_vol
calm_vol    <- if (low_vol  == 2) r2_vol  else r1_vol
stress_mean <- if (high_vol == 2) r2_mean else r1_mean
calm_mean   <- if (low_vol  == 2) r2_mean else r1_mean
stress_n    <- sum(regime_df$Regime == high_vol)
calm_n      <- sum(regime_df$Regime == low_vol)
message("\n--- Regime Characterisation (Stress / Calm) ---")
message(sprintf("  STRESS (high-vol): mean = %+.4f%% | sigma = %.4f%% | N = %d (%.1f%%)",
                stress_mean*100, stress_vol*100, stress_n, 100*stress_n/nrow(regime_df)))
message(sprintf("  CALM  (low-vol):   mean = %+.4f%% | sigma = %.4f%% | N = %d (%.1f%%)",
                calm_mean*100, calm_vol*100, calm_n, 100*calm_n/nrow(regime_df)))
message("\n--- Transition Probability Matrix ---"); print(ms_model@transMat)
diag_tm <- diag(ms_model@transMat)
message(sprintf("  Expected durations: %.0f and %.0f trading days",
                1/(1-diag_tm[1]), 1/(1-diag_tm[2])))

# 8.2  Regime probability plot (smoothed, for visualisation)
highvol_col <- paste0("P_Regime", high_vol)
regime_df$P_HighVol <- regime_df[[highvol_col]]
p_regime <- ggplot(regime_df, aes(Date, P_HighVol)) +
  geom_ribbon(aes(ymin = 0, ymax = P_HighVol), fill = "firebrick", alpha = 0.25) +
  geom_line(color = "firebrick3", linewidth = 0.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey30", linewidth = 0.4) +
  geom_vline(data = events, aes(xintercept = date), linetype = "dotted",
             color = "grey30", linewidth = 0.6, inherit.aes = FALSE) +
  geom_text(data = events, aes(x = date + 30, y = y_pos, label = label),
            size = 2.3, hjust = 0, color = "grey20", inherit.aes = FALSE) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(title = sprintf("Smoothed Probability of the Stress (High-Volatility) Regime (sigma = %.2f%%)",
                       stress_vol*100),
       subtitle = sprintf("Above the 50%% line = stress regime | %.0f%% of trading days (smoothed)",
                          100*stress_n/nrow(regime_df)),
       x = NULL, y = "P(Stress regime)")
print(p_regime)

# 8.3  Returns coloured by regime  (merged sample)
ret_regime_df <- merge(
  data.frame(Date = index(portfolio_ret_merged), Return = as.numeric(portfolio_ret_merged)*100),
  regime_df[, c("Date", "Regime")], by = "Date", all.x = TRUE)
ret_regime_df$Regime[is.na(ret_regime_df$Regime)] <- low_vol
p_ret_regime <- ggplot(ret_regime_df, aes(Date, Return, color = factor(Regime))) +
  geom_line(linewidth = 0.35, alpha = 0.85) +
  scale_color_manual(values = setNames(c("steelblue","firebrick"),
                                        as.character(c(low_vol, high_vol))),
                     labels = setNames(c("Calm (low-vol)","Stress (high-vol)"),
                                       as.character(c(low_vol, high_vol))), name = "") +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "Portfolio Returns Coloured by Market Regime",
       subtitle = "Blue = Calm (low-vol) | Red = Stress (high-vol)",
       x = NULL, y = "Log Return (%)")
print(p_ret_regime)


# =============================================================================
# SECTION 9 — REGIME-CONDITIONAL NARDL (FULL-SAMPLE DUMMY INTERACTION)
# =============================================================================
message("=== SECTION 9: REGIME-CONDITIONAL NARDL (interaction) ===")

# 9.1  D_stress from FILTERED probabilities (no look-ahead bias)
nardl_regime <- merge(nardl_data,
  regime_df[, c("Date","Regime","F_Regime1","F_Regime2")], by = "Date", all.x = TRUE)
nardl_regime$Regime[is.na(nardl_regime$Regime)] <- low_vol
filt_highvol_col <- paste0("F_Regime", high_vol)
nardl_regime$P_stress_filtered <- nardl_regime[[filt_highvol_col]]
nardl_regime$P_stress_filtered[is.na(nardl_regime$P_stress_filtered)] <- 0
nardl_regime$D_stress <- as.integer(nardl_regime$P_stress_filtered > 0.5)
message(sprintf("Filtered stress dummy: %d stress (%.1f%%) | %d calm (%.1f%%) | no look-ahead.",
                sum(nardl_regime$D_stress), 100*mean(nardl_regime$D_stress),
                sum(1-nardl_regime$D_stress), 100*mean(1-nardl_regime$D_stress)))

# Build lags / differences on the CONTINUOUS calendar series
nr <- nardl_regime
nr$lnPort_L1 <- c(NA, nr$lnPortfolio[-nrow(nr)])
nr$POS_IO_L1 <- c(NA, nr$POS_IO[-nrow(nr)])
nr$NEG_IO_L1 <- c(NA, nr$NEG_IO[-nrow(nr)])
nr$AUDUSD_L1 <- c(NA, nr$lnAUDUSD[-nrow(nr)])
nr$VIX_L1    <- c(NA, nr$lnVIX[-nrow(nr)])
nr$SP500_L1  <- c(NA, nr$lnSP500[-nrow(nr)])
nr$d_Port   <- c(NA, diff(nr$lnPortfolio))
nr$d_POS_IO <- c(NA, diff(nr$POS_IO))
nr$d_NEG_IO <- c(NA, diff(nr$NEG_IO))
for (j in 1:p_sel)
  nr[[paste0("d_Port_L", j)]] <- c(rep(NA, j), nr$d_Port[-((nrow(nr)-j+1):nrow(nr))])
# q lagged commodity differences (kept consistent with the main spec)
if (q_sel > 0) for (j in 1:q_sel) {
  nr[[paste0("d_POS_IO_L", j)]] <- c(rep(NA, j), nr$d_POS_IO[-((nrow(nr)-j+1):nrow(nr))])
  nr[[paste0("d_NEG_IO_L", j)]] <- c(rep(NA, j), nr$d_NEG_IO[-((nrow(nr)-j+1):nrow(nr))])
}
# interactions
nr$D_x_POS_IO_L1 <- nr$D_stress * nr$POS_IO_L1
nr$D_x_NEG_IO_L1 <- nr$D_stress * nr$NEG_IO_L1
nr$D_x_d_POS_IO  <- nr$D_stress * nr$d_POS_IO
nr$D_x_d_NEG_IO  <- nr$D_stress * nr$d_NEG_IO

need_cols <- c("d_Port","lnPort_L1","POS_IO_L1","NEG_IO_L1","AUDUSD_L1","VIX_L1",
               "SP500_L1","D_x_POS_IO_L1","D_x_NEG_IO_L1","D_x_d_POS_IO","D_x_d_NEG_IO",
               paste0("d_Port_L", 1:p_sel))
if (q_sel > 0) need_cols <- c(need_cols, paste0("d_POS_IO_L", 1:q_sel),
                              paste0("d_NEG_IO_L", 1:q_sel))
nr_clean <- nr[complete.cases(nr[, need_cols]), ]
message(sprintf("Interaction model sample: %d obs", nrow(nr_clean)))

# 9.2  Estimate interaction NARDL
ar_cols  <- paste(paste0("d_Port_L", 1:p_sel), collapse = " + ")
q_cols   <- if (q_sel > 0)
  paste("+", paste(c(paste0("d_POS_IO_L", 1:q_sel), paste0("d_NEG_IO_L", 1:q_sel)),
                   collapse = " + ")) else ""
interaction_formula <- as.formula(paste0(
  "d_Port ~ lnPort_L1 + POS_IO_L1 + NEG_IO_L1 + AUDUSD_L1 + VIX_L1 + SP500_L1 +",
  "D_x_POS_IO_L1 + D_x_NEG_IO_L1 + D_stress +", ar_cols,
  " + d_POS_IO + d_NEG_IO + D_x_d_POS_IO + D_x_d_NEG_IO", q_cols))
nardl_interact <- lm(interaction_formula, data = nr_clean)
hac_interact   <- coeftest(nardl_interact, vcov = NeweyWest(nardl_interact, lag = 4))
message("\n--- Interaction NARDL: HAC-Robust Coefficients ---"); print(hac_interact)

# 9.3  Regime-conditional multipliers
b_int <- coef(nardl_interact); rho_int <- b_int["lnPort_L1"]
theta_pos <- b_int["POS_IO_L1"]; theta_neg <- b_int["NEG_IO_L1"]
delta_pos <- b_int["D_x_POS_IO_L1"]; delta_neg <- b_int["D_x_NEG_IO_L1"]
pi_pos <- b_int["d_POS_IO"]; pi_neg <- b_int["d_NEG_IO"]
gamma_pos <- b_int["D_x_d_POS_IO"]; gamma_neg <- b_int["D_x_d_NEG_IO"]
Lpos_calm   <- -theta_pos/rho_int;               Lneg_calm   <- -theta_neg/rho_int
Lpos_stress <- -(theta_pos+delta_pos)/rho_int;   Lneg_stress <- -(theta_neg+delta_neg)/rho_int
message(sprintf("\n  rho = %.4f", rho_int))
message(sprintf("  CALM   : L+ = %.4f | L- = %.4f | gap = %.4f",
                Lpos_calm, Lneg_calm, Lpos_calm - Lneg_calm))
message(sprintf("  STRESS : L+ = %.4f | L- = %.4f | gap = %.4f",
                Lpos_stress, Lneg_stress, Lpos_stress - Lneg_stress))

# 9.4  Wald tests (long-run + short-run; both regimes; regime difference)
vcov_int <- NeweyWest(nardl_interact, lag = 4)
test_lr_calm   <- linearHypothesis(nardl_interact, "POS_IO_L1 = NEG_IO_L1", vcov = vcov_int)
test_lr_stress <- linearHypothesis(nardl_interact,
  "POS_IO_L1 + D_x_POS_IO_L1 = NEG_IO_L1 + D_x_NEG_IO_L1", vcov = vcov_int)
test_lr_diff   <- linearHypothesis(nardl_interact, "D_x_POS_IO_L1 = D_x_NEG_IO_L1", vcov = vcov_int)
test_sr_calm   <- linearHypothesis(nardl_interact, "d_POS_IO = d_NEG_IO", vcov = vcov_int)
test_sr_stress <- linearHypothesis(nardl_interact,
  "d_POS_IO + D_x_d_POS_IO = d_NEG_IO + D_x_d_NEG_IO", vcov = vcov_int)
test_sr_diff   <- linearHypothesis(nardl_interact, "D_x_d_POS_IO = D_x_d_NEG_IO", vcov = vcov_int)
p_lr_calm <- test_lr_calm$Pr[2]; p_lr_stress <- test_lr_stress$Pr[2]; p_lr_diff <- test_lr_diff$Pr[2]
p_sr_calm <- test_sr_calm$Pr[2]; p_sr_stress <- test_sr_stress$Pr[2]; p_sr_diff <- test_sr_diff$Pr[2]
message("\n--- Wald tests (HAC) ---")
message(sprintf("  LR asymmetry  CALM:   F=%.3f p=%.4f", test_lr_calm$F[2],   p_lr_calm))
message(sprintf("  LR asymmetry  STRESS: F=%.3f p=%.4f", test_lr_stress$F[2], p_lr_stress))
message(sprintf("  LR regime-diff:       F=%.3f p=%.4f", test_lr_diff$F[2],   p_lr_diff))
message(sprintf("  SR asymmetry  CALM:   F=%.3f p=%.4f", test_sr_calm$F[2],   p_sr_calm))
message(sprintf("  SR asymmetry  STRESS: F=%.3f p=%.4f", test_sr_stress$F[2], p_sr_stress))
message(sprintf("  SR regime-diff:       F=%.3f p=%.4f", test_sr_diff$F[2],   p_sr_diff))


# =============================================================================
# === v8 CHANGE (4): GENERATED-REGRESSOR BOOTSTRAP
# -----------------------------------------------------------------------------
# D_stress is an ESTIMATE from the Markov stage; the HAC SEs above treat it as
# fixed and known, understating uncertainty in the interaction terms. Two
# bootstraps quantify the missing uncertainty:
#
#   (A) CLASSIFICATION bootstrap (fast, always runs):
#       D_stress(b)[t] ~ Bernoulli(filtered stress prob[t]); refit the linear
#       interaction model; collect the regime-conditional quantities.
#       SCOPE: propagates uncertainty in WHICH days are stress, holding the MS
#       parameter estimates fixed. Captures the dominant practical channel.
#
#   (B) FULL TWO-STAGE block bootstrap (rigorous; slow; flag below):
#       moving-block resample the increments, re-integrate the I(1) levels,
#       REFIT the Markov model, rebuild D_stress, refit the interaction NARDL.
#       SCOPE: also propagates MS PARAMETER uncertainty. This is the version to
#       report in the full paper. Runtime ~ B * (MS fit time); set B accordingly.
#
# Report (B) if you run it; otherwise report (A) and state its scope honestly.
# =============================================================================
set.seed(20260530)

# --- (A) classification bootstrap -------------------------------------------
message("\n=== GENERATED-REGRESSOR BOOTSTRAP (A: classification uncertainty) ===")
B_class <- 999
p_stress_vec <- nr_clean$P_stress_filtered
boot_one_A <- function() {
  d <- nr_clean
  d$D_stress      <- rbinom(nrow(d), 1, pmin(pmax(p_stress_vec, 0), 1))
  d$D_x_POS_IO_L1 <- d$D_stress * d$POS_IO_L1
  d$D_x_NEG_IO_L1 <- d$D_stress * d$NEG_IO_L1
  d$D_x_d_POS_IO  <- d$D_stress * d$d_POS_IO
  d$D_x_d_NEG_IO  <- d$D_stress * d$d_NEG_IO
  m <- tryCatch(lm(interaction_formula, data = d), error = function(e) NULL)
  if (is.null(m)) return(rep(NA, 7))
  b <- coef(m); r <- b["lnPort_L1"]
  c(delta_pos   = unname(b["D_x_POS_IO_L1"]),
    delta_neg   = unname(b["D_x_NEG_IO_L1"]),
    regime_diff = unname(b["D_x_POS_IO_L1"] - b["D_x_NEG_IO_L1"]),
    Lpos_calm   = unname(-b["POS_IO_L1"]/r),
    Lneg_calm   = unname(-b["NEG_IO_L1"]/r),
    Lpos_stress = unname(-(b["POS_IO_L1"]+b["D_x_POS_IO_L1"])/r),
    Lneg_stress = unname(-(b["NEG_IO_L1"]+b["D_x_NEG_IO_L1"])/r))
}
boot_A <- replicate(B_class, boot_one_A())
boot_A <- boot_A[, colSums(is.na(boot_A)) == 0, drop = FALSE]
summ_A <- data.frame(
  quantity = rownames(boot_A),
  point    = c(delta_pos, delta_neg, delta_pos - delta_neg,
               Lpos_calm, Lneg_calm, Lpos_stress, Lneg_stress),
  boot_SE  = apply(boot_A, 1, sd),
  ci_lo    = apply(boot_A, 1, quantile, 0.025),
  ci_hi    = apply(boot_A, 1, quantile, 0.975), row.names = NULL)
message(sprintf("  (A) %d successful replications", ncol(boot_A)))
print(summ_A, digits = 4)

# --- (B) full two-stage moving-block bootstrap (optional / rigorous) --------
RUN_FULL_MS_BOOTSTRAP <- TRUE     # set FALSE to skip the slow rigorous version
B_full     <- 200                 # increase for the final paper (e.g. 999); slow
block_len  <- 21                  # ~1 trading month; preserves serial dependence

if (RUN_FULL_MS_BOOTSTRAP) {
  message("\n=== GENERATED-REGRESSOR BOOTSTRAP (B: full two-stage MBB) ===")
  message(sprintf("  WARNING: refits the Markov model %d times; this is SLOW.", B_full))

  base <- nardl_data[order(nardl_data$Date), ]
  Tn   <- nrow(base)
  # stationary building blocks aligned by row:
  inc <- data.frame(
    dPort = c(NA, diff(base$lnPortfolio)),
    dIO   = c(NA, diff(base$lnIronOre)),
    dFX   = c(NA, diff(base$lnAUDUSD)),
    lnVIX = base$lnVIX,    # I(0) levels carried directly
    lnSP  = base$lnSP500)[-1, ]
  l0_port <- base$lnPortfolio[1]; l0_io <- base$lnIronOre[1]; l0_fx <- base$lnAUDUSD[1]
  Ti      <- nrow(inc)
  n_blocks <- ceiling(Ti / block_len)

  mbb_indices <- function() {
    starts <- sample.int(Ti - block_len + 1, n_blocks, replace = TRUE)
    idx <- as.vector(sapply(starts, function(s) s:(s + block_len - 1)))
    idx[1:Ti]
  }

  one_full_boot <- function() {
    idx <- mbb_indices()
    bi  <- inc[idx, ]
    # re-integrate I(1) levels from resampled increments
    lnPort_b <- cumsum(c(l0_port, bi$dPort))
    lnIO_b   <- cumsum(c(l0_io,   bi$dIO))
    lnFX_b   <- cumsum(c(l0_fx,   bi$dFX))
    lnVIX_b  <- c(bi$lnVIX[1], bi$lnVIX)   # length match
    lnSP_b   <- c(bi$lnSP[1],  bi$lnSP)
    Tb <- length(lnPort_b)
    # portfolio returns for the MS stage
    pr <- diff(lnPort_b)
    msdf <- data.frame(ret = pr[-1], lag1 = pr[-length(pr)])
    lmb  <- lm(ret ~ lag1, data = msdf)
    msb  <- tryCatch(MSwM::msmFit(lmb, k = 2, sw = c(TRUE, TRUE, TRUE), p = 0,
                       control = list(parallel = FALSE, maxiter = 500, tol = 1e-6)),
                     error = function(e) NULL)
    if (is.null(msb)) return(rep(NA, 3))
    fp <- msb@Fit@filtProb
    reg <- ifelse(fp[,2] > 0.5, 2, 1)
    v1 <- sd(msdf$ret[reg == 1], na.rm = TRUE); v2 <- sd(msdf$ret[reg == 2], na.rm = TRUE)
    hv <- if (v2 > v1) 2 else 1
    D  <- as.integer(fp[, hv] > 0.5)
    # partial sums from resampled iron-ore increments
    dIO_b <- diff(lnIO_b)
    POSb <- cumsum(pmax(dIO_b, 0)); NEGb <- cumsum(pmin(dIO_b, 0))
    # assemble model frame (align lengths: differenced series lose 1 obs)
    m <- min(length(diff(lnPort_b)), length(POSb), nrow(fp))
    dfb <- data.frame(
      d_Port    = diff(lnPort_b)[1:m],
      lnPort_L1 = lnPort_b[1:m],
      POS_IO_L1 = POSb[1:m], NEG_IO_L1 = NEGb[1:m],
      AUDUSD_L1 = lnFX_b[1:m], VIX_L1 = lnVIX_b[1:m], SP500_L1 = lnSP_b[1:m],
      d_POS_IO = c(NA, diff(POSb))[1:m], d_NEG_IO = c(NA, diff(NEGb))[1:m],
      D_stress = D[1:m])
    for (j in 1:p_sel) dfb[[paste0("d_Port_L", j)]] <- c(rep(NA, j), dfb$d_Port[1:(m-j)])
    dfb$D_x_POS_IO_L1 <- dfb$D_stress * dfb$POS_IO_L1
    dfb$D_x_NEG_IO_L1 <- dfb$D_stress * dfb$NEG_IO_L1
    dfb$D_x_d_POS_IO  <- dfb$D_stress * dfb$d_POS_IO
    dfb$D_x_d_NEG_IO  <- dfb$D_stress * dfb$d_NEG_IO
    f_b <- as.formula(paste0(
      "d_Port ~ lnPort_L1 + POS_IO_L1 + NEG_IO_L1 + AUDUSD_L1 + VIX_L1 + SP500_L1 +",
      "D_x_POS_IO_L1 + D_x_NEG_IO_L1 + D_stress +",
      paste(paste0("d_Port_L", 1:p_sel), collapse = " + "),
      "+ d_POS_IO + d_NEG_IO + D_x_d_POS_IO + D_x_d_NEG_IO"))
    mb <- tryCatch(lm(f_b, data = dfb[complete.cases(dfb), ]), error = function(e) NULL)
    if (is.null(mb)) return(rep(NA, 3))
    bb <- coef(mb)
    c(delta_pos   = unname(bb["D_x_POS_IO_L1"]),
      delta_neg   = unname(bb["D_x_NEG_IO_L1"]),
      regime_diff = unname(bb["D_x_POS_IO_L1"] - bb["D_x_NEG_IO_L1"]))
  }

  boot_B <- matrix(NA, nrow = 3, ncol = B_full,
                   dimnames = list(c("delta_pos","delta_neg","regime_diff"), NULL))
  ok <- 0
  for (b in 1:B_full) {
    res <- tryCatch(one_full_boot(), error = function(e) rep(NA, 3))
    boot_B[, b] <- res
    if (all(is.finite(res))) ok <- ok + 1
    if (b %% 25 == 0) message(sprintf("    ... %d / %d reps (%d usable)", b, B_full, ok))
  }
  boot_B <- boot_B[, colSums(is.na(boot_B)) == 0, drop = FALSE]
  summ_B <- data.frame(
    quantity = rownames(boot_B),
    point    = c(delta_pos, delta_neg, delta_pos - delta_neg),
    boot_SE  = apply(boot_B, 1, sd),
    ci_lo    = apply(boot_B, 1, quantile, 0.025),
    ci_hi    = apply(boot_B, 1, quantile, 0.975), row.names = NULL)
  message(sprintf("  (B) %d usable replications (block length %d)", ncol(boot_B), block_len))
  print(summ_B, digits = 4)
} else {
  message("\n(Full two-stage MS bootstrap skipped: RUN_FULL_MS_BOOTSTRAP = FALSE)")
}


# =============================================================================
# SECTION 10 — INDIVIDUAL FIRM HETEROGENEITY   [now BIC-consistent]
# =============================================================================
# === v8 CHANGE (1): firms use the SAME (p_sel, q_sel) as the portfolio model ==
message("\n=== SECTION 10: FIRM-LEVEL HETEROGENEITY (BIC-consistent) ===")
firm_ln <- na.omit(merge(ln_bhp, ln_rio, ln_fmg, join = "inner"))
colnames(firm_ln) <- c("lnBHP", "lnRIO", "lnFMG")
firm_data <- merge(nardl_data,
  data.frame(Date = as.Date(index(firm_ln)),
             lnBHP = as.numeric(firm_ln[,"lnBHP"]),
             lnRIO = as.numeric(firm_ln[,"lnRIO"]),
             lnFMG = as.numeric(firm_ln[,"lnFMG"])), by = "Date", all.x = TRUE)
firm_data <- na.omit(firm_data)

firms <- list(list(name = "BHP Group", var = "lnBHP"),
              list(name = "Rio Tinto", var = "lnRIO"),
              list(name = "Fortescue", var = "lnFMG"))

firm_results <- lapply(firms, function(f) {
  message(sprintf("\n--- NARDL: %s ---", f$name))
  ts_firm <- ts(firm_data[, c(f$var, "POS_IO","NEG_IO","lnAUDUSD","lnVIX","lnSP500")],
                frequency = 252)
  colnames(ts_firm)[1] <- "lnEquity"
  f_form <- build_nardl_formula("lnEquity","POS_IO","NEG_IO", ctrl_main, p_sel, q_sel)
  mod <- tryCatch(dynlm(f_form, data = ts_firm),
                  error = function(e) { message("Error: ", e$message); NULL })
  if (is.null(mod)) return(NULL)
  V <- NeweyWest(mod, lag = 4); b <- coef(mod); r <- b["L(lnEquity, 1)"]
  asym <- tryCatch(linearHypothesis(mod, "L(POS_IO, 1) = L(NEG_IO, 1)", vcov = V),
                   error = function(e) NULL)
  Lpos <- -b["L(POS_IO, 1)"]/r; Lneg <- -b["L(NEG_IO, 1)"]/r
  p_val <- if (!is.null(asym)) asym$Pr[2] else NA
  message(sprintf("  L+ = %.4f | L- = %.4f | Wald p = %.4f", Lpos, Lneg, p_val))
  list(firm = f$name, Lpos = Lpos, Lneg = Lneg, p_value = p_val)
})

message("\n--- Firm-Level Comparison ---")
message(sprintf("%-12s %-9s %-9s %-10s %-9s", "Firm", "L+", "L-", "gap", "p"))
message(strrep("-", 52))
for (r in firm_results) if (!is.null(r))
  message(sprintf("%-12s %-9.4f %-9.4f %-10.4f %-9.4f",
                  r$firm, r$Lpos, r$Lneg, r$Lpos - r$Lneg, r$p_value))
message("\nNOTE: report the BHP-vs-Rio ordering honestly — the pure-play (FMG) shows")
message("3-4x the multipliers, but BHP and Rio are close and need not be monotonic")
message("in iron-ore share. Avoid over-claiming a strictly monotonic gradient.\n")


# =============================================================================
# SECTION 11 — PUBLICATION PLOTS
# =============================================================================
message("=== GENERATING FINAL PLOTS ===")

# 11.1  Main long-run asymmetry
asym_main_df <- data.frame(
  Direction = factor(c("L+ (positive shocks)","L- (negative shocks)"),
                     levels = c("L+ (positive shocks)","L- (negative shocks)")),
  Multiplier = c(L_pos, L_neg))
p_asym_main <- ggplot(asym_main_df, aes(Direction, Multiplier, fill = Direction)) +
  geom_col(width = 0.5, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.3f", Multiplier)), vjust = -0.4, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("L+ (positive shocks)" = "steelblue",
                               "L- (negative shocks)" = "firebrick3")) +
  scale_y_continuous(limits = c(0, max(L_pos, L_neg)*1.2)) +
  annotate("segment", x = 1, xend = 2, y = max(L_pos, L_neg)*1.08,
           yend = max(L_pos, L_neg)*1.08, linewidth = 0.6, color = "grey30") +
  annotate("text", x = 1.5, y = max(L_pos, L_neg)*1.13,
           label = sprintf("Gap = %.3f\n(Wald p = %.3g)", abs(L_pos-L_neg), asym_test_iron$Pr[2]),
           size = 3.2, color = "grey20") +
  labs(title = "Long-Run Multipliers: Iron Ore Shock Transmission (Portfolio)",
       subtitle = sprintf("NARDL (n = %d) | S&P 500 control | HAC SEs", nrow(nardl_data)),
       x = NULL, y = "Long-run multiplier")
print(p_asym_main)

# 11.2  Combined regime chart (long-run + short-run)
combined_regime_df <- rbind(
  data.frame(Regime = rep(c("Calm","Stress"), each = 2),
             Direction = rep(c("Positive shock","Negative shock"), 2),
             Value = c(Lpos_calm, Lneg_calm, Lpos_stress, Lneg_stress), Run = "Long-run"),
  data.frame(Regime = rep(c("Calm","Stress"), each = 2),
             Direction = rep(c("Positive shock","Negative shock"), 2),
             Value = c(pi_pos, pi_neg, pi_pos+gamma_pos, pi_neg+gamma_neg), Run = "Short-run"))
combined_regime_df$Direction <- factor(combined_regime_df$Direction,
                                       levels = c("Positive shock","Negative shock"))
combined_regime_df$Regime <- factor(combined_regime_df$Regime, levels = c("Calm","Stress"))
combined_regime_df$Run    <- factor(combined_regime_df$Run, levels = c("Long-run","Short-run"))
p_regime_combined <- ggplot(combined_regime_df, aes(Regime, Value, fill = Direction)) +
  geom_col(position = position_dodge(width = 0.55), width = 0.45) +
  geom_text(aes(label = sprintf("%.3f", Value)), position = position_dodge(width = 0.55),
            vjust = -0.4, size = 2.8, fontface = "bold") +
  scale_fill_manual(values = c("Positive shock" = "steelblue", "Negative shock" = "firebrick3"),
                    name = "Iron ore shock direction") +
  facet_wrap(~ Run, scales = "free_y") +
  labs(title = "Regime-Conditional Asymmetric Transmission: Long-Run and Short-Run",
       subtitle = "D_stress from filtered Markov probabilities — no look-ahead bias",
       x = NULL, y = "Multiplier") +
  theme(legend.position = "top", strip.text = element_text(face = "bold", size = 10))
print(p_regime_combined)

# 11.3  Market-control robustness
if (!is.null(m1) && !is.null(m2) && !is.null(m3)) {
  robust_df <- data.frame(
    Model = rep(c("No market\ncontrol","S&P 500\n(main)","ASX 200\n(comparison)"), each = 2),
    Direction = rep(c("L+ (positive)","L- (negative)"), 3),
    Value = c(m1$Lpos, m1$Lneg, m2$Lpos, m2$Lneg, m3$Lpos, m3$Lneg))
  robust_df$Model <- factor(robust_df$Model,
    levels = c("No market\ncontrol","S&P 500\n(main)","ASX 200\n(comparison)"))
  robust_df$Direction <- factor(robust_df$Direction, levels = c("L+ (positive)","L- (negative)"))
  p_asym_robust <- ggplot(robust_df, aes(Model, Value, fill = Direction)) +
    geom_col(position = position_dodge(width = 0.55), width = 0.45) +
    geom_text(aes(label = sprintf("%.3f", Value)), position = position_dodge(width = 0.55),
              vjust = -0.4, size = 3, fontface = "bold") +
    scale_fill_manual(values = c("L+ (positive)" = "steelblue", "L- (negative)" = "firebrick3"),
                      name = "") +
    scale_y_continuous(limits = c(0, max(robust_df$Value)*1.2)) +
    labs(title = "Robustness: Long-Run Multipliers Across Market Control Specifications",
         subtitle = "L+ > L- confirmed across controls", x = NULL, y = "Long-run multiplier") +
    theme(legend.position = "top")
  print(p_asym_robust)
}

# 11.4  Firm comparison
firm_plot_df <- data.frame(
  Firm = rep(c("BHP\n(diversified)","Rio Tinto\n(diversified)",
               "Fortescue\n(pure iron ore)","Portfolio\n(equal weight)"), each = 2),
  Direction = rep(c("L+ (positive)","L- (negative)"), 4),
  Value = c(firm_results[[1]]$Lpos, firm_results[[1]]$Lneg,
            firm_results[[2]]$Lpos, firm_results[[2]]$Lneg,
            firm_results[[3]]$Lpos, firm_results[[3]]$Lneg, L_pos, L_neg))
firm_plot_df$Firm <- factor(firm_plot_df$Firm,
  levels = c("BHP\n(diversified)","Rio Tinto\n(diversified)",
             "Fortescue\n(pure iron ore)","Portfolio\n(equal weight)"))
firm_plot_df$Direction <- factor(firm_plot_df$Direction, levels = c("L+ (positive)","L- (negative)"))
p_firm_comparison <- ggplot(firm_plot_df, aes(Firm, Value, fill = Direction)) +
  geom_col(position = position_dodge(width = 0.55), width = 0.45) +
  geom_text(aes(label = sprintf("%.3f", Value)), position = position_dodge(width = 0.55),
            vjust = -0.4, size = 3, fontface = "bold") +
  scale_fill_manual(values = c("L+ (positive)" = "steelblue", "L- (negative)" = "firebrick3"),
                    name = "") +
  scale_y_continuous(limits = c(0, max(firm_plot_df$Value, na.rm = TRUE)*1.2)) +
  labs(title = "Long-Run Multipliers by Firm: Iron Ore Shock Transmission",
       subtitle = "Pure-play (Fortescue) shows the largest multipliers and widest gap",
       x = NULL, y = "Long-run multiplier") +
  theme(legend.position = "top")
print(p_firm_comparison)

p_final <- (p_returns | p_ironore) / (p_regime | p_ret_regime) +
  patchwork::plot_annotation(
    title = "Asymmetric Commodity Shock Transmission in Australian Mining Equities",
    subtitle = "BHP / RIO / FMG equally weighted portfolio | 2011 – present",
    theme = theme(plot.title = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(size = 10, color = "grey40")))
print(p_final)


# =============================================================================
# SECTION 12 — CUSUM ON THE EXACT MAIN MODEL  + SAVE OUTPUTS
# =============================================================================
# === v8 CHANGE (7): CUSUM is run on the SAME model as Table 3 (the main NARDL,
# including the q lagged-difference terms), not a reduced reconstruction.
message("\n=== CUSUM STABILITY (exact main NARDL model) ===")

# Rebuild the main model as an lm on an explicit data frame so efp() matches it.
cusum_df <- nardl_data[order(nardl_data$Date), ]
cusum_df$d_Port    <- c(NA, diff(cusum_df$lnPortfolio))
cusum_df$lnPort_L1 <- c(NA, cusum_df$lnPortfolio[-nrow(cusum_df)])
cusum_df$POS_IO_L1 <- c(NA, cusum_df$POS_IO[-nrow(cusum_df)])
cusum_df$NEG_IO_L1 <- c(NA, cusum_df$NEG_IO[-nrow(cusum_df)])
cusum_df$AUDUSD_L1 <- c(NA, cusum_df$lnAUDUSD[-nrow(cusum_df)])
cusum_df$VIX_L1    <- c(NA, cusum_df$lnVIX[-nrow(cusum_df)])
cusum_df$SP500_L1  <- c(NA, cusum_df$lnSP500[-nrow(cusum_df)])
cusum_df$d_POS_IO  <- c(NA, diff(cusum_df$POS_IO))
cusum_df$d_NEG_IO  <- c(NA, diff(cusum_df$NEG_IO))
for (j in 1:p_sel)
  cusum_df[[paste0("d_Port_L", j)]] <- c(rep(NA, j), cusum_df$d_Port[-((nrow(cusum_df)-j+1):nrow(cusum_df))])
if (q_sel > 0) for (j in 1:q_sel) {
  cusum_df[[paste0("d_POS_IO_L", j)]] <- c(rep(NA, j), cusum_df$d_POS_IO[-((nrow(cusum_df)-j+1):nrow(cusum_df))])
  cusum_df[[paste0("d_NEG_IO_L", j)]] <- c(rep(NA, j), cusum_df$d_NEG_IO[-((nrow(cusum_df)-j+1):nrow(cusum_df))])
}
ar_c <- paste(paste0("d_Port_L", 1:p_sel), collapse = " + ")
q_c  <- if (q_sel > 0) paste("+", paste(c(paste0("d_POS_IO_L", 1:q_sel),
          paste0("d_NEG_IO_L", 1:q_sel)), collapse = " + ")) else ""
cusum_formula <- as.formula(paste0(
  "d_Port ~ lnPort_L1 + POS_IO_L1 + NEG_IO_L1 + AUDUSD_L1 + VIX_L1 + SP500_L1 + ",
  ar_c, " + d_POS_IO + d_NEG_IO", q_c))
cusum_df_clean <- cusum_df[complete.cases(cusum_df[, all.vars(cusum_formula)]), ]
cusum_test <- strucchange::efp(cusum_formula, data = cusum_df_clean, type = "OLS-CUSUM")
print(sctest(cusum_test))
plot(cusum_test)

# Save plots + key tables
pdf("mining_asymmetry_results_v8.pdf", width = 12, height = 8, onefile = TRUE)
  print(p_ironore); print(p_returns); print(p_rollcor); print(p_decomp)
  print(p_asym_main); print(p_regime_combined)
  if (exists("p_asym_robust")) print(p_asym_robust)
  print(p_regime); print(p_ret_regime); print(p_firm_comparison); print(p_final)
dev.off()
write.csv(desc_stats,  "01_descriptive_statistics.csv", row.names = FALSE)
write.csv(adf_results, "02_adf_unit_root_tests.csv",   row.names = FALSE)
write.csv(regime_df,   "03_regime_classification.csv",  row.names = FALSE)
write.csv(lag_robust,  "04_lag_robustness.csv",         row.names = FALSE)
write.csv(summ_A,      "05_bootstrap_classification.csv", row.names = FALSE)
if (RUN_FULL_MS_BOOTSTRAP && exists("summ_B"))
  write.csv(summ_B,    "06_bootstrap_full_twostage.csv", row.names = FALSE)

message("\n=== ANALYSIS COMPLETE (v8) ===")
message("New/updated outputs to feed back: lag-robustness grid, bounds (F+t/pssbounds),")
message("BIC-consistent copper/control/firm multipliers, regime table (merged sample,")
message("Stress/Calm), and the generated-regressor bootstrap SEs (A, and B if run).")
