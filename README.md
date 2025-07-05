# GARCH-MIDAS Modeling of S&P 500 Realized Volatility

This repository contains an R-based pipeline for analyzing the realized volatility of the S&P 500 index using the **GARCH-MIDAS** (Mixed Data Sampling) framework. The project pulls data from the Federal Reserve Economic Data (FRED) API, processes daily returns, constructs realized volatility measures, and estimates a GARCH-MIDAS model to capture both short- and long-term components of market volatility.

---

## Objectives

- Retrieve and process **S&P 500 daily price data**
- Calculate **daily returns** and **realized volatility** (fixed-time and rolling-window)
- Test statistical properties: stationarity, normality, ARCH effects
- Fit a **GARCH-MIDAS** model using monthly realized volatility as the long-run risk driver
- Analyze residuals and diagnostic statistics

---

## Dependencies

Make sure the following R packages are installed:

```r
install.packages(c(
  "quantmod", "xts", "zoo", "fredr", "dplyr", "mfGARCH", "lubridate", "tseries",
  "statmod", "MASS", "moments", "FinTS", "goftest"
))
````

---

## Setup

1. Get a [FRED API key](https://fred.stlouisfed.org/docs/api/api_key.html)
2. Set your API key in the script:

   ```r
   fredr_set_key("YOUR_API_KEY")
   ```

---

## Methodology Overview

### 1. **Data Collection**

* Uses `fredr` to download daily S\&P 500 index levels from 2015 to 2025.
* Non-trading days are automatically removed via `na.omit()`.

### 2. **Return Calculation**

* Computes daily log returns in percentage form.
* Basic statistical summaries and tests (Jarque-Bera, ADF, ARCH-LM).

### 3. **Realized Volatility Construction**

* **Fixed-time (monthly)** realized volatility via sum of squared returns.
* **Rolling-window** realized volatility via sliding window aggregation.

### 4. **GARCH-MIDAS Estimation**

* Estimated using `mfGARCH::fit_mfgarch()` with monthly fixed RV as the low-frequency regressor.
* Extracts short-run (G), long-run (Tau), and conditional volatility components.

### 5. **Model Diagnostics**

* Standardized residuals analyzed using ACF, PACF, histograms, Q-Q plots.
* Normality and uniformity tests: Kolmogorov-Smirnov, Cramér–von Mises.
* Box-Ljung tests on residuals and squared residuals.

---

## Visualizations

The analysis includes plots of:

* S\&P 500 price levels and returns
* Histogram of returns and realized volatility
* Estimated long- and short-run GARCH-MIDAS components
* Residual diagnostics (time series, ACF, PACF, Q-Q, histogram)

---

## Output

The script prints:

* Summary statistics for returns and volatility
* Model estimation results
* Diagnostic test outputs
* Plots of fitted components and residuals

---

## References

* Engle, R. F., Ghysels, E., & Sohn, B. (2013). *Stock Market Volatility and Macroeconomic Fundamentals*. Review of Economics and Statistics.
* [FRED API Documentation](https://fred.stlouisfed.org/docs/api/fred/)

---

## ✅ License

This project is licensed under the MIT License. See the `LICENSE` file for details.
