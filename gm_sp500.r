# load packages
library(quantmod)
library(xts)
library(zoo)
library(fredr)
library(dplyr)
library(mfGARCH)
library(lubridate)
library(tseries)
library(statmod)
library(MASS)
library(moments)
library(FinTS)
library(goftest)

# set FRED API key
fredr_set_key("API_KEY")

# retriev daily S&P 500 data and exclude non-trading days
sp500_data <- fredr(
  series_id = "SP500",
  frequency = "d",
  observation_start = as.Date("2015-04-16"),
  observation_end = as.Date("2025-04-16")
) %>% na.omit()

# convert data to xts object
sp500_xts <- xts(sp500_data$value, order.by = as.Date(sp500_data$date))

# calculate daily returns (as %)
sp500_returns <- diff(log(sp500_xts)) * 100
returns_mean <- mean(sp500_returns, na.rm = TRUE)   # mean of daily returns
returns_sd <- sd(sp500_returns, na.rm = TRUE)       # standard deviation of daily returns

# plot data
plot(sp500_xts, main = "S&P 500 Daily Close")
plot(sp500_returns, main = "S&P 500 Daily Returns (in %)")

hist(sp500_returns, breaks = 50, probability = TRUE,
     main = "Histogram of S&P 500 Daily Returns",
     xlab = "Daily Returns (%)", col = "lightblue", border = "black")
curve(dnorm(x, mean = 0, sd = 1),
      col = "red", lwd = 2, add = TRUE)
legend("topleft", legend = c("Standard Normal Distribution"),
       col = c("red"), lwd = 2, bty = "n")

# run Jarque-Bera test on daily returns
jarque.bera.test(na.omit(sp500_returns))

# run ARCH-LM test on daily returns
ArchTest(na.omit(sp500_returns))

# run ADF test on daily returns
adf.test(na.omit(sp500_returns), alternative = "stationary")

# function to calculate fixed-time realized volatility
calculate_fixed_rv <- function(returns, period = "months") {
  period_endpoints <- endpoints(returns, on = period)                 # period grouping
  rv <- period.apply(returns^2, INDEX = period_endpoints, FUN = sum)  # sum of squared returns per period
  return(rv)
}

# function to calculate n-day rolling-window realized volatility
calculate_rw_rv <- function(returns, window_size) {
  rv <- rollapply(returns^2, width = window_size, FUN = sum, align = "right", fill = NA)
  return(rv)
}

# calculate and plot monthly-fixed realized volatility
monthly_fixed_rv <- calculate_fixed_rv(sp500_returns, period = "months")
plot(monthly_fixed_rv, main = "Monthly-Fixed Realized Volatility (RV)")
hist(monthly_fixed_rv, breaks = 50, probability = TRUE,
     main = "Histogram of Monthly-fixed Realized Volatility",
     xlab = "Realized Volatility", col = "lightblue", border = "black")

# merge into data frame and fill backward monthly RV
df <- merge(sp500_returns, monthly_fixed_rv, fill = NA)
df$monthly_fixed_rv <-  na.locf(df$monthly_fixed_rv, fromLast = TRUE)
df <- df[!(format(index(df), "%Y-%m") == "2015-04")]

# convert to data frame and omit NAs
df <- data.frame(
  date = index(df$sp500_returns),
  return = coredata(df$sp500_returns),
  monthly_fixed_rv = coredata(df$monthly_fixed_rv)
) %>% na.omit()

# add 'month' column to data frame
df$month <- floor_date(df$date, unit = "month")

# print summary statistic of returns
returns_summary <- data.frame(
  Min. = min(df$sp500_returns, na.rm = TRUE),
  Median = median(df$sp500_returns, na.rm = TRUE),
  Mean = mean(df$sp500_returns, na.rm = TRUE),
  Max. = max(df$sp500_returns, na.rm = TRUE),
  SD = sd(df$sp500_returns, na.rm = TRUE),
  Skewness = skewness(df$sp500_returns, na.rm = TRUE),
  Kurtosis = kurtosis(df$sp500_returns, na.rm = TRUE)
)
print("Summary Statistic of S&P 500 Returns:")
print(returns_summary)

# print summary statistic of RV
rv_summary <- data.frame(
  Min. = min(df$monthly_fixed_rv, na.rm = TRUE),
  Median = median(df$monthly_fixed_rv, na.rm = TRUE),
  Mean = mean(df$monthly_fixed_rv, na.rm = TRUE),
  Max. = max(df$monthly_fixed_rv, na.rm = TRUE),
  SD = sd(df$monthly_fixed_rv, na.rm = TRUE),
  Skewness = skewness(df$monthly_fixed_rv, na.rm = TRUE),
  Kurtosis = kurtosis(df$monthly_fixed_rv, na.rm = TRUE)
)
print("Summary Statistic of Monthly-Fixed Realized Volatility:")
print(rv_summary)

# fit GARCH-MIDAS model with monthly fixed-time span realized volatility
gm_model_fit <- fit_mfgarch(data = df,
                            y = "sp500_returns",
                            x = "monthly_fixed_rv",
                            low.freq = "month",
                            K = 12,
                            gamma = FALSE,
                            weighting = "beta.restricted")

# print parameter estimations
print(gm_model_fit)

# plot fitted values of Tau component
tau_hat <- na.omit(gm_model_fit$tau)
plot(tau_hat, lwd = 1.5, type = "l", xlab = "Time-Index", ylab = "Fitted Long-run Component")

# plot fitted values of G component and conditional volatility
g_hat <- na.omit(gm_model_fit$g)
cond_vola <- g_hat * tau_hat
plot(
  g_hat,
  lwd = 1.5,
  type = "l",
  col = "blue",
  xlab = "Time Index",
  ylab = "Fitted Short-run Component and Conditional Volatility"
)
lines(cond_vola, lwd = 1.5, col = "red")
legend("topright", legend = c("Fitted Short-run Component", "Conditional Volatility"),
       col = c("blue", "red"), lwd = 1.5)

# plot weighting scheme
plot_weighting_scheme(gm_model_fit)

# compute (squared) standardized residuals from fitted model
std_res <- na.omit(gm_model_fit$df.fitted$residuals)
sqrd_std_res <- std_res^2

# plot time series of standardized residuals from fitted model
plot(
  std_res,
  lwd = 1.5,
  type = "l",
  xlab = "Time Index",
  ylab = "Standardized Residuals",
  main = "Time Series of Standardized Residuals"
)

# plot the (partial) autocorrelation function of standardized residuals
acf(std_res, main = "ACF of Standardized Residuals", lag.max = 20)
pacf(std_res, main = "PACF of Standardized Residuals", lag.max = 20)

# create QQ plot standardized residuals from fitted model
qqnorm(std_res, main = "Q-Q Plot of Residuals")
qqline(std_res, col = "red")

# plot histogram of standardized residuals from fitted model
hist(std_res, breaks = 50, probability = TRUE,
     main = "Histogram of Standardized Residuals",
     xlab = "Standardized Residuals", col = "lightblue", border = "black")
curve(dnorm(x, mean = 0, sd = 1), col = "red", lwd = 2, add = TRUE)
legend("topleft", legend = c("Standard Normal Distribution"),
       col = c("red"), lwd = 2, bty = "n")

# run Kolmogorov-Smirnov and Cramér–von Mises test on standardized residuals
ks.test(std_res, "pnorm", mean = 0, sd = 1)               # Kolmogorov-Smirnov test
cvm.test(std_res, null = "pnorm", mean = 0, sd = 1)       # Cramér–von Mises test

# run Box-Ljung test on (squared) standardized residuals
Box.test(std_res, lag = 10, type = "Ljung-Box")           # standardized residuals
Box.test(sqrd_std_res, lag = 10, type = "Ljung-Box")      # squared standardized residuals

# calculate uniform values using Gaussian distribution
u <- pnorm(std_res, mean = mean(std_res), sd = sd(std_res))

# run Kolmogorov-Smirnov and Cramér–von Mises test for checking uniformity
ks.test(u, "punif", min = 0, max = 1)
cvm.test(u, null = "punif")