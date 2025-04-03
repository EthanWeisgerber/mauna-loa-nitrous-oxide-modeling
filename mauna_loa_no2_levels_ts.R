library(tseries)
library(forecast)

mauna_loa <- read.table("mauna_loa.txt", header = TRUE, sep = "", fill = TRUE)
head(mauna_loa)
tail(mauna_loa)
length(mauna_loa)

values <- ts(mauna_loa[,4], start=c(1997,5), frequency=12)
values


plot.ts(values, main='Monthly Nitrous Oxide Levels 1997-2023')


values_last_10 <- window(values, start=c(2014, 1))
values_last_10

plot(values_last_10, ylab='Value', main='Monthly Nitrous Oxide Levels 2014-2023', col = "black")
lines(lowess(values_last_10, f=0.10), lwd=2, lty=1, col="blue")
lines(lowess(values_last_10, f=0.05), lwd=2, lty=2, col="red")

# Decompose time series into seasonality, trend, and remainders
decomp <- stl(values_last_10, "per")
plot(decomp)
names(decomp)
head(decomp$time.series)



# Save and plot these features
seasonal <- decomp$time.series[,1]
trend <- decomp$time.series[,2]
residuals <- decomp$time.series[,3]

par(mfrow=c(1,1))
plot(seasonal, main='Seasonality')
plot(trend, main='Trend')
plot(residuals, main="Residuals")

# Use acf to see if residuals are white noise
acf(residuals, lag.max = 500, main='Mauna Loa Time Series Residuals')
# We can see that it is not quite white noise, as values are outside blue lines

# Let's try to model the residuals autoregressively to see where this dependency lies
# We consider the function ar.yw to fit the AR(p) model where the order p is 
# determined by the information criterion AIC.
fit_ar <- ar.yw(residuals, order.max = NULL)
fit_ar
# This fit suggests the order is 17, AR(17), with Wt ~ WN(0, 0.002344)

# Let's study the residuals from this AR(17) model to check for autocorrelation
names(fit_ar)
model_resid <- ts(fit_ar$resid)
model_resid
# Remove first 17 values as they are NA (No autocorrelation until after lag p=17)
acf(model_resid[-c(1:17)], main="ACF of residuals from AR(17) fit", lag.max = 500)
# This resembles white noise!
checkresiduals(model_resid[-c(1:17)])

# We have a good model for our residuals to be white noise, now let's focus on trend
# and seasonality

plot(trend)

# Extract time information from the time series object
time_index <- seq_along(values_last_10)
time_index

# Linear trend
linear_fit <- lm(as.numeric(values_last_10) ~ time_index)
linear_trend <- ts(predict(linear_fit), start = start(values_last_10), frequency = frequency(values_last_10))

plot(trend, type = "l", col = "black", lwd = 2,
     main = "Linear Model with Nitrous Oxide Levels Trend",
     ylab = "Nitrous Oxide Levels (ppm)", xlab = "Time")
lines(linear_trend, col = "blue", lwd = 2, lty = 2)
legend("topleft", legend = c("Observed Trend", "Modeled Trend"),
       col = c("black", "blue"), lty = c(1, 2), lwd = 2)

# Polynomial trend (degree 3)
poly_fit <- lm(as.numeric(values_last_10) ~ poly(time_index, 3))
poly_fit
poly_trend <- ts(predict(poly_fit), start = start(values_last_10), frequency = frequency(values_last_10))

plot(trend, type = "l", col = "black", lwd = 2,
     main = "Polynomial Model (degree 3) with Nitrous Oxide Levels Trend",
     ylab = "Nitrous Oxide Levels (ppm)", xlab = "Time")
lines(poly_trend, col = "red", lwd = 2, lty = 2)
legend("topleft", legend = c("Observed Trend", "Modeled Trend"),
       col = c("black", "red"), lty = c(1, 2), lwd = 2)

# Exponential trend
a_use <- (max(values_last_10) - min(values_last_10)) / 10
b_use <- 0.001
c_use <- mean(values_last_10)

exp_fit <- nls(as.numeric(values_last_10) ~ a * exp(b * time_index) + c, 
               start = list(a = a_use, b = b_use, c = c_use))
exp_trend <- ts(predict(exp_fit), start = start(values_last_10), frequency = frequency(values_last_10))

plot(trend, type = "l", col = "black", lwd = 2,
     main = "Exponential Model with Nitrous Oxide Levels Trend",
     ylab = "Nitrous Oxide Levels (ppm)", xlab = "Time")
lines(exp_trend, col = "purple", lwd = 2, lty = 2)
legend("topleft", legend = c("Observed Trend", "Modeled Trend"),
       col = c("black", "purple"), lty = c(1, 2), lwd = 2)


# We have an okay linear model, and pretty accurate looking exponential and polynomial models!
# They are saved as linear_trend, poly_trend, and exp_trend, respectively.

# Now we will turn our attention to modeling seasonality

plot(seasonal, main='Seasonality')

time_index <- seq_along(seasonal)

a <- 2 * pi/12
dd <- data.frame(season = seasonal, time = time_index)
head(dd)

fit_ss <- lm(season ~ I(cos(a * time)) + I(sin(a * time)) +
               I(cos(2 * a * time)) + I(sin(2 * a * time)) +
               I(cos(3 * a * time)) + I(sin(3 * a * time)), data = dd)
fit_ss
coef_ss <- summary(fit_ss)$coefficients

# Generate fitted values
fitted_seasonality <- predict(fit_ss)

plot(time_index, seasonal, type = "l", col = "black", lwd = 2,
     main = "Seasonality Fit Using Trig Functions",
     xlab = "Time Index", ylab = "Seasonality")
lines(time_index, fitted_seasonality, col = "blue", lwd = 2, lty = 2)

summary(fit_ss)

# Now we have a trigonometric trend for the seasonality, which means we have seasonality, trend,
# and white noise residuals. Let's create some models!

# Additive models are denoted by Yt = Trend(t) + Seasonality(t) + Residuals(t)

### MODEL 1 ###

# This will be an additive model with linear trend, seasonality, and white noise

# Ensure they are numeric for adding
class(linear_fit)  # Returns 'lm' - need to convert
class(fitted_seasonality)  # Returns 'numeric' - no need to convert
linear_fit <- as.numeric(linear_trend)
class(linear_fit)

set.seed(12345)
white_noise <- rnorm(length(time_index), mean=0, sd=sd(residuals)) # Create basic white noise for testing purposes

# Make our model
yhat1 <- linear_fit + fitted_seasonality + white_noise

# Plot
plot(time_index, values_last_10, type = "l", col = "black", lwd = 2,
     main = "Linear Trend + Seasonality + White Noise",
     xlab = "Time Index", ylab = "Nitrous Oxide Levels (ppm)")
lines(time_index, yhat1, col = "blue", lwd = 2, lty = 2)
legend("topleft", legend = c("Observed Data", "Modeled Data (yhat1)"),
       col = c("black", "blue"), lty = c(1, 2), lwd = 2)

# Plot Residuals
residuals_yhat1 <- values_last_10 - yhat1
plot(time_index, residuals_yhat1, type = "p", col = "blue", pch = 20,
     main = "Residuals of yhat1 Model",
     xlab = "Time Index", ylab = "Residuals")
abline(h = 0, col = "red", lty = 2)


### MODEL 2 ###

# This will be an additive model with exponential trend, seasonality, and white noise

# Extract fitted values from exponential model
exp_trend <- predict(exp_fit)

exp_trend <- as.numeric(exp_trend)

yhat2 <- exp_trend + fitted_seasonality + white_noise

plot(time_index, values_last_10, type = 'l', col = 'black', lwd = 2,
     main = "Exponential Trend + Seasonality + White Noise",
     xlab = 'Time Index', ylab = 'Nitrous Oxide Levels (ppm)')
lines(time_index, yhat2, col='red', lwd=2, lty=2)
legend("topleft", legend = c("Observed Data", "Modeled Data (yhat2)"),
       col = c("black", "red"), lty = c(1, 2), lwd = 2)

# Plot Residuals
residuals_yhat2 <- values_last_10 - yhat2
plot(time_index, residuals_yhat2, type = "p", col = "red", pch = 20,
     main = "Residuals of yhat1 Model",
     xlab = "Time Index", ylab = "Residuals")
abline(h = 0, col = "blue", lty = 2)


### MODEL 3 ###

# This will be an additive model with polynomial trend, seasonality, and white noise

class(poly_trend)
poly_trend <- as.numeric(poly_trend)

yhat3 <- poly_trend + fitted_seasonality + white_noise

plot(time_index, values_last_10, type = 'l', col = 'black', lwd = 2,
     main = "Polynomial (Degree 3) Trend + Seasonality + White Noise",
     xlab = 'Time Index', ylab = 'Nitrous Oxide Levels (ppm)')
lines(time_index, yhat3, col='green', lwd=2, lty=2)
legend("topleft", legend = c("Observed Data", "Modeled Data (yhat3)"),
       col = c("black", "green"), lty = c(1, 2), lwd = 2)

# Plot Residuals
residuals_yhat3 <- values_last_10 - yhat3
plot(time_index, residuals_yhat3, type = "p", col = "green", pch = 20,
     main = "Residuals of yhat1 Model",
     xlab = "Time Index", ylab = "Residuals")
abline(h = 0, col = "black", lty = 2)

# It looks like our second and third models are much better than our first - we can throw 
# that model out. Let's repeat yhat2 and yhat3 with the residual AR(17) model.


# Let's first fit our AR(17)

fit.w1 <- arima(residuals, order=c(17,0,0), include.mean=FALSE)
fit.w1
fit.w1$aic # -407.1485 is very low - this is good!

# AR(17) seems to be good, but let's try AR(15) and AR(20) to compare:
test_fit20 <- arima(residuals, order=c(20,0,0), include.mean=FALSE)
test_fit20$aic # -405.9633 - AR(20) is higher than AR(17), worse
test_fit30 <- arima(residuals, order=c(15,0,0), include.mean=FALSE)
test_fit30$aic # -393.767 - AR(15) is higher than AR(17), worse

coef_w1 <- coef(fit.w1)
coef_w1

# We plot the ACF for the initial residuals and the residuals from the AR(17)

par(mfrow=c(1,1))
acf(residuals, main="ACF Plot of MLO Data")
res_w1 <- ts(fit.w1$residuals)
acf(res_w1, main="ACF Plot of AR(17) Model") # This one clearly better

# Now let's recreate our two models using AR(17)

### MODEL 2 with AR(17) ###

yhat2.2 <- exp_trend + fitted_seasonality + res_w1

plot(time_index, values_last_10, type = 'l', col = 'black', lwd = 2,
     main = "Exponential Trend + Seasonality + AR(17)",
     xlab = 'Time Index', ylab = 'Nitrous Oxide Levels (ppm)')
lines(time_index, yhat2.2, col='red', lwd=2, lty=2)
legend("topleft", legend = c("Observed Data", "Modeled Data (yhat2.2)"),
       col = c("black", "red"), lty = c(1, 2), lwd = 2)


### MODEL 3 with AR(17) ###

yhat3.2 <- poly_trend + fitted_seasonality + res_w1

plot(time_index, values_last_10, type = 'l', col = 'black', lwd = 2,
     main = "Polynomial (Degree 3) Trend + Seasonality + AR(17)",
     xlab = 'Time Index', ylab = 'Nitrous Oxide Levels (ppm)')
lines(time_index, yhat2.2, col='green', lwd=2, lty=2)
legend("topleft", legend = c("Observed Data", "Modeled Data (yhat2.3)"),
       col = c("black", "green"), lty = c(1, 2), lwd = 2)


#MSE

y <- as.numeric(values_last_10)

sum((y-yhat1)^2) / length(y)     # 0.09144829
sum((y-yhat2)^2) / length(y)     # 0.04082456
sum((y-yhat3)^2) / length(y)     # 0.0399147 
sum((y-yhat2.2)^2) / length(y)   # 0.02597242
sum((y-yhat3.2)^2) / length(y)   # 0.02507394

# Compute MAPE for each model
mean(abs((y - yhat1) / y)) * 100     # 0.07777365
mean(abs((y - yhat2) / y)) * 100     # 0.05045046
mean(abs((y - yhat3) / y)) * 100     # 0.04934287
mean(abs((y - yhat2.2) / y)) * 100   # 0.04028163
mean(abs((y - yhat3.2) / y)) * 100   # 0.03879621

# yhat3.2 is the best model, save it as final_model
final_model <- yhat3.2

plot(time_index, values_last_10, type = 'l', col = 'grey', lwd = 2,
     main = "Final Model",
     xlab = 'Time Index', ylab = 'Nitrous Oxide Levels (ppm)')
lines(time_index, final_model, col='blue', lwd=2, lty=2)
legend("topleft", legend = c("Observed Data", "Modeled Data"),
       col = c("grey", "blue"), lty = c(1, 2), lwd = 2)

# Check residuals
acf(res_w1, main='Final Model Residuals', lag.max = 50)
checkresiduals(res_w1)
#No autocorrelation!





#### FORECASTING ####


# Train/Test Split to test model with actual data (2022-23)


# Split data: Training set (2014-2021) and test set (2022-2023)
training_data <- window(values_last_10, end = c(2021, 12))
test_data <- window(values_last_10, start = c(2022, 1))

# Plot the split
x_range <- range(time(training_data), time(test_data))
y_range <- range(training_data, test_data)


final_model <- ts(final_model, start = c(2014,1), frequency = 12)
final_model


final_model_last_years <- window(final_model, start = c(2022,1))

# Plot the forecast
plot(training_data, type = "l", col = "black", lwd = 2, 
     main = "Forecast for 2022-2023", xlab = "Time Index", ylab = "Nitrous Oxide Levels", xlim = x_range, ylim = y_range)
lines(test_data, col = "red", lwd = 2)
lines(final_model_last_years, col = "blue", lwd = 2, lty = 2)
legend("topleft", legend = c("Training Data", "Actual Test Data", "Forecasted Test Data"), 
       col = c("black", "red", "blue"), lty = c(1, 1, 2), lwd = 2)




# Calculate errors
forecast_errors <- test_data - final_model_last_years
mse <- mean(forecast_errors^2)  # Mean Squared Error
mape <- mean(abs(forecast_errors / test_data)) * 100  # Mean Absolute Percentage Error

# Print evaluation metrics
cat("MSE:", mse, "\n")
cat("MAPE:", mape, "%\n")

# MSE=0.0246: This is a very low value, indicating that the average squared 
# differences between the forecasted and actual values are minimal. A smaller 
# MSE generally reflects better model performance.



# MAPE=0.04%: This is an excellent result, indicating that, on average, 
# the forecasted values deviate from the actual values by only 0.04% relative 
# to their magnitude. A MAPE below 5% is typically considered very accurate for 
# most forecasting applications.


### FORECASTING (2024-25) with predict ###


plot(final_model)

# Forecast future values, let's say for 12 future periods
future_values <- forecast(final_model, h=24)

# Plot the forecast
plot(future_values, main='Future Nitrous Oxide Levels')

# Print the forecasted values
print(future_values$mean)













# Forecasting manually with model

# Generate future indices for 2024-2025
future_index_24 <- seq(from = length(values_last_10) + 1, to = length(values_last_10) + 24)

# Extend the polynomial trend for 2024-2025
future_trend_24 <- predict(poly_fit, newdata = data.frame(time_index = future_index_24))

# Extend the seasonality for 2024-2025
future_seasonality_24 <- predict(fit_ss, newdata = data.frame(time = future_index_24))

# Generate residuals for 2024-2025 using the AR(17) model
future_residuals_24 <- arima.sim(model = list(ar = coef(fit.w1)), n = 24, sd = sqrt(fit.w1$sigma2))

# Combine the components to generate the forecast
forecast_24 <- future_trend_24 + future_seasonality_24 + future_residuals_24

# Generate x-axis values matching the forecasted period (2024-2025)
future_time <- seq(from = 2024 + (1 - 1) / 12, to = 2025 + 11 / 12, by = 1 / 12)

# Check lengths for debugging
length(future_time)  # Should be 24
length(forecast_24)  # Should also be 24

# Determine x-axis range (observed + forecasted data)
x_range <- range(time(values_last_10), future_time)

# Determine y-axis range (observed + forecasted data)
y_range <- range(values_last_10, forecast_24)

# Plot the observed data
plot(values_last_10, type = "l", col = "black", lwd = 2, 
     main = "Nitrous Oxide Levels Forecast (2024–2025)", 
     xlab = "Year", ylab = "Nitrous Oxide Levels", 
     xlim = x_range, ylim = y_range)

# Add the forecasted data
lines(future_time, forecast_24, col = "blue", lwd = 2, lty = 2)

# Add a legend
legend("topleft", legend = c("Observed Data", "Forecasted Data"), 
       col = c("black", "blue"), lty = c(1, 2), lwd = 2)

forecast_24


# Holt-Winters

values_last_10

# Apply Holt-Winters exponential smoothing
hw_model <- HoltWinters(values_last_10)

# Display the model summary
print(hw_model)

# Plot the fitted model
plot(hw_model, main = "Holt-Winters Fitted Model for 2014-2024")

# Forecast the next 24 months (2 years)
hw_forecast <- forecast(hw_model, h = 24)

# Plot the forecast
plot(hw_forecast, main = "Holt-Winters Forecast for 2024-2026")

# Print forecasted values
print(hw_forecast$mean)

# Evaluate residuals of the Holt-Winters model
hw_residuals <- residuals(hw_model)
acf(hw_residuals, main = "ACF of Holt-Winters Residuals")

# Perform Ljung-Box test on residuals to check for autocorrelation
Box.test(hw_residuals, lag = 20, type = "Ljung-Box")



# Predictions from Holt-Winters for the training period
hw_predictions <- fitted(hw_model)[, "xhat"]

# Predictions from your custom model
custom_predictions <- final_model

# Compute MSE for both models
mse_custom <- mean((values_last_10 - custom_predictions)^2)
mse_hw <- mean((values_last_10 - hw_predictions)^2)

cat("Custom Model MSE:", mse_custom, "\n")
cat("Holt-Winters Model MSE:", mse_hw, "\n")

# Compute MAPE for both models
mape_custom <- mean(abs((values_last_10 - custom_predictions) / values_last_10)) * 100
mape_hw <- mean(abs((values_last_10 - hw_predictions) / values_last_10)) * 100

cat("Custom Model MAPE:", mape_custom, "%\n")
cat("Holt-Winters Model MAPE:", mape_hw, "%\n")




### ADDITIONAL TESTS ###

# Test to see which month results in highest N2O values

# Identify the month with the highest average N2O levels
library(dplyr)

# Convert the time series data to a data frame for easier manipulation
values_df <- data.frame(
  Year = floor(time(values_last_10)),
  Month = cycle(values_last_10),
  Value = as.numeric(values_last_10)
)

# Calculate the average value per month
monthly_avg <- values_df %>%
  group_by(Month) %>%
  summarise(Average = mean(Value, na.rm = TRUE)) %>%
  arrange(desc(Average))

print(monthly_avg)

# Highlight the month with the highest N2O levels
highest_month <- monthly_avg[1, ]
cat("Month with highest N2O levels:", highest_month$Month, "\n")


boxplot(Value ~ Month, data = values_df, 
        main = "Monthly Distribution of N2O Levels",
        xlab = "Month", ylab = "Nitrous Oxide Levels (ppm)",
        col = "lightblue")


# Create monthly time index
months <- 1:12




par(mfrow = c(1, 1))
plot(seasonal)
seasonal


# Extract raw N₂O values for each month
monthly_avg <- tapply(values, cycle(values), mean)

# Create a bar chart of raw N₂O averages
barplot(monthly_avg, names.arg = month.abb, col = ifelse(monthly_avg == max(monthly_avg), "blue", "gray"),
        main = "Average Monthly N2O Levels (Raw Data)",
        ylab = "N2O Levels (ppm)", xlab = "Month",
        ylim = c(min(monthly_avg) - 0.5, max(monthly_avg) + 0.5))
abline(h = mean(monthly_avg), col = "red", lty = 2)  # Add average line

# Highlight the peak month
legend("topright", legend = c("Peak Month", "Other Months"),
       fill = c("blue", "gray"), bty = "n")






