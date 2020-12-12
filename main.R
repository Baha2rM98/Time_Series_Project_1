# Created by: Bahador Mirzazadeh
# Created on: 11/30/2020

# In this section we use student number as seed due to generating three parted sets of data.
Code <- 956321025

FUN = function(x) { if (!require(x, character.only = TRUE)) {
  install.packages(x, dependencies = TRUE)
  library(x, character.only = TRUE) } }

#-----------------------------
packages = c("smooth", "purrr")
package.check <- lapply(packages, FUN)
#-----------------------------
{ set.seed(Code)
  Data1 <- sim.ces("p", b = 0.4, frequency = rdunif(1, 12, 6), obs = 90)
  plot(Data1)
  Data2 <- sim.es("MAdM", frequency = rdunif(1, 12, 6), obs = 90, phi = 0.95, persistence = c(0.1, 0.05, 0.01), randomizer = "rlnorm", meanlog = 0, sdlog = 0.015)
  plot(Data2)
  Data3 <- ts(rnorm(120, 0, 5) + rep(runif(12, -50, 50), 10) * rep(c(1:10), each = 12), frequency = rdunif(1, 12, 6))
  plot(Data3)
}
Data1$data # Data set 1
Data2$data # Data set 2
Data3 # Data set 3
#=========================================================================================

# We assign each of data sets to a variable.
dataset1 <- ts(Data1$data, frequency = 9, start = c(1, 1), end = c(10, 9))
dataset1

dataset2 <- ts(Data2$data, frequency = 12, start = c(1, 1), end = c(8, 6))
dataset2

dataset3 <- Data3
dataset3
#============================================

# Classical decomposition approach implementation.
# Additive mode
DEC.ADD <- function(DATA, s, l) {
  n <- length(DATA); W <- rep(NA, times = n); D <- rep(NA, times = n)
  s0 <- rep(0, times = n); m <- rep(0, times = s); y <- 0
  forecast <- 0; shat <- rep(0, times = n); seaseffect <- 0
  Time <- seq(1:n)
  q <- floor(s / 2)
  if (s %% 2 == 0)
    for (i in (q + 1):(n - q)) {
      W[i] <- (0.5 * (DATA[i - q] + DATA[i + q]) + sum(DATA[(i - q + 1):(i + q - 1)])) / s }

  else

    for (i in (q + 1):(n - q)) {
      W[i] <- mean(DATA[(i - q):(i + q)]) }

  for (i in (q + 1):(n - q)) {
    D[i] <- DATA[i] - W[i]
    if (i %% s == 0) {
      s0[s] <- s0[s] + D[i]
      m[s] <- m[s] + 1 }
    else
      for (j in 1:(s - 1)) {
        if (i %% s == j) {
          s0[j] <- s0[j] + D[i]
          m[j] <- m[j] + 1 }
      }
  }
  for (j in 1:s) { seaseffect[j] <- s0[j] / m[j] - mean(s0[1:s] / m) }

  shat <- rep(seaseffect, length.out = n)
  y <- DATA - shat
  betahat <- lm(y ~ Time)$coefficients
  yhat <- predict(lm(y ~ Time))
  forecast <- as.vector(yhat + shat)
  error <- as.vector(DATA - forecast)
  FORECAST <- 0
  for (i in 1:l) {
    if ((n + i) %% s == 0) {
      FORECAST[i] <- (betahat[1] + betahat[2] * (n + i)) + seaseffect[s]
    }
    else FORECAST[i] <- (betahat[1] + betahat[2] * (n + i)) + seaseffect[(n + i) %% s]
  }
  DECAD <- data.frame(DATA, W, D, shat, y, yhat, forecast, error, seaseffect, betahat)
  return(DECAD)
}

# =============
# Product mode
DEC.MULT <- function(DATA, s, l) {
  n <- length(DATA); W <- rep(NA, times = n); D <- rep(NA, times = n)
  s0 <- rep(0, times = n); m <- rep(0, times = s); y <- 0
  Time <- seq(1:n); forecast <- 0; shat <- rep(0, times = n)
  seaseffect <- 0; q <- floor(s / 2);

  if (s %% 2 == 0)
    for (i in (q + 1):(n - q)) {
      W[i] <- (0.5 * (DATA[i - q] + DATA[i + q]) + sum(DATA[(i - q + 1):(i + q - 1)])) / s }
  else
    for (i in (q + 1):(n - q)) {
      W[i] <- mean(DATA[(i - q):(i + q)]) }

  for (i in (q + 1):(n - q)) {
    D[i] <- DATA[i] / W[i]
    if (i %% s == 0) {
      s0[s] <- s0[s] + D[i]
      m[s] <- m[s] + 1 }
    else for (j in 1:(s - 1)) {
      if (i %% s == j) {
        s0[j] <- s0[j] + D[i]
        m[j] <- m[j] + 1 }
    }
  }
  for (j in 1:s) {
    seaseffect[j] <- (s0[j] / m[j]) * (s / sum(s0[1:s] / m)) }

  shat <- rep(seaseffect, length.out = n)
  y <- DATA / shat
  betahat <- lm(y ~ Time)$coefficients
  yhat <- predict(lm(y ~ Time))
  forecast <- as.vector(yhat * shat)
  error <- as.vector(DATA - forecast)
  FORECAST <- 0
  for (i in 1:l) {
    if ((n + i) %% s == 0) {
      FORECAST[i] <- (betahat[1] + betahat[2] * (n + i)) * seaseffect[s]
    }
    else FORECAST[i] <- (betahat[1] + betahat[2] * (n + i)) * seaseffect[(n + i) %% s]
  }
  DECMU <- data.frame(DATA, W, D, shat, y, yhat, forecast, error, seaseffect, betahat)
  return(DECMU)
}

# =============
# ===============================================

# Dataset 1
# Portray the plot of Dataset 1.
plot(dataset1, ylab = 'Dataset 1', type = 'o', pch = 19, lwd = 3, xlim = c(1, 11.5))

# Portray the regression line.
regression_line <- lm(dataset1 ~ time(dataset1))
abline(regression_line, lwd = 3, col = 6)
segments(1, 950, 12, 990, lwd = 3, col = 4)
segments(1, 250, 12, 300, lwd = 3, col = 4)

# We choose the best frequency for Dataset1.
freq1 <- 9

# Prediction
TRSE <- NULL
TRSE <- DEC.ADD(dataset1, freq1, 1)
TRSE

# Portray prediction points.
n <- length(dataset1)
points(time(dataset1), TRSE$forecast[1:n], col = 3, type = "b", lwd = 3, pch = 19)

betahat <- TRSE$betahat[1:2]
seaseffect <- TRSE$seaseffect

xhatnltrs <- function(l) {
  FORECAST <- 0
  if ((n + l) %% freq1 == 0) {
    FORECAST <- (betahat[1] + betahat[2] * (n + l)) + seaseffect[freq1]
  }
  else
    FORECAST <- (betahat[1] + betahat[2] * (n + l)) + seaseffect[(n + l) %% freq1]
  return(FORECAST)
}

last_time <- max(time(dataset1))
points(last_time, xhatnltrs(1), col = 2)
points(last_time + 1 / freq1, xhatnltrs(2), col = 2)
points(last_time + 2 / freq1, xhatnltrs(3), col = 2)
points(last_time + 3 / freq1, xhatnltrs(4), col = 2)
points(last_time + 4 / freq1, xhatnltrs(5), col = 2)
points(last_time + 5 / freq1, xhatnltrs(6), col = 2)
points(last_time + 6 / freq1, xhatnltrs(7), col = 2)
points(last_time + 7 / freq1, xhatnltrs(8), col = 2)
points(last_time + 8 / freq1, xhatnltrs(9), col = 2)
points(last_time + 9 / freq1, xhatnltrs(10), col = 2)
lines(seq(last_time, last_time + 9 / freq1, 1 / freq1), lapply(1:10, xhatnltrs), col = 2, lwd = 3)

# Prediction details
par(mfrow = c(1, 1))
plot(TRSE$shat, ylab = 'Seasonal Efect', type = 'o', pch = 19, lwd = 3)
plot(TRSE$yhat, ylab = 'Trend', type = 'o', pch = 19, lwd = 3)
plot(TRSE$error, ylab = 'Errors', type = 'o', pch = 19, lwd = 3)
abline(a = 0, b = 0, col = 'blue', lwd = 3)
par(mfrow = c(1, 1))
# ==============================================================================================

# Dataset 2
# Portray the plot of Dataset 2.
plot(dataset2, ylab = 'Dataset 2', type = 'o', pch = 19, lwd =, xlim = c(1, 9))

# Portray the regression line.
regression_line <- lm(dataset2 ~ time(dataset2))
abline(regression_line, lwd = 3, col = 6)
segments(1, 10200, 9.5, 10200, lwd = 3, col = 4)
segments(1, 1700, 10, 1700, lwd = 3, col = 4)

# We choose the best frequency for Dataset2.
#freq2 <- 6

# Prediction
TRSE <- NULL
TRSE <- DEC.ADD(dataset2, freq2, 1)
TRSE

# Portray prediction points.
n <- length(dataset2)
points(time(dataset2), TRSE$forecast[1:n], col = 3, type = "b", lwd = 3, pch = 19)

betahat <- TRSE$betahat[1:2]
seaseffect <- TRSE$seaseffect

xhatnltrs <- function(l) { FORECAST <- 0
  if ((n + l) %% freq2 == 0) {
    FORECAST <- (betahat[1] + betahat[2] * (n + l)) + seaseffect[freq2]
  }
  else
    FORECAST <- (betahat[1] + betahat[2] * (n + l)) + seaseffect[(n + l) %% freq2]
  return(FORECAST)
}

last_time <- max(time(dataset2))
points(last_time, xhatnltrs(1), col = 2)
points(last_time + 1 / freq2, xhatnltrs(2), col = 2)
points(last_time + 2 / freq2, xhatnltrs(3), col = 2)
points(last_time + 3 / freq2, xhatnltrs(4), col = 2)
points(last_time + 4 / freq2, xhatnltrs(5), col = 2)
points(last_time + 5 / freq2, xhatnltrs(6), col = 2)
points(last_time + 6 / freq2, xhatnltrs(7), col = 2)
lines(seq(last_time, last_time + 6 / freq2, 1 / freq2), lapply(1:7, xhatnltrs), col = 2, lwd = 3)

# Prediction details
par(mfrow = c(1, 1))
plot(TRSE$shat, ylab = 'Seasonal Efect', type = 'o', pch = 19, lwd = 3)
plot(TRSE$yhat, ylab = 'Trend', type = 'o', pch = 19, lwd = 3)
plot(TRSE$error, ylab = 'Errors', type = 'o', pch = 19, lwd = 3)
abline(a = 0, b = 0, col = 'blue', lwd = 3)
par(mfrow = c(1, 1))
# ==============================================================================================

# Dataset 3
# Portray the plot of Dataset 3.
plot(dataset3, ylab = 'Dataset 3', type = 'o', pch = 19, lwd = 3, xlim = c(1, 19))

# Portray the regression line.
regression_line <- lm(dataset3 ~ time(dataset3))
abline(regression_line, lwd = 3, col = 6)
segments(1, 80, 18, 400, lwd = 3, col = 4)
segments(1, -80, 18, -475, lwd = 3, col = 4)

# We choose the best frequency for Dataset3.
#freq3 <- 12

# Prediction
TRSE <- NULL
TRSE <- DEC.MULT(dataset3, freq3, 1)
TRSE

# Portray prediction points.
n <- length(dataset3)
points(time(dataset3), TRSE$forecast[1:n], col = 3, type = "b", lwd = 3, pch = 19)

betahat <- TRSE$betahat[1:2]
seaseffect <- TRSE$seaseffect

xhatnltrsp <- function(l) {
  FORECAST <- 0
  if ((n + l) %% freq3 == 0) {
    FORECAST <- (betahat[1] + betahat[2] * (n + l)) * seaseffect[freq3]
  }
  else
    FORECAST <- (betahat[1] + betahat[2] * (n + l)) * seaseffect[(n + l) %% freq3]
  return(FORECAST)
}

last_time <- max(time(dataset3))
points(last_time, xhatnltrsp(1), col = 2)
points(last_time + 1 / freq3, xhatnltrsp(2), col = 2)
points(last_time + 2 / freq3, xhatnltrsp(3), col = 2)
points(last_time + 3 / freq3, xhatnltrsp(4), col = 2)
points(last_time + 4 / freq3, xhatnltrsp(5), col = 2)
points(last_time + 5 / freq3, xhatnltrsp(6), col = 2)
points(last_time + 6 / freq3, xhatnltrsp(7), col = 2)
points(last_time + 7 / freq3, xhatnltrsp(8), col = 2)
points(last_time + 8 / freq3, xhatnltrsp(9), col = 2)
points(last_time + 9 / freq3, xhatnltrsp(10), col = 2)
points(last_time + 10 / freq3, xhatnltrsp(11), col = 2)
points(last_time + 11 / freq3, xhatnltrsp(12), col = 2)
points(last_time + 12 / freq3, xhatnltrsp(13), col = 2)
lines(seq(last_time, last_time + 12 / freq3, 1 / freq3), lapply(1:13, xhatnltrsp), col = 2, lwd = 3)

# Prediction details
par(mfrow = c(1, 1))
plot(TRSE$shat, ylab = 'Seasonal Efect', type = 'o', pch = 19, lwd = 3)
plot(TRSE$yhat, ylab = 'Trend', type = 'o', pch = 19, lwd = 3)
plot(TRSE$error, ylab = 'Errors', type = 'o', pch = 19, lwd = 3)
abline(a = 0, b = 0, col = 'blue', lwd = 3)
par(mfrow = c(1, 1))
# ==============================================================================================
# The End :)