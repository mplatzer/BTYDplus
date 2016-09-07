
x <- readline("Estimate BG/CNBD-k for CDNow data (press Enter)")

#' load CDNow data
cdnow <- cdnow.sample()
elog  <- cdnow$elog
cbs   <- cdnow$cbs

#' estimate BG/CNBD-k parameters
(params <- bgcnbd.EstimateParameters(cbs))
#' k=1 -> no regularity detected for CDNow

#' estimate transactions in holdout period
cbs$xstar.est <- bgcnbd.ConditionalExpectedTransactions(params, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)

#' estimate probability of being alive at end of calibration period
cbs$palive <- bgcnbd.PAlive(params, cbs$x, cbs$t.x, cbs$T.cal)

#' plot estimated incremental transactions
inc <- elog2inc(elog)
nil <- bgcnbd.PlotTrackingInc(params, cbs$T.cal, T.tot = 78, inc)


x <- readline("Estimate BG/CNBD-k for self generated data (press Enter)")

#' generate artificial BG/CNBD-3 data
n <- 4000 # number of customers
T.cal <- round(runif(n, 24, 32))  # 24-32 weeks of calibration period
T.star <- 32  # 32 weeks of hold-out period
params <- c(k = 3, r = 0.85, alpha = 1.45, a = 0.79, b = 2.42)
#' regularity in interpurchase-times (Erlang-k)
#' purchase frequency lambda_i ~ Gamma(r, alpha) 
#' dropout probability p_i ~ Beta(a, b)
data <- bgcnbd.GenerateData(n, T.cal, T.star, params, return.elog = TRUE)
cbs <- data$cbs  # CBS summary - one record per customer
elog <- data$elog  # Event log - one row per event/purchase

#' estimate regularity from event log
(k.est <- estimateRegularity(elog))
#' -> Wheat-Morrison estimator correctly detects Erlang-3

#' estimate parameters, and compare to true parameters
params_k <- bgcnbd.EstimateParameters(cbs[, c("x", "t.x", "T.cal", "litt")])
params_1 <- BTYD::bgnbd.EstimateParameters(cbs[, c("x", "t.x", "T.cal", "litt")])
rbind(actual = params, 
      `bg/cnbd-k` = round(params_k, 2), 
      `bg/nbd` = c(1, round(params_1, 2)))
#' -> underlying parameters are successfully identified via Maximum Likelihood Estimation

#' estimate future transactions in holdout-period
cbs$xstar_k <- bgcnbd.ConditionalExpectedTransactions(params_k, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)
cbs$xstar_1 <- BTYD::bgnbd.ConditionalExpectedTransactions(params_1, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)

#' estimate P(alive)
cbs$palive_k <- bgcnbd.PAlive(params_k, cbs$x, cbs$t.x, cbs$T.cal)
cbs$palive_1 <- BTYD::bgnbd.PAlive(params_1, cbs$x, cbs$t.x, cbs$T.cal)

#' compare forecast accuracy to bg/nbd
rbind(`bg/cnbd-k` = mean(abs(cbs$x.star - cbs$xstar_k)), 
      `bg/nbd` = mean(abs(cbs$x.star - cbs$xstar_1)))
#' -> BG/CNBD-k forecast more accurate than BG/NBD

#' compare forecast bias to bg/nbd
rbind(`bg/cnbd-k` = 1 - sum(cbs$xstar_k)/sum(cbs$x.star), 
      `bg/nbd` = 1 - sum(cbs$xstar_1)/sum(cbs$x.star))
#' -> Unbiased estimate for BG/CNBD-k


x <- readline("Compare aggregate fit in calibration to BG/NBD (press Enter)")

op <- par(mfrow = c(1, 2))
nil <- bgcnbd.PlotFrequencyInCalibration(params_k, cbs, censor = 7)
nil <- BTYD::bgnbd.PlotFrequencyInCalibration(params_1, cbs, censor = 7)
par(op)


x <- readline("Compare incremental transactions in holdout to BG/NBD (press Enter)")

op <- par(mfrow = c(1, 2))
inc <- elog2inc(elog, by = 1)
T.tot <- max(cbs$T.cal + cbs$T.star)
inc_k <- bgcnbd.PlotTrackingInc(params_k, cbs$T.cal, T.tot, inc)
inc_1 <- BTYD::bgnbd.PlotTrackingInc(params_1, cbs$T.cal, T.tot, inc)
par(op)


x <- readline("Compare estimated distributions (press Enter)")

#' compare estimated with actual distributions in lambda & churn probability
par(mfrow = c(2, 1), mar = c(2, 1, 2, 1))
xlim <- 1.5
x <- seq(0, xlim, len = 1000)[-1]
y <- dgamma(x, shape = params_k[2], rate = params_k[3] * params_k[1])
plot(x, y, typ = "l", col = "black", lwd = 2, main = "Heterogeneity in Intertransaction Times", ylab = "", xlab = "", 
  axes = FALSE, xlim = c(0, xlim))
lines(x, dgamma(x, shape = params_k[2], rate = params_k[3] * params_k[1]), col = "red", lwd = 2)
lines(x, dgamma(x, shape = params_1[1], rate = params_1[2]), col = "blue", lwd = 2)
axis(1, pos = 0, labels = round(1/(xlim * (0:10/10)), 1), at = xlim * (0:10/10))
legend(xlim * 0.6, max(y) * 0.9, c("Actual", "BG/CNBD-k", "BG/NBD"), col = c("black", "red", "blue"), pch = 15, 
  bty = "n")

xlim <- 1
x <- seq(0.05, xlim, len = 1000)
y <- dbeta(x, params_k[4], params_k[5])
plot(x, y, typ = "l", col = "black", lwd = 2, main = "Heterogeneity in Churn Probability", ylab = "", xlab = "", 
  axes = FALSE, xlim = c(0, xlim))
lines(x, dbeta(x, params_k[4], params_k[5]), col = "red", lwd = 2)
lines(x, dbeta(x, params_1[3], params_1[4]), col = "blue", lwd = 2)
axis(1, pos = 0)
legend(xlim * 0.6, max(y) * 0.9, c("Actual", "BG/CNBD-k", "BG/NBD"), col = c("black", "red", "blue"), pch = 15, 
  bty = "n") 
