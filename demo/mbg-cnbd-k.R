
x <- readline("Estimate MBG/CNBD-k for CDNow data (press Enter)")

#' load CDNow data
cdnow <- cdnow.sample()
elog  <- cdnow$elog
cbs   <- cdnow$cbs

#' estimate MBG/CNBD-k parameters
(params <- mbgcnbd.EstimateParameters(cbs))
#' k=1 -> no regularity detected for CDNow

#' estimate transactions in holdout period
cbs$xstar.est <- mbgcnbd.ConditionalExpectedTransactions(params, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)

#' estimate probability of being alive at end of calibration period
cbs$palive <- mbgcnbd.PAlive(params, cbs$x, cbs$t.x, cbs$T.cal)

#' plot estimated incremental transactions
inc <- elog2inc(elog)
nil <- mbgcnbd.PlotTrackingInc(params, cbs$T.cal, T.tot = 78, inc)


x <- readline("Estimate MBG/CNBD-k for self generated data (press Enter)")

#' generate artificial MBG/CNBD-3 data
n <- 4000 # number of customers
T.cal <- round(runif(n, 24, 32)/4)*4  # 24-32 weeks of calibration period
T.star <- 32  # 32 weeks of hold-out period
params <- c(k = 3, r = 0.85, alpha = 1.45, a = 0.79, b = 2.42)
#' regularity in interpurchase-times (Erlang-k)
#' purchase frequency lambda_i ~ Gamma(r, alpha) 
#' dropout probability p_i ~ Beta(a, b)
data <- mbgcnbd.GenerateData(n, T.cal, T.star, params, return.elog = TRUE)
cbs <- data$cbs  # CBS summary - one record per customer
elog <- data$elog  # Event log - one row per event/purchase

#' estimate regularity from event log
(k.est <- estimateRegularity(elog))
#' -> Wheat-Morrison estimator correctly detects Erlang-3

#' estimate parameters, and compare to true parameters
params_k <- mbgcnbd.EstimateParameters(cbs[, c("x", "t.x", "T.cal", "litt")])
params_1 <- BTYD::bgnbd.EstimateParameters(cbs[, c("x", "t.x", "T.cal", "litt")])
rbind(actual = params, 
      `MBG/CNBD-k` = round(params_k, 2), 
      `BG/NBD` = c(1, round(params_1, 2)))
#' -> underlying parameters are successfully identified via Maximum Likelihood Estimation

#' estimate future transactions in holdout-period
cbs$xstar_k <- mbgcnbd.ConditionalExpectedTransactions(params_k, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)
cbs$xstar_1 <- BTYD::bgnbd.ConditionalExpectedTransactions(params_1, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)

#' estimate P(alive)
cbs$palive_k <- mbgcnbd.PAlive(params_k, cbs$x, cbs$t.x, cbs$T.cal)
cbs$palive_1 <- BTYD::bgnbd.PAlive(params_1, cbs$x, cbs$t.x, cbs$T.cal)

#' compare forecast accuracy to bg/nbd
rbind(`MBG/CNBD-k` = mean(abs(cbs$x.star - cbs$xstar_k)), 
      `BG/NBD` = mean(abs(cbs$x.star - cbs$xstar_1)))
#' -> MBG/CNBD-k forecast more accurate than MBG/NBD

#' compare forecast bias to bg/nbd
rbind(`MBG/CNBD-k` = 1 - sum(cbs$xstar_k)/sum(cbs$x.star), 
      `BG/NBD` = 1 - sum(cbs$xstar_1)/sum(cbs$x.star))
#' -> Unbiased estimate for MBG/CNBD-k


x <- readline("Compare aggregate fit in calibration to MBG/NBD (press Enter)")

op <- par(mfrow = c(1, 2))
nil <- mbgcnbd.PlotFrequencyInCalibration(params_k, cbs, censor = 7, title = "MBG/CNBD-k")
nil <- BTYD::bgnbd.PlotFrequencyInCalibration(params_1, cbs, censor = 7, title = "BG/NBD")
par(op)


x <- readline("Compare incremental transactions in holdout to MBG/NBD (press Enter)")

op <- par(mfrow = c(1, 2))
inc <- elog2inc(elog, by = 1)
T.tot <- max(cbs$T.cal + cbs$T.star)
inc_k <- mbgcnbd.PlotTrackingInc(params_k, cbs$T.cal, T.tot, inc, title = "MBG/CNBD-k")
inc_1 <- BTYD::bgnbd.PlotTrackingInc(params_1, cbs$T.cal, T.tot, inc, title = "BG/NBD")
par(op)
