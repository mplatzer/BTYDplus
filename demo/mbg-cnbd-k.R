
#' Simulate artificial MBG/CNBD-3 data.
set.seed(1234)
n      <- 4000                         # number of customers
T.cal  <- round(runif(n, 24, 32)/4)*4  # 24-32 weeks of calibration period
T.star <- 32                           # 32 weeks of hold-out period
params <- c(k = 3,                  # regularity in interpurchase-times (Erlang-k)
            r = 0.85, alpha = 1.45, # purchase frequency lambda_i ~ Gamma(r, alpha) 
            a = 0.79, b = 2.42)     # dropout probability p_i ~ Beta(a, b)
data <- mbgcnbd.GenerateData(n, T.cal, T.star, params, return.elog = TRUE)

cbs  <- data$cbs   # CBS summary - one record per customer
head(cbs)

elog <- data$elog  # Event log - one row per event/purchase
head(elog)


x <- readline("Estimate regularity via Wheat/Morrison estimator (press Enter)")

(k.est <- estimateRegularity(elog))
#' -> Wheat-Morrison estimator correctly detects Erlang-3.


x <- readline("Estimate MBG/CNBD-k model (press Enter)")

params.mbgcnbd <- mbgcnbd.EstimateParameters(cbs)
params.pnbd    <- BTYD::pnbd.EstimateParameters(cbs)
rbind(`Actual`     = params, 
      `MBG/CNBD-k` = round(params.mbgcnbd, 2), 
      `Pareto/NBD` = c(1, round(params.pnbd, 2)))
#' -> Underlying parameters are successfully recovered.


x <- readline("Estimate transactions during holdout period (press Enter)")

cbs$xstar.mbgcnbd <- mbgcnbd.ConditionalExpectedTransactions(
                 params = params.mbgcnbd, 
                 T.star = cbs$T.star, 
                 x      = cbs$x, 
                 t.x    = cbs$t.x, 
                 T.cal  = cbs$T.cal)
cbs$xstar.pnbd <- BTYD::pnbd.ConditionalExpectedTransactions(
                 params = params.pnbd, 
                 T.star = cbs$T.star, 
                 x      = cbs$x, 
                 t.x    = cbs$t.x, 
                 T.cal  = cbs$T.cal)

#' compare forecast accuracy to Pareto/NBD
(mae <- c(`MBG/CNBD-k` = mean(abs(cbs$x.star - cbs$xstar.mbgcnbd)), 
         `Pareto/NBD` = mean(abs(cbs$x.star - cbs$xstar.pnbd))))
(lift <- 1 - mae[1]/mae[2])
#' -> 11% lift in customer-level accuracy when taking regularity into account


x <- readline("Estimate probabilty of being still alive at end of calibration (press Enter)")

cbs$palive.mbgcnbd <- mbgcnbd.PAlive(
                  params = params.mbgcnbd, 
                  x      = cbs$x, 
                  t.x    = cbs$t.x, 
                  T.cal  = cbs$T.cal)
cbs$palive.pnbd <- BTYD::pnbd.PAlive(
                  params = params.pnbd, 
                  x      = cbs$x, 
                  t.x    = cbs$t.x, 
                  T.cal  = cbs$T.cal)
rbind(`Actual` = mean(cbs$alive),
      `MBG/CNBD-k` = mean(cbs$palive.mbgcnbd),
      `Pareto/NBD` = mean(cbs$palive.pnbd))


x <- readline("Compare aggregate fit in calibration to Pareto/NBD (press Enter)")

op <- par(mfrow = c(1, 2))
nil <- mbgcnbd.PlotFrequencyInCalibration(params.mbgcnbd, cbs, censor = 7, title = "MBG/CNBD-k")
nil <- BTYD::pnbd.PlotFrequencyInCalibration(params.pnbd, cbs, censor = 7, title = "Pareto/NBD")
par(op)


x <- readline("Compare incremental transactions in holdout to Pareto/NBD (press Enter)")

op <- par(mfrow = c(1, 2))
inc <- elog2inc(elog, by = 1)
T.tot <- max(cbs$T.cal + cbs$T.star)
nil <- mbgcnbd.PlotTrackingInc(params.mbgcnbd, cbs$T.cal, T.tot, inc, title = "MBG/CNBD-k")
nil <- BTYD::pnbd.PlotTrackingInc(params.pnbd, cbs$T.cal, T.tot, inc, title = "Pareto/NBD")
par(op)
