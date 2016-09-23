
#' Load transaction records of 2357 CDNow customers.
cdnowElog <- read.csv(system.file("data/cdnowElog.csv", package = "BTYD"),
                      stringsAsFactors = FALSE,
                      col.names = c("cust", "sampleid", "date", "cds", "sales"))
cdnowElog$date <- as.Date(as.character(cdnowElog$date), format = "%Y%m%d")
head(cdnowElog)
range(cdnowElog$date)

#' Convert from event log to customer-by-sufficient-statistic summary.
#' Split into 39 weeks calibration, and 39 weeks holdout period.
cbs <- elog2cbs(cdnowElog, T.cal = "1997-09-30", T.tot = "1998-06-30")
head(cbs)

#' Report some basic summary stats for calibration and holdout period.
summary(cbs[, c("x", "sales", "x.star", "sales.star")])


x <- readline("Estimate NBD model (press Enter)")

# Estimate NBD model parameters.
(params.nbd <- nbd.EstimateParameters(cbs))

# Predict transactions at customer level with NBD model.
cbs$xstar.nbd <- nbd.ConditionalExpectedTransactions(
                   params = params.nbd,
                   T.star = cbs$T.star,
                   x      = cbs$x,
                   T.cal  = cbs$T.cal)

# Estimate total transactions during holdout, based on NBD model.
sum(cbs$xstar.nbd)


x <- readline("Estimate MBG/CNBD-k model (press Enter)")

# Estimate MBG/CNBD-k model parameters.
(params.mbgcnbd <- mbgcnbd.EstimateParameters(cbs))
#' k=1 -> no regularity detected for CDNow

# Predict transactions at customer level with MBG/CNBD-k model.
cbs$xstar.mbgcnbd <- mbgcnbd.ConditionalExpectedTransactions(
                       params = params.mbgcnbd,
                       T.star = cbs$T.star,
                       x      = cbs$x,
                       t.x    = cbs$t.x,
                       T.cal  = cbs$T.cal)

# Estimate total transactions during holdout, based on MBG/CNBD-k model.
sum(cbs$xstar.mbgcnbd)

# Estimate probabilty of being still a customer at end of calibration period.
cbs$palive.mbgcnbd <- mbgcnbd.PAlive(
                        params = params.mbgcnbd,
                        x      = cbs$x,
                        t.x    = cbs$t.x,
                        T.cal  = cbs$T.cal)

# Estimate share of retained customers at end of calibration period.
mean(cbs$palive.mbgcnbd)


x <- readline("Compare Log-Likelihoods of various models (press Enter)")

params.pnbd <- BTYD::pnbd.EstimateParameters(cbs) # estimate Pareto/NBD
params.bgcnbd <- bgcnbd.EstimateParameters(cbs)  # estimate BG/CNBD-k

(ll <- c(`NBD`        = nbd.cbs.LL(params.nbd, cbs),
        `Pareto/NBD` = BTYD::pnbd.cbs.LL(params.pnbd, cbs),
        `BG/CNBD-k`  = bgcnbd.cbs.LL(params.bgcnbd, cbs),
        `MBG/CNBD-k` = mbgcnbd.cbs.LL(params.mbgcnbd, cbs)))
names(which.max(ll))
# -> MBG/CNBD-k provides best fit according to log-likelihood


x <- readline("Plot Frequency in Calibration (press Enter)")

op <- par(mfrow = c(1, 2))
nil <- mbgcnbd.PlotFrequencyInCalibration(params.mbgcnbd, cbs, censor = 7, title = "MBG/CNBD-k")
nil <- BTYD::pnbd.PlotFrequencyInCalibration(params.pnbd, cbs, censor = 7, title = "Pareto/NBD")
par(op)


x <- readline("Plot Incremental Transactions (press Enter)")

inc <- elog2inc(cdnowElog)
op <- par(mfrow = c(1, 2))
nil <- mbgcnbd.PlotTrackingInc(params.mbgcnbd, cbs$T.cal, T.tot = 78, inc, title = "MBG/CNBD-k")
nil <- BTYD::pnbd.PlotTrackingInc(params.pnbd, cbs$T.cal, T.tot = 78, inc, title = "Pareto/NBD")
par(op)


x <- readline("Compare Forecasting Accuracy (press Enter)")

cbs$xstar.pnbd <- BTYD::pnbd.ConditionalExpectedTransactions(
                    params  = params.pnbd,
                    T.star  = cbs$T.star,
                    x       = cbs$x,
                    t.x     = cbs$t.x,
                    T.cal   = cbs$T.cal)
cbs$xstar.bgcnbd <- bgcnbd.ConditionalExpectedTransactions(
                      params  = params.bgcnbd,
                      T.star  = cbs$T.star,
                      x       = cbs$x,
                      t.x     = cbs$t.x,
                      T.cal   = cbs$T.cal)

measures <- c(
  "MAE" = function(a, f) mean(abs(a - f)),
  "MSLE" = function(a, f) mean(((log(a + 1) - log(f + 1)))^2),
  "BIAS" = function(a, f) sum(f)/sum(a) - 1)
models <- c(
  "NBD" = "nbd",
  "Pareto/NBD" = "pnbd",
  "BG/CNBD-k" = "bgcnbd",
  "MBG/CNBD-k" = "mbgcnbd")

sapply(measures, function(measure) {
  sapply(models, function(model) {
    err <- do.call(measure, list(a = cbs$x.star, f = cbs[[paste0("xstar.", model)]]))
    round(err, 3)
  })
})
#' -> Pareto/NBD and MBG/CNBD-k provide best forecasts


x <- readline("Calculate CLV (press Enter)")

#' calculate average spends per transaction
cbs$sales.avg <- cbs$sales / (cbs$x + 1)

#' Note: in CDNow some customers have sales.avg = 0. We substitute the zeros
#' with the minimum non-zero spend, as the estimation fails otherwise
cbs$sales.avg[cbs$sales.avg == 0] <- min(cbs$sales.avg[cbs$sales.avg > 0])

#' Estimate expected average transaction value based on gamma-gamma spend model
spend.params <- BTYD::spend.EstimateParameters(cbs$sales.avg, cbs$x + 1)
cbs$sales.avg.est <- BTYD::spend.expected.value(spend.params, cbs$sales.avg, cbs$x + 1)

#' Estimated total sales during holdout for MBG/CNBD-k
cbs$sales.mbgcnbd <- cbs$sales.avg.est * cbs$xstar.mbgcnbd
c("Estimated Sales" = sum(cbs$sales.mbgcnbd),
  "Actual Sales"    = sum(cbs$sales.star))


x <- readline("For a demo of Pareto/NBD (Abe) with CDNow see `demo(\"pareto-abe\")`.")
