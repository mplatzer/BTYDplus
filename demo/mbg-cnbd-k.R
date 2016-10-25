#' Load transaction records of 1525 grocery customers.
data("groceryElog", envir = environment())
head(groceryElog)
range(groceryElog$date)

#' Convert from event log to customer-by-sufficient-statistic summary.
#' Split into 52 weeks calibration, and 52 weeks holdout period.
cbs <- elog2cbs(groceryElog, T.cal = "2006-12-31", T.tot = "2007-12-30")
head(cbs)


x <- readline("Estimate regularity via Wheat/Morrison estimator (press Enter)")

(k.est <- estimateRegularity(groceryElog))
#' -> Wheat-Morrison estimator detects Erlang-2.

#' Plot Timing Patterns of a few sampled customers
plotTimingPatterns(groceryElog, T.cal = "2006-12-31")


x <- readline("Estimate MBG/CNBD-k model (press Enter)")

# Estimate MBG/CNBD-k model parameters.
(params.mbgcnbd <- mbgcnbd.EstimateParameters(cbs))
#' k=2 -> regularity also detected via MBG/CNBD-k model

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


x <- readline("Compare log-likelihoods of various models (press Enter)")

params.nbd <- nbd.EstimateParameters(cbs)            # estimate NBD
params.pnbd <- BTYD::pnbd.EstimateParameters(cbs)    # estimate Pareto/NBD
params.bgnbd <- BTYD::bgnbd.EstimateParameters(cbs)  # estimate BG/NBD
params.mbgnbd <- mbgnbd.EstimateParameters(cbs)      # estimate MBG/NBD

(ll <- c(`NBD`        = nbd.cbs.LL(params.nbd, cbs),
         `Pareto/NBD` = BTYD::pnbd.cbs.LL(params.pnbd, cbs),
         `BG/NBD`     = BTYD::bgnbd.cbs.LL(params.bgnbd, cbs),
         `MBG/NBD`    = mbgcnbd.cbs.LL(params.mbgnbd, cbs),
         `MBG/CNBD-k` = mbgcnbd.cbs.LL(params.mbgcnbd, cbs)))
names(which.max(ll))
# -> MBG/CNBD-k provides best data fit according to log-likelihood


x <- readline("Compare forecast accuracies of various models (press Enter)")

cbs$xstar.nbd <- nbd.ConditionalExpectedTransactions(
  params.nbd, cbs$T.star, cbs$x, cbs$T.cal)
cbs$xstar.pnbd <- BTYD::pnbd.ConditionalExpectedTransactions(
  params.pnbd, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)
cbs$xstar.bgnbd <- BTYD::bgnbd.ConditionalExpectedTransactions(
  params.bgnbd, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)
cbs$xstar.mbgnbd <- mbgcnbd.ConditionalExpectedTransactions(
  params.mbgnbd, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)

measures <- c(
  "MAE" = function(a, f) mean(abs(a - f)),
  "MSLE" = function(a, f) mean(((log(a + 1) - log(f + 1)))^2),
  "BIAS" = function(a, f) sum(f)/sum(a) - 1)
models <- c(
  "NBD" = "nbd",
  "Pareto/NBD" = "pnbd",
  "BG/NBD"     = "bgnbd",
  "MBG/NBD"    = "mbgnbd",
  "MBG/CNBD-k" = "mbgcnbd")

sapply(measures, function(measure) {
  sapply(models, function(model) {
    err <- do.call(measure, list(a = cbs$x.star, f = cbs[[paste0("xstar.", model)]]))
    round(err, 3)
  })
})
#' -> MBG/CNBD-k provides best customer-level forecast accuracy


x <- readline("Plot Frequency in Calibration (press Enter)")

op <- par(mfrow = c(1, 2))
nil <- mbgcnbd.PlotFrequencyInCalibration(params.mbgcnbd, cbs, censor = 7, title = "MBG/CNBD-k")
nil <- BTYD::pnbd.PlotFrequencyInCalibration(params.pnbd, cbs, censor = 7, title = "Pareto/NBD")
par(op)


x <- readline("Plot Incremental Transactions (press Enter)")

inc <- elog2inc(groceryElog)
T.tot <- max(cbs$T.cal+cbs$T.star)
op <- par(mfrow = c(1, 2))
nil <- mbgcnbd.PlotTrackingInc(params.mbgcnbd, cbs$T.cal, T.tot = T.tot, inc, title = "MBG/CNBD-k")
nil <- BTYD::pnbd.PlotTrackingInc(params.pnbd, cbs$T.cal, T.tot = T.tot, inc, title = "Pareto/NBD")
par(op)
