
#' load CDNow data
cdnow <- cdnow.sample()
elog  <- cdnow$elog
cbs   <- cdnow$cbs

x <- readline("Estimate Models via MLE (press Enter)")

# NBD
(params.nbd <- nbd.EstimateParameters(cbs))

# Pareto/NBD (from BTYD package)
(params.pnbd <- BTYD::pnbd.EstimateParameters(cbs))

# BG/NBD (from BTYD package)
(params.bgnbd <- BTYD::bgnbd.EstimateParameters(cbs))

# MBG/CNBD-k
(params.mbgcnbd <- mbgcnbd.EstimateParameters(cbs))
# -> MBG/CNBD-k is identical to MBG/NBD, as no regularity is detected, hence k=1


x <- readline("Compare Log-Likelihoods (press Enter)")

rbind(`NBD` = nbd.cbs.LL(params.nbd, cbs), 
      `Pareto/NBD` = BTYD::pnbd.cbs.LL(params.pnbd, cbs), 
      `BG/NBD` = BTYD::bgnbd.cbs.LL(params.bgnbd, cbs), 
      `MBG/NBD` = mbgnbd.cbs.LL(params.mbgnbd, cbs), 
      `MBG/CNBD-k` = mbgcnbd.cbs.LL(params.mbgcnbd, cbs))
# -> MBG/NBD provides best fit according to LL


x <- readline("Estimate Transactions in Holdout period (press Enter)")

cbs$xstar.nbd <- nbd.ConditionalExpectedTransactions(params.nbd, cbs$T.star, cbs$x, cbs$T.cal)
cbs$xstar.pnbd <- BTYD::pnbd.ConditionalExpectedTransactions(params.pnbd, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)
cbs$xstar.bgnbd <- BTYD::bgnbd.ConditionalExpectedTransactions(params.bgnbd, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)
cbs$xstar.mbgnbd <- mbgnbd.ConditionalExpectedTransactions(params.mbgnbd, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)
cbs$xstar.mbgcnbd <- mbgcnbd.ConditionalExpectedTransactions(params.mbgcnbd, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)


x <- readline("Estimate P(alive) (press Enter)")

cbs$palive.nbd <- 1
cbs$palive.pnbd <- BTYD::pnbd.PAlive(params = params.pnbd, cbs$x, cbs$t.x, cbs$T.cal)
cbs$palive.bgnbd <- BTYD::bgnbd.PAlive(params = params.bgnbd, cbs$x, cbs$t.x, cbs$T.cal)
cbs$palive.mbgnbd <- mbgnbd.PAlive(params = params.mbgnbd, cbs$x, cbs$t.x, cbs$T.cal)
cbs$palive.mbgcnbd <- mbgcnbd.PAlive(params = params.mbgcnbd, cbs$x, cbs$t.x, cbs$T.cal)


x <- readline("Compare Forecasting Accuracy (press Enter)")

MAE <- function(a, f) mean(abs(a - f))
RMSE <- function(a, f) sqrt(mean((a - f)^2))
MSLE <- function(a, f) mean(((log(a + 1) - log(f + 1)))^2)
BIAS <- function(a, f) sum(f)/sum(a) - 1
bench <- function(cbs, models) {
  acc <- t(sapply(models, function(model) {
    est <- cbs[[paste0("xstar.", model)]]
    c(MAE(cbs$x.star, est),
      RMSE(cbs$x.star, est), 
      MSLE(cbs$x.star, est), 
      BIAS(cbs$x.star, est))
  }))
  colnames(acc) <- c("MAE", "RMSE", "MSLE", "BIAS")
  round(acc, 3)
}

bench(cbs, c("nbd", "pnbd", "bgnbd", "mbgnbd", "mbgcnbd"))
#' Pareto/NBD and MBG/NBD provide best forecasts


x <- readline("Calculate CLV (press Enter)")

#' calculate average spends per transaction
cbs$sales.avg <- cbs$sales / (cbs$x + 1)

#' Note: in CDNow some customers have sales.avg = 0. We substitute the zeros
#' with the minimum non-zero spend, as the estimation fails otherwise
cbs$sales.avg[cbs$sales.avg == 0] <- min(cbs$sales.avg[cbs$sales.avg > 0])

#' Estimate expected average transaction value based on gamma-gamma spend model
spend.params <- BTYD::spend.EstimateParameters(cbs$sales.avg, cbs$x + 1)
cbs$sales.avg.est <- BTYD::spend.expected.value(spend.params, cbs$sales.avg, cbs$x + 1)

#' Calculate CLV for Pareto/NBD
cbs$sales.pnbd <- cbs$sales.avg.est * cbs$xstar.pnbd
cat("Estimated Sales:", round(sum(cbs$sales.pnbd), 1), "\n") 
cat("Actual Sales:", round(sum(cbs$sales.star), 1), "\n") 


x <- readline("Estimate Pareto/NBD (HB) model via MCMC (press Enter)")

#' draw parameter estimates
params.draws <- pnbd.mcmc.DrawParameters(cbs, mcmc = 1500, burnin = 1500, chains = 4)
params.mcmc <- as.list(summary(params.draws$level_2)$quantiles[, "50%"])
rbind(`Pareto/NBD (MCMC)` = round(unlist(params.mcmc), 3), 
      `Pareto/NBD (MLE)` = round(params.pnbd, 3))
#' -> Parameter Estimates between MLE and MCMC match closely

#' draw future transaction estimates
xstar.draws <- mcmc.DrawFutureTransactions(cbs, params.draws)
cbs$xstar.mcmc <- apply(xstar.draws, 2, mean)

#' estimate P(alive)
cbs$palive.mcmc <- mcmc.PAlive(params.draws)

#' estimate P(active)
cbs$pactive.mcmc <- mcmc.PActive(xstar.draws)

#' plot P(active) diagnostic plot
mcmc.plotPActiveDiagnostic(cbs, xstar.draws)


x <- readline("For a demo of Pareto/GGG see `demo(\"pareto-ggg\")`")
