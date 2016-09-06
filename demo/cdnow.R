
# Load CDNow event log
library(BTYD)
data(cdnowElog, package = "BTYD", envir = environment())
elog <- data.frame(t(sapply(2:nrow(cdnowElog), 
                            function(i) strsplit(as.character(cdnowElog[i, ]), split = ",")[[1]])),
                   stringsAsFactors = FALSE)
names(elog) <- c("cust", "sampleid", "date", "cds", "sales")
elog$date <- as.Date(elog$date, "%Y%m%d")
elog$sales <- as.numeric(elog$sales)

# Transform to CBS (including extra summary statistic 'litt' for estimating regularity)
cutoff <- as.Date("1997-09-30")
cbs <- elog2cbs(elog, per = "week", T.cal = cutoff)

nrow(cbs)
head(cbs)


x <- readline("Estimate Models via MLE (press Enter)")

# NBD
(params.nbd <- nbd.EstimateParameters(cbs))

# Pareto/NBD (from BTYD package)
(params.pnbd <- BTYD::pnbd.EstimateParameters(cbs))

# BG/NBD (from BTYD package)
(params.bgnbd <- BTYD::bgnbd.EstimateParameters(cbs))

# MBG/NBD
(params.mbgnbd <- mbgnbd.EstimateParameters(cbs))

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
    est <- cbs[[paste0('xstar.', model)]]
    c(MAE(cbs$x.star, est),
      RMSE(cbs$x.star, est), 
      MSLE(cbs$x.star, est), 
      BIAS(cbs$x.star, est))
  }))
  colnames(acc) <- c("MAE", "RMSE", "MSLE", "BIAS")
  round(acc, 3)
}

bench(cbs, c("nbd", "pnbd", "bgnbd", "mbgnbd", "mbgcnbd"))


x <- readline("Calculate CLV (press Enter)")

# calculate average spends per transaction
cbs$sales.avg      <- cbs$sales / (cbs$x + 1)

# Note: we substitute customers with 0 spend with minimum non-zero spend, as the estimation fails otherwise
cbs$sales.avg[cbs$sales.avg == 0] <- min(cbs$sales.avg[cbs$sales.avg > 0])

# Estimate expected average transaction value based on gamma-gamma spend model
spend.params <- BTYD::spend.EstimateParameters(cbs$sales.avg, cbs$x + 1)
cbs$sales.avg.est <- spend.expected.value(spend.params, cbs$sales.avg, cbs$x + 1)

# Calculate CLV for Pareto/NBD
cbs$sales.pnbd <- cbs$sales.avg.est * cbs$xstar.pnbd
sprintf('Estimated Sales: %.1f, Actual Sales: %.1f',
        sum(cbs$sales.pnbd), sum(cbs$sales.star))


x <- readline("Estimate Pareto/NBD model via MCMC (press Enter)")

set.seed(1)
pnbd.draws <- pnbd.mcmc.DrawParameters(cbs, mcmc = 1500, burnin = 1500, chains = 4)

summary(pnbd.draws$level_2)
plot(pnbd.draws$level_2)

params.pnbd.mcmc <- as.list(summary(pnbd.draws$level_2)$quantiles[, "50%"])
rbind(MCMC = round(unlist(params.pnbd.mcmc), 3), MLE = round(params.pnbd, 3))
# -> Parameter Estimates between MLE and MCMC match closely


x <- readline("Estimate Pareto/GGG model via MCMC - this will take several minutes (press Enter)")

# Note: For keeping the runtime of this demo short, we limit the MCMC to 500 iterations for both chains. In
# fact, the chains should be run longer to ensure convergence, and collect enough samples

set.seed(1)

pggg.draws <- pggg.mcmc.DrawParameters(cbs, mcmc = 500, burnin = 100, chains = 2, thin = 10)

plot(pggg.draws$level_2, density = FALSE)
plot(pggg.draws$level_2, trace = FALSE)

coda::gelman.diag(pggg.draws$level_2)
# -> MCMC chains have not converged yet

(summary(pggg.draws$level_2)$quantiles[, "50%"])

pggg.mcmc.plotRegularityRateHeterogeneity(pggg.draws)
# -> very narrow distribution around k=1; CDNow customers purchases indeed seem to follow Poisson process

round(coda::effectiveSize(pggg.draws$level_2))
# -> effective sample size are small for such a short chain


x <- readline("Estimate Future Transactions & P(active) for MCMC models - this will take several minutes (press Enter)")

# draw future transaction
pnbd.xstar <- mcmc.DrawFutureTransactions(cbs, pnbd.draws, T.star = cbs$T.star)
pggg.xstar <- mcmc.DrawFutureTransactions(cbs, pggg.draws, T.star = cbs$T.star)

# calculate mean over future transaction draws for each customer
cbs$pnbd.mcmc <- apply(pnbd.xstar, 2, mean)
cbs$pggg.mcmc <- apply(pggg.xstar, 2, mean)

# forecasting accuracy
bench(cbs, c("pnbd", "pnbd.mcmc", "pggg.mcmc"))

# calculate P(active)
cbs$pactive.pnbd.mcmc <- apply(pnbd.xstar, 2, function(x) mean(x > 0))
cbs$pactive.pggg.mcmc <- apply(pggg.xstar, 2, function(x) mean(x > 0))

# calculate P(alive)
cbs$palive.pnbd.mcmc <- mcmc.PAlive(cbs, pnbd.draws)
cbs$palive.pggg.mcmc <- mcmc.PAlive(cbs, pggg.draws) 

