
# Load CDNow event log
library(BTYD)
data(cdnowElog, package = "BTYD")
elog <- data.frame(t(sapply(2:nrow(cdnowElog), function(i) strsplit(as.character(cdnowElog[i, ]), split = ",")[[1]])))
names(elog) <- c("cust", "sampleid", "date", "cds", "sales")
elog$date <- as.Date(elog$date, "%Y%m%d")

# Transform to CBS (including extra column 'litt')
cbs <- elog2cbs(elog, per = "week", T.cal = as.Date("1997-09-30"))

nrow(cbs)
head(cbs)


x <- readline("Estimate Models via MLE (press Enter)")

# NBD
(params.nbd <- nbd.EstimateParameters(cbs))

# Pareto/NBD (from BTYD package)
(params.pnbd <- BTYD::pnbd.EstimateParameters(cbs))

# Gamma/Gompertz/NBD (params.ggnbd <- ggnbd.EstimateParameters(cbs, trace=10)) # would take several minutes
(params.ggnbd <- c(r = 0.552256, alpha = 10.568448, b = 3.5e-05, s = 0.609494, beta = 0.000417))

# BG/NBD (from BTYD package)
(params.bgnbd <- BTYD::bgnbd.EstimateParameters(cbs))

# MBG/NBD
(params.mbgnbd <- mbgnbd.EstimateParameters(cbs))

# MBG/CNBD-k
(params.mbgcnbd <- mbgcnbd.EstimateParameters(cbs))
# -> MBG/CNBD-k is identical to MBG/NBD, as no regularity is detected, hence k=1


x <- readline("Compare Log-Likelihoods (press Enter)")

rbind(NBD = nbd.cbs.LL(params.nbd, cbs), `Pareto/NBD` = BTYD::pnbd.cbs.LL(params.pnbd, cbs), `GG/NBD` = ggnbd.cbs.LL(params.ggnbd, 
    cbs), `BG/NBD` = BTYD::bgnbd.cbs.LL(params.bgnbd, cbs), `MBG/NBD` = mbgnbd.cbs.LL(params.mbgnbd, cbs), `MBG/CNBD-k` = mbgcnbd.cbs.LL(params.mbgcnbd, 
    cbs))
# -> MBG/NBD provides best fit according to LL



x <- readline("Estimate Transactions in Holdout period (press Enter)")

cbs$nbd <- nbd.ConditionalExpectedTransactions(params.nbd, cbs$T.star, cbs$x, cbs$T.cal)
cbs$pnbd <- BTYD::pnbd.ConditionalExpectedTransactions(params.pnbd, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)
cbs$ggnbd <- ggnbd.ConditionalExpectedTransactions(params.ggnbd, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)
cbs$bgnbd <- BTYD::bgnbd.ConditionalExpectedTransactions(params.bgnbd, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)
cbs$mbgnbd <- mbgnbd.ConditionalExpectedTransactions(params.mbgnbd, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)
cbs$mbgcnbd <- mbgcnbd.ConditionalExpectedTransactions(params.mbgcnbd, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)


x <- readline("Estimate P(alive) (press Enter)")

cbs$palive.nbd <- 1
cbs$palive.pnbd <- BTYD::pnbd.PAlive(params = params.pnbd, cbs$x, cbs$t.x, cbs$T.cal)
cbs$palive.ggnbd <- ggnbd.PAlive(params = params.ggnbd, cbs$x, cbs$t.x, cbs$T.cal)
cbs$palive.bgnbd <- BTYD::bgnbd.PAlive(params = params.bgnbd, cbs$x, cbs$t.x, cbs$T.cal)
cbs$palive.mbgnbd <- mbgnbd.PAlive(params = params.mbgnbd, cbs$x, cbs$t.x, cbs$T.cal)
cbs$palive.mbgcnbd <- mbgcnbd.PAlive(params = params.mbgcnbd, cbs$x, cbs$t.x, cbs$T.cal)


x <- readline("Compare Forecasting Accuracy (press Enter)")

MAPE <- function(a, f) sum(abs(a - f)/sum(a))
RMSE <- function(a, f) sqrt(mean((a - f)^2))
MSLE <- function(a, f) mean(((log(a + 1) - log(f + 1)))^2)
BIAS <- function(a, f) sum(f)/sum(a) - 1
bench <- function(cbs, models) {
    acc <- t(sapply(models, function(model) c(MAPE(cbs$x.star, cbs[, model]), RMSE(cbs$x.star, cbs[, model]), 
        MSLE(cbs$x.star, cbs[, model]), BIAS(cbs$x.star, cbs[, model]))))
    colnames(acc) <- c("MAPE", "RMSE", "MSLE", "BIAS")
    round(acc, 3)
}

bench(cbs, c("nbd", "pnbd", "ggnbd", "bgnbd", "mbgnbd", "mbgcnbd"))


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
