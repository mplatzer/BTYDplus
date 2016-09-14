
x <- readline("Estimate Pareto/GGG for for CDNow (press Enter)")

#' Load subset of CDNow data
#' Note: we only use a small set of customers, and a short MCMC chain to keep 
#' this demo finish in reasonable time. Main purpose of this code is to
#' demonstrate the basic workflow of running Pareto/GGG on top of your data.
cdnow <- cdnow.sample()
elog  <- cdnow$elog
cbs   <- cdnow$cbs
set.seed(1234)
cust.subset <- sample(cbs$cust, size = 500)
elog <- elog[elog$cust %in% cust.subset, ]
cbs  <- cbs[cbs$cust %in% cust.subset, ]

#' draw Pareto/GGG parameters and report median estimates
param.draws <- pggg.mcmc.DrawParameters(
                 cal.cbs = cbs, 
                 mcmc = 1500, burnin = 500, chains = 2, thin = 50, 
                 mc.cores = 1)
round(summary(param.draws$level_2)$quantiles[, "50%"], 2)

#' estimate distribution for future transactions
xstar.draws <- mcmc.DrawFutureTransactions(
                 cal.cbs     = cbs, 
                 draws       = param.draws, 
                 T.star      = cbs$T.star)
cbs$xstar.pggg <- apply(xstar.draws, 2, mean)

#' calculate P(alive)
cbs$palive <- mcmc.PAlive(param.draws)

#' calculate P(active)
cbs$pactive <- apply(xstar.draws, 2, function(x) mean(x > 0))

#mcmc.plotPActiveDiagnostic()
#mcmc.PlotTrackingInc()
#mcmc.PlotFrequencyInCalibration()


x <- readline("Estimate Pareto/GGG and Pareto/NBD (HB) for simulated data (press Enter)")

#' generate artificial Pareto/GGG data
params <- list(t = 4.5, gamma = 1.5, 
               r = 5, alpha = 30, 
               s = 0.8, beta = 12)
set.seed(1234)
data <- pggg.GenerateData(
          n = 1000, 
          T.cal = 21, 
          T.star = 21, 
          params = params, 
          return.elog = TRUE)
cbs <- data$cbs
elog <- data$elog

#' estimate Pareto/NBD (HB) and Pareto/GGG
pnbd.draws <- pnbd.mcmc.DrawParameters(cbs, mcmc = 1500, burnin = 500, chains = 2, thin = 50, mc.cores = 1)
pggg.draws <- pggg.mcmc.DrawParameters(cbs, mcmc = 1500, burnin = 500, chains = 2, thin = 50, mc.cores = 1)
rbind(`Actual`     = params, 
      `Pareto/GGG` = summary(pggg.draws$level_2)$quantiles[, "50%"],
      `Pareto/NBD` = c(t=NA, gamma=NA, summary(pnbd.draws$level_2)$quantiles[, "50%"]))

#' plot estimated parameter distributions
plot(param.draws$level_2, density = TRUE, trace = FALSE)
pggg.plotRegularityRateHeterogeneity(pggg.draws)

#' check MCMC convergence
plot(param.draws$level_2, density = FALSE, trace = TRUE)
coda::effectiveSize(pggg.draws$level_2)

#' draw future transaction
pnbd.xstar.draws <- mcmc.DrawFutureTransactions(cbs, pnbd.draws, T.star = cbs$T.star, sample_size = 400)
pggg.xstar.draws <- mcmc.DrawFutureTransactions(cbs, pggg.draws, T.star = cbs$T.star, sample_size = 400)

#' calculate mean over future transaction draws for each customer
cbs$xstar.pnbd <- apply(pnbd.xstar.draws, 2, mean)
cbs$xstar.pggg <- apply(pggg.xstar.draws, 2, mean)

#' compare forecast accuracy to Pareto/NBD
(mae <- c(`Pareto/GGG` = mean(abs(cbs$x.star - cbs$xstar.pggg)), 
          `Pareto/NBD` = mean(abs(cbs$x.star - cbs$xstar.pnbd))))
(lift <- 1 - mae[1]/mae[2])
#' -> 18% lift in customer-level accuracy when taking regularity into account

# P(active) diagnostic plot
nil <- mcmc.plotPActiveDiagnostic(cbs, pggg.xstar.draws)

#' Compare incremental transactions in holdout to Pareto/NBD
inc <- elog2inc(elog, by = 1)
T.tot <- max(cbs$T.cal + cbs$T.star)
op <- par(mfrow = c(1, 2))
nil <- mcmc.PlotTrackingInc(pggg.draws, cbs$T.cal, T.tot, inc, title = "Pareto/GGG", xlab="Days")
nil <- mcmc.PlotTrackingInc(pnbd.draws, cbs$T.cal, T.tot, inc, title = "Pareto/NBD", xlab="Days")
par(op)
