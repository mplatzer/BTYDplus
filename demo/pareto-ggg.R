
#' Note: we only use a small set of customers, and a short MCMC chain to keep 
#' this demo finish in reasonable time. Main purpose of this code is to
#' demonstrate the basic workflow of running Pareto/GGG on top of your data.

x <- readline("Estimate Pareto/GGG for subset of CDNow data (press Enter)")

#' load subset of CDNow data
cdnow <- cdnow.sample()
elog  <- cdnow$elog
cbs   <- cdnow$cbs
cust.subset <- sample(cbs$cust, size = 100)
elog <- elog[elog$cust %in% cust.subset, ]
cbs  <- cbs[cbs$cust %in% cust.subset, ]

#' draw Pareto/GGG parameters
pggg.draws <- pggg.mcmc.DrawParameters(cbs, mcmc = 1500, burnin = 500, chains = 2, thin = 50)
round(params <- summary(pggg.draws$level_2)$quantiles[, "50%"], 2)

#' check MCMC chain convergence and parameter distributions
plot(pggg.draws$level_2, density = T, trace = F)
plot(pggg.draws$level_2, density = F, trace = T)

#' estimate distribution for future transactions
xstar.pggg.draws <- mcmc.DrawFutureTransactions(cbs, pggg.draws, T.star = cbs$T.star)
cbs$xstar.pggg <- apply(xstar.pggg.draws, 2, mean)

#' calculate P(active)
cbs$pactive <- apply(xstar.pggg.draws, 2, function(x) mean(x > 0))

#' calculate P(alive)
cbs$palive <- mcmc.PAlive(cbs, pggg.draws)


x <- readline("Estimate Pareto/GGG and Pareto/NBD (HB) for self generated data (press Enter)")

#' generate artificial Pareto/GGG data
params <- list(t = 4.5, gamma = 1.5, r = 5, alpha = 10, s = 0.8, beta = 12)
data <- pggg.GenerateData(n = 1000, T.cal = 32, T.star = 32, params, return.elog = TRUE)
cbs <- data$cbs
elog <- data$elog

#' estimate Pareto/NBD MCMC and Pareto/GGG MCMC
pggg.draws <- pggg.mcmc.DrawParameters(cbs, mcmc = 1500, burnin = 500, chains = 2, thin = 50)
pnbd.draws <- pnbd.mcmc.DrawParameters(cbs, mcmc = 1500, burnin = 500, chains = 2, thin = 50)
rbind(actual = params, 
      `Pareto/GGG` = round(summary(pggg.draws$level_2)$quantiles[, "50%"], 2),
      `Pareto/NBD` = c(t=NA, gamma=NA, round(summary(pnbd.draws$level_2)$quantiles[, "50%"], 2)))

#' plot estimated distribution of regularity rate `k`
pggg.mcmc.plotRegularityRateHeterogeneity(pggg.draws)

#' calculate effective sample size for MCMC chain
round(coda::effectiveSize(pggg.draws$level_2))

#' draw future transaction
pnbd.xstar <- mcmc.DrawFutureTransactions(cbs, pnbd.draws, T.star = cbs$T.star)
pggg.xstar <- mcmc.DrawFutureTransactions(cbs, pggg.draws, T.star = cbs$T.star)

#' calculate mean over future transaction draws for each customer
cbs$xstar.pnbd <- apply(pnbd.xstar, 2, mean)
cbs$xstar.pggg <- apply(pggg.xstar, 2, mean)

#' compare Mean Absolute Error
mae <- function(a, f) mean(abs(a-f))
cat('MAE Pareto/GGG:', mae(cbs$x.star, cbs$xstar.pggg), '\n') 
cat('MAE Pareto/NBD:', mae(cbs$x.star, cbs$xstar.pnbd), '\n') 
