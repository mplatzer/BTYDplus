#' Load transaction records of 1525 grocery customers.
data("groceryElog", envir = environment())
head(groceryElog)
range(groceryElog$date)

#' Convert from event log to customer-by-sufficient-statistic summary.
#' Split into 52 weeks calibration, and 52 weeks holdout period.
cbs <- elog2cbs(groceryElog, T.cal = "2006-12-31", T.tot = "2007-12-30")
head(cbs)


x <- readline("Estimate Pareto/NBD (HB) and Pareto/GGG (press Enter)")

#' Draw Pareto/NBD parameters and compare median estimates with Pareto/NBD
#' implementation from BTYD package.
pnbd.draws <- pnbd.mcmc.DrawParameters(cal.cbs = cbs, mc.cores = 1)

#' Draw Pareto/GGG parameters and report median estimates. Note, that we only
#' use a relatively short MCMC chain to keep this demo finish reasonably fast.
pggg.draws <- pggg.mcmc.DrawParameters(cal.cbs = cbs,
                 mcmc = 500, burnin = 500, chains = 2, thin = 20,
                 mc.cores = 1)

round(rbind(`Pareto/GGG`= summary(pggg.draws$level_2)$quantiles[, "50%"],
      `Pareto/NBD (HB)` = c(NA, NA, summary(pnbd.draws$level_2)$quantiles[, "50%"]),
      `Pareto/NBD`      = c(NA, NA, BTYD::pnbd.EstimateParameters(cbs))), 2)

#' plot estimated parameter distributions
plot(pggg.draws$level_2, density = TRUE, trace = FALSE)
pggg.plotRegularityRateHeterogeneity(pggg.draws)
# -> regularity detected in grocery dataset

#' check MCMC convergence
plot(pggg.draws$level_2, density = FALSE, trace = TRUE)
coda::effectiveSize(pggg.draws$level_2)


x <- readline("Estimate future transactions (press Enter)")

#' draw future transaction
xstar.pnbd.draws <- mcmc.DrawFutureTransactions(cbs, pnbd.draws, T.star = cbs$T.star, sample_size = 400)
xstar.pggg.draws <- mcmc.DrawFutureTransactions(cbs, pggg.draws, T.star = cbs$T.star, sample_size = 400)

#' calculate mean over future transaction draws for each customer
cbs$xstar.pnbd <- apply(xstar.pnbd.draws, 2, mean)
cbs$xstar.pggg <- apply(xstar.pggg.draws, 2, mean)

#' calculate P(alive)
cbs$palive.pnbd <- mcmc.PAlive(pnbd.draws)
cbs$palive.pggg <- mcmc.PAlive(pggg.draws)

#' calculate P(active)
cbs$pactive.pnbd <- mcmc.PActive(xstar.pnbd.draws)
cbs$pactive.pggg <- mcmc.PActive(xstar.pggg.draws)

#' compare forecast accuracy to Pareto/NBD
(mae <- c(`Pareto/GGG` = mean(abs(cbs$x.star - cbs$xstar.pggg)),
          `Pareto/NBD` = mean(abs(cbs$x.star - cbs$xstar.pnbd))))
(lift <- 1 - mae[1]/mae[2])
#' -> 9% lift in customer-level accuracy when taking regularity into account

# P(active) diagnostic plot
nil <- mcmc.plotPActiveDiagnostic(cbs, xstar.pggg.draws)
