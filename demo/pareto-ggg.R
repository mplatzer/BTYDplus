
set.seed(1)

# generate artificial Pareto/GGG data
n <- 1000  # no. of customers
T.cal <- 32  # length of calibration period
T.star <- 32  # length of hold-out period
params <- list(t = 4.5, gamma = 1.5, r = 5, alpha = 10, s = 0.8, beta = 12)
# regularity parameter k_i ~ Gamma(t, gamma) purchase frequency lambda_i ~ Gamma(r, alpha) dropout probability
# mu_i ~ Gamma(s, beta)

data <- pggg.GenerateData(n, T.cal, T.star, params, return.elog = TRUE)
cbs <- data$cbs
elog <- data$elog

# estimate Pareto/NBD MCMC
pnbd.draws <- pnbd.mcmc.DrawParameters(cbs, mcmc = 1500, burnin = 500, chains = 2, thin = 50)
plot(pnbd.draws$level_2)

# estimate Pareto/GGG MCMC
pggg.draws <- pggg.mcmc.DrawParameters(cbs, mcmc = 1500, burnin = 500, chains = 2, thin = 50)
plot(pggg.draws$level_2, density = FALSE)
rbind(actual = params, estimated = round(summary(pggg.draws$level_2)$quantiles[, "50%"], 2))
# t gamma r alpha s beta actual 4.5 1.5 5 10 0.8 12 estimated 3.7 1.2 4.85 9.72 0.71 9.2

coda::gelman.diag(pggg.draws$level_2)
# -> MCMC chains have not converged yet

pggg.mcmc.plotRegularityRateHeterogeneity(pggg.draws)

round(coda::effectiveSize(pggg.draws$level_2))
# -> effective sample size are small for such a short chain

# draw future transaction
pnbd.xstar <- mcmc.DrawFutureTransactions(cbs, pnbd.draws, T.star = cbs$T.star)
pggg.xstar <- mcmc.DrawFutureTransactions(cbs, pggg.draws, T.star = cbs$T.star)

# calculate mean over future transaction draws for each customer
cbs$pnbd.mcmc <- apply(pnbd.xstar, 2, mean)
cbs$pggg.mcmc <- apply(pggg.xstar, 2, mean)

MAPE <- function(a, f) {
  return(sum(abs(a - f)/sum(a)))
}
RMSE <- function(a, f) {
  return(sqrt(mean((a - f)^2)))
}
MSLE <- function(a, f) {
  return(mean(((log(a + 1) - log(f + 1)))^2))
}
BIAS <- function(a, f) {
  return(mean(a) - mean(f))
}
bench <- function(cbs, models) {
  acc <- t(sapply(models, function(model) c(MAPE(cbs$x.star, cbs[, model]), RMSE(cbs$x.star, cbs[, model]), MSLE(cbs$x.star, 
    cbs[, model]), BIAS(cbs$x.star, cbs[, model]))))
  colnames(acc) <- c("MAPE", "RMSE", "MSLE", "BIAS")
  round(acc, 3)
}

bench(cbs, c("pnbd.mcmc", "pggg.mcmc"))
# MAPE RMSE MSLE BIAS pnbd.mcmc 0.532 4.485 0.477 -0.21 pggg.mcmc 0.488 4.437 0.379 -0.02

# calculate P(active)
cbs$pactive.pnbd.mcmc <- apply(pnbd.xstar, 2, function(x) mean(x > 0))
cbs$pactive.pggg.mcmc <- apply(pggg.xstar, 2, function(x) mean(x > 0))

# Brier score
c(pnbd.mcmc = mean((cbs$pactive.pnbd.mcmc - (cbs$x.star > 0))^2), pggg.mcmc = mean((cbs$pactive.pggg.mcmc - (cbs$x.star > 
  0))^2))
# pnbd.mcmc pggg.mcmc 0.04344944 0.03808333

# calculate P(alive)
cbs$palive.pnbd.mcmc <- mcmc.PAlive(cbs, pnbd.draws)
cbs$palive.pggg.mcmc <- mcmc.PAlive(cbs, pggg.draws) 
