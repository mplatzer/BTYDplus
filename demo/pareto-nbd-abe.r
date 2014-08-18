
# load CDnow data
data(cdnowSummary, package="BTYD")
data <- as.data.frame(cdnowSummary$cbs)

set.seed(1)
# estimate parameters
draws <- abe.mcmc.DrawParameters(data, burnin=5000, mcmc=5000, thin=100, chains=2)

agg <- function(x) round(c(quantile(x, 0.025), "mean"=mean(x), quantile(x, 0.975)), 2)
agg(as.matrix(draws$level_2)[, "log_lambda"])
agg(as.matrix(draws$level_2)[, "log_mu"])
agg((as.matrix(draws$level_2)[, "sigma_lambda_lambda"]))
agg((as.matrix(draws$level_2)[, "sigma_mu_mu"]))
agg((as.matrix(draws$level_2)[, "sigma_lambda_mu"]))
# -> parameter estimates match results shown in Abe's paper

xstar <- abe.mcmc.DrawFutureTransactions(data, draws, T.star=39)
mcmc.palive    <- abe.mcmc.PAlive(data, draws)
mcmc.est       <- colMeans(xstar)
#mcmc.est       <- apply(xstar, 2, median)
mcmc.est.pos   <- colMeans(xstar>0)

cor(data$x.star, mcmc.est)
# 0.62

mse <- function(x1, x2) sqrt(mean((x1-x2)^2))
mse(data$x.star, mcmc.est)
# 1.61

mae <- function(x1, x2) mean(abs(x1-x2))
mae(data$x.star, mcmc.est)
# 0.75
