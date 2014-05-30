
# TODO

# 
# ## CDNow ##
# 
# data(cdnowSummary, package="BTYD")
# data <- as.data.frame(cdnowSummary$cbs)
# 
# # Maximum Likelihood Estimate
# params <- BTYD::pnbd.EstimateParameters(data, c(10,10,10,10))
# 
# # MCMC Estimate with augmented parameter space; 2.5secs for 1000 users x 1000 steps
# draws1 <- pnbd.mcmc.DrawParameters(data, mcmc=10000, burnin=10000, thin=10, chains=2)
# plot(draws1$level_2)
# 
# # MCMC Estimate without augmented parameter space; 25secs for 1000 users x 1000 steps
# draws2 <- pnbd.mcmc.DrawParameters(data, mcmc=10000, burnin=10000, thin=10, chains=2, use_data_augmentation=F)
# 
# round(rbind(params, 
#             draws1=summary(draws1$level_2, quantiles=0.5)$quantiles,
#             draws2=summary(draws2$level_2, quantiles=0.5)$quantiles), 2)
# #           r alpha    s  beta
# # params 0.55 10.58 0.61 11.65
# # draws1 0.55 10.50 0.63 12.96
# # draws2 0.56 10.61 0.51  8.88
# 
# 
# ## Simulated Data ##
# 
# N      <- 5000  # no. of customers
# T.cal  <- 10    # length of calibration period
# T.star <- 5     # length of hold-out period
# params <- list(r=1.4, alpha=1.3, s=0.7, beta=7)
# 
# # generate data according to Pareto/NBD assumptions
# set.seed(1)
# data <- pcnbd.GenerateData(N, T.cal=T.cal, T.star=T.star, params=params)$cbs
# 
# # start MCMC and draw convergence
# set.seed(1)
# draws <- pnbd.mcmc.DrawParameters(data, mcmc=10000, burnin=10000, thin=10)
# 
# # Predict future transactions
# mcmc.est.draws <- pcnbd.mcmc.DrawFutureTransactions(data, draws, T.star)
# mcmc.palive    <- pcnbd.mcmc.PAlive(data, draws)
# mcmc.est       <- colMeans(mcmc.est.draws)
# mcmc.est.pos   <- colMeans(mcmc.est.draws>0)
# 
# table(data$x.star>0, mcmc.est.pos>.5)
# #       FALSE TRUE
# # FALSE  2268  412
# # TRUE    305 2015
# 
# # Benchmark with Maximum Likelihood Estimates
# params.est <- BTYD::pnbd.EstimateParameters(data)
# mle.est    <- BTYD::pnbd.ConditionalExpectedTransactions(params.est, T.star, data$x, data$t.x, data$T.cal)
# mle.palive <- BTYD::pnbd.PAlive(params.est, data$x, data$t.x, data$T.cal)
# 
# mae <- function(est, act) mean(abs(est-act))
# c(mle=mae(mle.est, data$x.star), mcmc=mae(mcmc.est, data$x.star))
# #      mle     mcmc 
# # 1.476639 1.477556 
# # -> Similar performance in terms of absolute number of future transactions
# 
# c(mle=mean((mle.palive>.5)==data$alive), mcmc=mean((mcmc.palive>.5)==data$alive))
# #    mle   mcmc 
# # 0.8746 0.8752 
# # -> Similar performance in terms of P(alive), i.e. detecting churn
# 
# 
# ## Analyse Convergence ##
# 
# effectiveSize(draws$level_2)
# plot(draws$level_2)
# plot(draws$level_1[[1]])
