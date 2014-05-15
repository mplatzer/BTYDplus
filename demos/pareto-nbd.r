# 
# ###   Step 0: simulate Pareto/NBD data   ###
# 
# n      <- 10000 # no. of customers
# T.cal  <- 9     # length of calibration period
# T.star <- 30    # length of hold-out period
# 
# # parameters of heterogeneity in intertransaction timings & churn probability
# params <- list(r=4.4, alpha=3.6, s=3, beta=35, k=1)
# 
# # generate data according to Pareto/NBD assumptions
# set.seed(1)
# cbs <- pcnbd.GenerateData(n, T.cal=T.cal, T.star=T.star, params=params)$cbs
# 
# mean(cbs$alive)
# # 0.4912
# 
# mean(cbs$x==0)
# # 0.0796
# 
# 
# ###   Step 1: Estimate Heterogeneity Parameters from Simulated Data   ###
# 
# params.est <- BTYD::pnbd.EstimateParameters(cbs, par.start=c(1,1,1,20))
# 
# # compare with actual parameters
# rbind("actual"=params, "Pareto/NBD"=round(params.est, 1))
# #            r   alpha s   beta
# # actual     4.4 3.6   3   35  
# # Pareto/NBD 4.2 3.5   3.4 40.3
# # --> we are able to re-estimate underlying heterogeneity parameters with maximum-likelihood-estimator
# 
# 
# ###   Step 2: Estimate Probability of a Customer still being Alive   ###
# 
# cbs$palive <- BTYD::pnbd.PAlive(params.est, cbs$x, cbs$t.x, cbs$T.cal)
# 
# mean((cbs$palive<.5)!=cbs$alive)
# # 0.898
# # --> 90% of customers are correctly classified
# 
# 
# ###   Step 3: Estimate No. of Transactions in Hold-Out Period   ###
# 
# cbs$est <- BTYD::pnbd.ConditionalExpectedTransactions(params.est, T.star, cbs$x, cbs$t.x, cbs$T.cal)
# 
# round(sum(cbs$x.star) / sum(cbs$est) - 1, 4)
# # 0.0359
# # --> Total No. of Transactions is missed by only 4%
# 
# cor(cbs$est, cbs$x.star)
# # 0.616
# # --> correlation between estimates and actuals
