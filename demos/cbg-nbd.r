# 
# ###   Step 0: simulate CBG/NBD data   ###
# 
# n      <- 5000 # no. of customers
# T.cal  <- 9    # length of calibration period
# T.star <- 30   # length of hold-out period
# 
# # parameters of heterogeneity in intertransaction timings & churn probability
# params <- list(r=1.9, alpha=4.2, a=1.3, b=3.5)
# 
# # generate data according to CBG/NBD assumptions
# set.seed(1)
# data <- cbgnbd.GenerateData(n=n, T.cal=T.cal, T.star=T.star, params=params)
# cbs <- data$cbs
# 
# sum(cbs$alive) / nrow(cbs)
# # 0.3224
# # --> 68% of customers have already churned before the end of calibration period
# 
# 
# ###   Step 1: Estimate Heterogeneity Parameters from Simulated Data   ###
# 
# params.cbg <- cbgnbd.EstimateParameters(cbs)
# params.bg <- bgnbd.EstimateParameters(cbs) # for comparison
# 
# # compare with actual parameters
# rbind("actual"=params, "CBG/NBD"=round(params.cbg, 1), "BG/NBD"=round(params.bg, 1))
# #         r   alpha a   b   
# # actual  1.9 4.2   1.3 3.5 
# # CBG/NBD 1.6 3.4   1.4 4.2 
# # BG/NBD  0.6 1.7   2.4 11.2
# # --> we are able to re-estimate underlying heterogeneity parameters with maximum-likelihood-estimator
# 
# 
# ###   Step 2: Estimate Probability of a Customer still being Alive   ###
# 
# cbs$palive.cbg <- cbgnbd.PAlive(params.cbg, cbs$x, cbs$t.x, cbs$T.cal)
# cbs$palive.bg <- bgnbd.PAlive(params.bg, cbs$x, cbs$t.x, cbs$T.cal)
# 
# # calculate accuracy 
# c("cbg/nbd"=sum((cbs$palive.cbg>=.5)==cbs$alive), 
#   "bg/nbd"=sum((cbs$palive.bg>=.5)==cbs$alive)) / nrow(cbs)
# # cbg/nbd  bg/nbd 
# #  0.7368  0.4902 
# 
# # --> accuracy of detected churn could be lifted from 49% to 74% by allowing drop-out opportunity at t=0
# 
# 
# ###   Step 3: Estimate No. of Transactions in Hold-Out Period   ###
# 
# cbs$est <- cbgnbd.ConditionalExpectedTransactions(params.cbg, T.star, cbs$x, cbs$t.x, cbs$T.cal)
# 
# round(sum(cbs$x.star) / sum(cbs$est) - 1, 4)
# # -0.0441
# # --> Total No. of Transactions is missed by only 4%
# 
# cor(cbs$est, cbs$x.star)
# # 0.487
# # --> correlation between estimates and actuals
