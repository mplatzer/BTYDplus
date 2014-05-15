
###   Step 0: simulate NBD data   ###

n      <- 5000 # no. of customers
T.cal  <- 9    # length of calibration period
T.star <- 30   # length of hold-out period

# parameters of heterogeneity in intertransaction timings
params <- list(r=1.9, alpha=4.2)

# generate data according to NBD assumptions
set.seed(1)
data <- nbd.GenerateData(n=n, T.cal=T.cal, T.star=T.star, params=params)
cbs <- data$cbs


###   Step 1: Estimate Heterogeneity Parameters from Simulated Data   ###

params.est <- nbd.EstimateParameters(cbs)

# compare with actual parameters
rbind("actual"=params, "NBD"=round(params.est, 1))
#        r   alpha
# actual 1.9 4.2  
# NBD    1.9 4.2 
# --> we are able to re-estimate underlying heterogeneity parameters with maximum-likelihood-estimator


###   Step 2: Estimate No. of Transactions in Hold-Out Period   ###

cbs$est <- nbd.ConditionalExpectedTransactions(params.est, T.star, cbs$x, cbs$T.cal)

round(sum(cbs$x.star) / sum(cbs$est) - 1, 4)
# 0.0061
# --> Total No. of Transactions is missed by less than 1%

cor(cbs$est, cbs$x.star)
# 0.7735612
# --> correlation between estimates and actuals
