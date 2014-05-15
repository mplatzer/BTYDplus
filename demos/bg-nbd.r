
###   Step 0: simulate BG/NBD data   ###

n      <- 5000 # no. of customers
T.cal  <- 9    # length of calibration period
T.star <- 30   # length of hold-out period

# parameters of heterogeneity in intertransaction timings & churn probability
params <- list(r=1.9, alpha=4.2, a=1.3, b=3.5)

# generate data according to BG/NBD assumptions
set.seed(1)
data <- bgnbd.GenerateData(n=n, T.cal=T.cal, T.star=T.star, params=params)
cbs <- data$cbs

sum(cbs$alive) / nrow(cbs)
# 0.391
# --> 61% of customers have already churned before the end of calibration period


###   Step 1: Estimate Heterogeneity Parameters from Simulated Data   ###

params.est <- bgnbd.EstimateParameters(cbs)

# compare with actual parameters
rbind("actual"=params, "BG/NBD"=round(params.est, 1))
#        r   alpha a   b  
# actual 1.9 4.2   1.3 3.5
# BG/NBD 2.1 4.6   1.4 3.9
# --> we are able to re-estimate underlying heterogeneity parameters with maximum-likelihood-estimator


###   Step 2: Estimate Probability of a Customer still being Alive   ###

cbs$palive <- bgnbd.PAlive(params.est, cbs$x, cbs$t.x, cbs$T.cal)

sum((cbs$palive<.5)!=cbs$alive)/nrow(cbs)
# 0.7154
# --> 72% of customers are correctly classified


###   Step 3: Estimate No. of Transactions in Hold-Out Period   ###

cbs$est <- bgnbd.ConditionalExpectedTransactions(params.est, T.star, cbs$x, cbs$t.x, cbs$T.cal)

round(sum(cbs$x.star) / sum(cbs$est) - 1, 4)
# 0.0101
# --> Total No. of Transactions is missed by only 1%

cor(cbs$est, cbs$x.star)
# 0.512
# --> correlation between estimates and actuals
