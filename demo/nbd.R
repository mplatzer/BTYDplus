
set.seed(1)

# generate artificial NBD data
n <- 1000  # no. of customers
T.cal <- 32  # length of calibration period
T.star <- 32  # length of hold-out period
params <- c(r = 0.85, alpha = 4.45)  # purchase frequency lambda_i ~ Gamma(r, alpha)

cbs <- nbd.GenerateData(n, T.cal, T.star, params)$cbs

# estimate parameters, and compare to true parameters
est <- nbd.EstimateParameters(cbs[, c("x", "T.cal")])
rbind(params, est = round(est, 2))
# r alpha params 0.85 4.45 est 0.84 4.56 -> underlying parameters are successfully identified via Maximum
# Likelihood Estimation

# estimate future transactions in holdout-period
cbs$x.est <- nbd.ConditionalExpectedTransactions(est, cbs$T.star, cbs$x, cbs$T.cal)

# compare forecast accuracy to naive forecast
c(nbd = sqrt(mean((cbs$x.star - cbs$x.est)^2)), naive = sqrt(mean((cbs$x.star - cbs$x)^2)))
# nbd naive 3.446776 3.469582 -> NBD forecast only marginally better than naive forecast 
