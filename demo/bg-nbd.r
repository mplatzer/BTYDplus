
set.seed(1)

# generate artificial BG/NBD data 
n      <- 1000 # no. of customers
T.cal  <- 32   # length of calibration period
T.star <- 32   # length of hold-out period
params <- c(r=0.85, alpha=4.45, # purchase frequency lambda_i ~ Gamma(r, alpha)
               a=0.79, b=2.42)     # dropout probability p_i ~ Beta(a, b)

cbs <- bgnbd.GenerateData(n, T.cal, T.star, params)$cbs

# estimate parameters, and compare to true parameters
est <- bgnbd.EstimateParameters(cbs[, c("x", "t.x", "T.cal")])
rbind(params, est=round(est, 2))
#        r    alpha a    b   
# params 0.85 4.45  0.79 2.42
# est    0.82 4.32  0.62 1.98
# -> underlying parameters are successfully identified via Maximum Likelihood Estimation

# estimate future transactions in holdout-period
cbs$x.est <- bgnbd.ConditionalExpectedTransactions(est, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)

# compare forecast accuracy to naive forecast
c("bgnbd"=sqrt(mean((cbs$x.star-cbs$x.est)^2)),
  "naive"=sqrt(mean((cbs$x.star-cbs$x)^2)))
#    bgnbd    naive 
# 2.386392 3.688767
# -> BG/NBD forecast better than naive forecast

# estimate P(alive)
cbs$palive <- bgnbd.PAlive(est, cbs$x, cbs$t.x, cbs$T.cal)

# compare to true (usually unobserved) alive status
prop.table(table(cbs$palive>.5, cbs$alive))
#       FALSE  TRUE
# FALSE 0.359 0.081
# TRUE  0.098 0.462
# -> 82% of customers are correctly classified
