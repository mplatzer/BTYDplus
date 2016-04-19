
set.seed(1)

# generate artificial GG/NBD data 
n      <- 5000 # no. of customers
T.cal  <- 32   # length of calibration period
T.star <- 32   # length of hold-out period
params <- c(r=0.85, alpha=4.45,        # purchase frequency lambda_i ~ Gamma(r, alpha)
            b=0.12, s=0.75, beta=9.55) # lifetime tau_i ~ Gompertz(b, eta_i); eta_i ~ Gamma(s, beta)

cbs <- ggnbd.GenerateData(n, T.cal, T.star, params)$cbs

# estimate parameters, and compare to true parameters
est <- ggnbd.EstimateParameters(cbs[, c("x", "t.x", "T.cal")], trace=10)
rbind(params, est=round(est, 2))
#           r alpha    b    s beta
# params 0.85  4.45 0.12 0.75 9.55
# est    0.86  4.45 0.13 0.57 8.19
# -> underlying parameters are successfully identified via Maximum Likelihood Estimation

# estimate future transactions in holdout-period
cbs$x.est <- ggnbd.ConditionalExpectedTransactions(est, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)

# compare forecast accuracy to naive forecast
c("ggnbd"=sqrt(mean((cbs$x.star-cbs$x.est)^2)),
  "naive"=sqrt(mean((cbs$x.star-cbs$x)^2)))
#    ggnbd    naive 
# 2.871434 6.104031 
# -> Gamma/Gompertz/NBD forecast better than naive forecast

# estimate P(alive)
cbs$palive <- ggnbd.PAlive(est, cbs$x, cbs$t.x, cbs$T.cal)

# compare to true (usually unobserved) alive status
prop.table(table(cbs$palive>.5, cbs$alive))
#        FALSE   TRUE
# FALSE 0.6878 0.1264
# TRUE  0.0420 0.1438
# -> 83% of customers are correctly classified
