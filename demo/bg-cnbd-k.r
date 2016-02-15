
set.seed(1)

# generate artificial BG/CNBD-k data 
n      <- 8000 # no. of customers
T.cal  <- 32   # length of calibration period
T.star <- 32   # length of hold-out period
params <- c(k=3,                # regularity in interpurchase-times (Erlang-k)
            r=0.85, alpha=1.45, # purchase frequency lambda_i ~ Gamma(r, alpha)
            a=0.79, b=2.42)     # dropout probability p_i ~ Beta(a, b)

data <- bgcnbd.GenerateData(n, T.cal, T.star, params, return.elog=TRUE)
cbs  <- data$cbs  # CBS summary - one record per customer
elog <- data$elog # Event log - one row per event/purchase

# estimate regularity from event log
(k.est <- estimateRegularity(elog))
# 2.927712; interpurchase-times indicate Erlang-3

# estimate parameters, and compare to true parameters
est  <- bgcnbd.EstimateParameters(cbs[, c("x", "t.x", "T.cal", "litt")])
est1 <- bgnbd.EstimateParameters(cbs[, c("x", "t.x", "T.cal", "litt")])
rbind("actual"=params, "bg/cnbd-k"=round(est, 2), "bg/nbd"=c(1, round(est1, 2)))
#           k    r alpha    a    b
# actual    3 0.85  1.45 0.79 2.42
# bg/cnbd-k 3 0.85  1.48 0.77 2.35
# bg/nbd    1 0.91  6.02 0.59 2.21
# -> underlying parameters are successfully identified via Maximum Likelihood Estimation

# estimate future transactions in holdout-period
cbs$x.est  <- bgcnbd.ConditionalExpectedTransactions(est, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)
cbs$x.est1 <- bgnbd.ConditionalExpectedTransactions(est1, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)

# compare forecast accuracy to bg/nbd and naive forecast
rbind("bg/cnbd-k" = mean(abs(cbs$x.star-cbs$x.est)),
      "bg/nbd"    = mean(abs(cbs$x.star-cbs$x.est1)),
      "naive"     = mean(abs(cbs$x.star-cbs$x)))
# bg/cnbd-k 1.056313
# bg/nbd    1.212169
# naive     2.018125
# -> BG/CNBD-k forecast better than BG/NBD and naive forecast

# estimate P(alive)
cbs$palive  <- bgcnbd.PAlive(est, cbs$x, cbs$t.x, cbs$T.cal)
cbs$palive1 <- bgnbd.PAlive(est1, cbs$x, cbs$t.x, cbs$T.cal)

# compare to true (usually unobserved) alive status
prop.table(table(cbs$palive>.5, cbs$alive))
#         FALSE    TRUE
#   FALSE 0.36800 0.04475
#   TRUE  0.08850 0.49875
# -> 87% of customers are correctly classified

# Brier score for P(alive)
rbind("bg/cnbd-k" = sqrt(mean((cbs$palive-cbs$alive)^2)),
      "bg/nbd"    = sqrt(mean((cbs$palive1-cbs$alive)^2)))
# bg/cnbd-k 0.3069095
# bg/nbd    0.3310709
# -> P(alive) is more accurate for BG/CNBD-k than for BG/NBD when regularity
# is present in the data


# compare estimated with actual distributions in lambda & churn probability
par(mfrow=c(2,1), mar=c(2,1,2,1))
xlim <- 1.5
x <- seq(0, xlim, len=1000)[-1]
y <- dgamma(x, shape=params[2], rate=params[3]*params[1])
plot(x, y, typ="l", col="black", lwd=2, main="Heterogeneity in Intertransaction Times", ylab="", xlab="", axes=FALSE, xlim=c(0, xlim))
lines(x, dgamma(x, shape=est[2], rate=est[3]*est[1]), col="red", lwd=2)
lines(x, dgamma(x, shape=est1[1], rate=est1[2]), col="blue", lwd=2)
axis(1, pos=0, labels=round(1/(xlim*(0:10/10)), 1), at=xlim*(0:10/10))
legend(xlim*.6, max(y)*.9, c("Actual", "BG/CNBD-k", "BG/NBD"), col=c("black", "red", "blue"), pch=15, bty="n")

xlim <- 1
x <- seq(0.05, xlim, len=1000)
y <- dbeta(x, params[4], params[5])
plot(x, y, typ="l", col="black", lwd=2, main="Heterogeneity in Churn Probability", ylab="", xlab="", axes=FALSE, xlim=c(0, xlim))
lines(x, dbeta(x, est[4], est[5]), col="red", lwd=2)
lines(x, dbeta(x, est1[3], est1[4]), col="blue", lwd=2)
axis(1, pos=0)
legend(xlim*.6, max(y)*.9, c("Actual", "BG/CNBD-k", "BG/NBD"), col=c("black", "red", "blue"), pch=15, bty="n")
