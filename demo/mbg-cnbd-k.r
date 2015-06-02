
set.seed(1)

# generate artificial MBG/CNBD-k data 
n      <- 1000 # no. of customers
T.cal  <- 32   # length of calibration period
T.star <- 32   # length of hold-out period
params <- c(k=3,                # regularity in interpurchase-times (Erlang-k)
            r=0.85, alpha=1.45, # purchase frequency lambda_i ~ Gamma(r, alpha)
            a=0.79, b=2.42)     # dropout probability p_i ~ Beta(a, b)

data <- mbgcnbd.GenerateData(n, T.cal, T.star, params, return.elog=TRUE)
cbs  <- data$cbs  # CBS summary - one record per customer
elog <- data$elog # Event log - one row per event/purchase

# estimate regularity from event log
(k.est <- estimateRegularity(elog))
# 2.784212; interpruchase-times indicate Erlang-3

# estimate parameters, and compare to true parameters
est <- mbgcnbd.EstimateParameters(cbs[, c("x", "t.x", "T.cal", "litt")])
est1 <- mbgnbd.EstimateParameters(cbs[, c("x", "t.x", "T.cal", "litt")])
rbind("actual"=params, "mbg/cnbd-k"=round(est, 2), "mbg/nbd"=c(1, round(est1, 2)))
#            k    r alpha    a    b
# actual     3 0.85  1.45 0.79 2.42
# mbg/cnbd-k 3 0.76  1.35 0.81 2.68
# mbg/nbd    1 0.56  4.22 0.86 5.51
# -> underlying parameters are successfully identified via Maximum Likelihood Estimation

# estimate future transactions in holdout-period
cbs$x.est <- mbgcnbd.ConditionalExpectedTransactions(est, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)
cbs$x.est1 <- mbgnbd.ConditionalExpectedTransactions(est1, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)

# compare forecast accuracy to mbg/nbd and naive forecast
rbind("mbg/cnbd-k" = sqrt(mean((cbs$x.star-cbs$x.est)^2)),
      "mbg/nbd" = sqrt(mean((cbs$x.star-cbs$x.est1)^2)),
      "naive"   = sqrt(mean((cbs$x.star-cbs$x)^2)))
# mbg/cnbd-k 1.814616
# mbg/nbd    1.870667
# naive      3.098225
# -> MBG/CNBD-k forecast better than MBG/NBD and naive forecast

# estimate P(alive)
cbs$palive <- mbgcnbd.PAlive(est, cbs$x, cbs$t.x, cbs$T.cal)
cbs$palive1 <- mbgnbd.PAlive(est1, cbs$x, cbs$t.x, cbs$T.cal)

# compare to true (usually unobserved) alive status
prop.table(table(cbs$palive>.5, cbs$alive))
#       FALSE  TRUE
# FALSE 0.479 0.165
# TRUE  0.057 0.299
# -> 78% of customers are correctly classified

# Brier score for P(alive)
rbind("mbg/cnbd-k" = sqrt(mean((cbs$palive-cbs$alive)^2)),
      "mbg/nbd" = sqrt(mean((cbs$palive1-cbs$alive)^2)))
# mbg/cnbd-k 0.3849790
# mbg/nbd    0.4528768
# -> P(alive) is more accurate for MBG/CNBD-k than for MBG/NBD when regularity
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
legend(xlim*.6, max(y)*.9, c("Actual", "MBG/CNBD-k", "MBG/NBD"), col=c("black", "red", "blue"), pch=15, bty="n")

xlim <- 1
x <- seq(0.05, xlim, len=1000)
y <- dbeta(x, params[4], params[5])
plot(x, y, typ="l", col="black", lwd=2, main="Heterogeneity in Churn Probability", ylab="", xlab="", axes=FALSE, xlim=c(0, xlim))
lines(x, dbeta(x, est[4], est[5]), col="red", lwd=2)
lines(x, dbeta(x, est1[3], est1[4]), col="blue", lwd=2)
axis(1, pos=0)
legend(xlim*.6, max(y)*.9, c("Actual", "MBG/CNBD-k", "MBG/NBD"), col=c("black", "red", "blue"), pch=15, bty="n")
