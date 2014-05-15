# 
# ###   Step 0: simulate CBG/CNBD-k data   ###
# 
# n      <- 5000  # no. of customers
# T.cal  <- 9     # length of calibration period
# T.star <- 30    # length of hold-out period
# k      <- 3     # regularity in interpurchase-times (Erlang-k)
# 
# # parameters of heterogeneity in intertransaction timings & churn probability
# params <- list(r=1.9, alpha=1.4, a=1.3, b=3.5)
# 
# # generate data according to CBG/CNBD-k assumptions
# set.seed(1)
# data <- cbgcnbd.GenerateData(n=n, k=k, T.cal=T.cal, T.star=T.star, params=params, return.elog=T)
# 
# cbs <- data$cbs
# elog <- data$elog
# 
# sum(cbs$alive) / nrow(cbs)
# # [1] 0.3328
# # --> 67% of customers have already churned before the end of calibration period
# 
# table(cbs$x)
# #    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   18 
# # 1699 1164  803  466  319  203  120   86   48   30   20   22    6   10    1    2    1 
# # --> about a third of customers never made a repurchase, whereas
# #       some customers had up to 18 transactions recorded 
# 
# 
# ###   Step 1: Estimate Regularity Parameter from Simulated Data   ###
# 
# # Estimate k by Wheat&Morrison method
# elog.cal <- subset(elog, t<=T.cal)
# (k.est <- cbgcnbd.EstimateRegularity(elog.cal))
# # [1] 3.02553
# k.est <- round(k.est)
# # k.est = 3
# # --> we are able to re-estimate correct regularity parameter with Wheat&Morrison method
# 
# 
# ###   Step 2: Estimate Heterogeneity Parameters from Simulated Data   ###
# 
# params.k <- cbgcnbd.EstimateParameters(cbs, k=k.est)
# # for comparison we also compute numbers based on CBD/NBD model, and other values for k
# params.1 <- cbgcnbd.EstimateParameters(cbs, k=1) 
# params.2 <- cbgcnbd.EstimateParameters(cbs, k=2) 
# params.5 <- cbgcnbd.EstimateParameters(cbs, k=5) 
# rbind("actual"=params, "CBG/CNBD-k"=round(params.k, 1), "CBG/NBD"=round(params.1, 1), "CBG/CNBD-3"=round(params.2, 1), "CBG/CNBD-5"=round(params.5, 1))
# #            r   alpha a   b  
# # actual     1.9 1.4   1.3 3.5
# # CBG/CNBD-k 2   1.4   1.2 3.1
# # CBG/NBD    1.8 6     0.4 1.7
# # CBG/CNBD-3 2   2.4   1   2.6
# # CBG/CNBD-5 1.8 0.7   1.5 3.8
# # --> we are able to re-estimate underlying heterogeneity parameters with maximum-likelihood-estimator
# 
# # compare estimated with actual distributions in lambda & churn probability
# par(mfrow=c(2,1), mar=c(2,1,2,1))
# 
# xlim <- 1.5
# ylim <- 2.4
# x <- seq(0, xlim, len=1000)
# y <- dgamma(x, shape=params$r, rate=params$alpha*k)
# plot(x, y, typ="l", col="black", lwd=2, main="Heterogeneity in Intertransaction Times", ylab="", xlab="", axes=F, xlim=c(0, xlim), ylim=c(0, ylim))
# #polygon(c(0, x, 5), c(0, y, 0), col="gray")
# lines(x, dgamma(x, shape=params.k[1], rate=params.k[2]*k.est), col="red", lwd=2)
# lines(x, dgamma(x, shape=params.1[1], rate=params.1[2]), col="blue", lwd=2)
# axis(1, pos=0, labels=round(1/(xlim*(0:10/10)), 1), at=xlim*(0:10/10))
# #axis(3, pos=ylim, labels=round(1/(xlim*(0:10/10)), 1), at=xlim*(0:10/10), cex=.5)
# segments(params$r/params$alpha/k, 0, params$r/params$alpha/k, ylim*0.1, lty=1, lwd=2, col="black")
# segments(params.k[1]/params.k[2]/k.est, 0, params.k[1]/params.k[2]/k.est, ylim*0.1,lty=1, lwd=2, col="red")
# segments(params.1[1]/params.1[2]/1, 0, params.1[1]/params.1[2]/1, ylim*0.1, lty=1, lwd=2, col="blue")
# legend(xlim*.6, ylim*.9, c("Actual", "CBG/CNBD-k", "CBG/NBD"), col=c("black", "red", "blue"), pch=15, bty="n")
# 
# xlim <- 1
# ylim <- 4
# x <- seq(0, xlim, len=1000)
# y <- dbeta(x, params$a, params$b)
# plot(x, y, typ="l", col="black", lwd=2, main="Heterogeneity in Churn Probability", ylab="", xlab="", axes=F, xlim=c(0, xlim), ylim=c(0, ylim))
# lines(x, dbeta(x, params.k[3], params.k[4]), col="red", lwd=2)
# lines(x, dbeta(x, params.1[3], params.1[4]), col="blue", lwd=2)
# axis(1, pos=0)
# segments(params$a/(params$a+params$b), 0, params$a/(params$a+params$b), ylim*0.1, lty=1, lwd=2, col="black")
# segments(params.k[3]/(params.k[3]+params.k[4]), 0, params.k[3]/(params.k[3]+params.k[4]), ylim*0.1,lty=1, lwd=2, col="red")
# segments(params.1[3]/(params.1[3]+params.1[4]), 0, params.1[3]/(params.1[3]+params.1[4]), ylim*0.1, lty=1, lwd=2, col="blue")
# legend(xlim*.6, ylim*.9, c("Actual", "CBG/CNBD-k", "CBG/NBD"), col=c("black", "red", "blue"), pch=15, bty="n")
# 
# 
# ###   Step 3: Estimate Probability of a Customer still being Alive   ###
# 
# cbs$palive.k <- cbgcnbd.PAlive(params.k, k=k.est, cbs$x, cbs$t.x, cbs$T.cal)
# cbs$palive.1 <- cbgcnbd.PAlive(params.1, k=1,     cbs$x, cbs$t.x, cbs$T.cal)
# 
# # Share of 'Still Alive'
# c("actual"=sum(cbs$alive), "cbg/cnbd-k"=sum(cbs$palive.k>=.5), "cbg/nbd"=sum(cbs$palive.1>=.5)) / nrow(cbs)
# #     actual cbg/cnbd-k    cbg/nbd 
# #     0.3328     0.3732     0.5470 
# 
# # calculate accuracy
# c("cbg/cnbd-k"=sum((cbs$palive.k>=.5)==cbs$alive), 
#   "cbg/nbd"=sum((cbs$palive.1>=.5)==cbs$alive)) / nrow(cbs)
# # cbg/cnbd-k    cbg/nbd 
# #     0.7688     0.6642 
# 
# # --> accuracy of detected churn could be lifted from 66% to 77% by taking regularity into account
# 
# 
# ###   Step 4: Estimate No. of Transactions in Hold-Out Period   ###
# 
# cbs$est.k <- cbgcnbd.ConditionalExpectedTransactions(params.k, k=k.est, T.star, cbs$x, cbs$t.x, cbs$T.cal)
# cbs$est.1 <- cbgcnbd.ConditionalExpectedTransactions(params.1, k=1,     T.star, cbs$x, cbs$t.x, cbs$T.cal)
# 
# # sum
# round(c("actual"     = sum(cbs$x.star), 
#         "cbg/cnbd-k" = sum(cbs$est.k), 
#         "cbg/nbd"    = sum(cbs$est.1)),0)
# #     actual cbg/cnbd-k    cbg/nbd 
# #      11643      11415      19935 
# 
# # --> Total No. of Transactions is missed by only 2%
# #     Not taking regularity into account would overestimate by 71%, 
# #        as many customers were still assumed to be active at end of calibration period
# 
# # root mean squared error
# rmse <- function(est, act) { return(sqrt(mean((est-act)^2))) }
# round(c("cbg/cnbd-k" = rmse(cbs$est.k, cbs$x.star), 
#         "cbg/nbd"    = rmse(cbs$est.1, cbs$x.star)), 3)
# # cbg/cnbd-k    cbg/nbd 
# #      3.800      4.346 
# 
# # mean squared logarithmic error
# msle <- function(est, act) { return(mean((log(est+1)-log(act+1))^2)) }
# round(c("cbg/cnbd-k" = msle(cbs$est.k, cbs$x.star), 
#         "cbg/nbd"    = msle(cbs$est.1, cbs$x.star)), 3)
# # cbg/cnbd-k    cbg/nbd 
# #      0.612      1.199
# 
# # --> Forecast error on individual level improves if regularity is taking into account
