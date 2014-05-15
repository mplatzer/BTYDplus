# 
# data(cdnowSummary)
# cbs <- as.data.frame(cdnowSummary$cbs)
# cbs$T.star <- 39
# 
# #################
# ### NBD Model ###
# #################
# 
# params <- nbd.EstimateParameters(cbs)
# round(params, 3)
# # 0.382 12.072
# 
# cbs$nbd <- nbd.ConditionalExpectedTransactions(params, cbs$T.star, cbs$x, cbs$T.cal)
# 
# 
# ########################
# ### Pareto-NBD Model ###
# ########################
# 
# # Note: these methods are part of BTYD-package
# params <- pnbd.EstimateParameters(cbs)
# round(params, 3)
# # 0.553 10.580  0.606 11.656
# 
# pnbd.cbs.LL(params, cbs)
# # -9594.976
# 
# p.alives <- pnbd.PAlive(params, cbs$x, cbs$t.x, cbs$T.cal)
# summary(p.alives)
# #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# #  0.0000  0.3174  0.3499  0.4462  0.4763  1.0000 
# plot(density(p.alives))
# 
# cbs$pnbd <- pnbd.ConditionalExpectedTransactions(params, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)
# 
# 
# ########################
# ### G/G-NBD Model ###
# ########################
# 
# params <- ggnbd.EstimateParameters(cbs, trace=10)
# round(params, 5)
# # 0.55259 10.56802  0.00013  0.60828  0.00148
# 
# ggnbd.cbs.LL(params, cbs)
# # -9594.976
# 
# cbs$ggnbd <- ggnbd.ConditionalExpectedTransactions(params, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal) 
# 
# 
# ####################
# ### BG/NBD model ###
# ####################
# 
# params <- bgnbd.EstimateParameters(cbs)
# round(params, 3)
# # 0.243 4.414 0.793 2.426
# 
# p.alives <- bgnbd.PAlive(params, cbs$x, cbs$t.x, cbs$T.cal)
# summary(p.alives)
# #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# #  0.0000  0.6604  1.0000  0.8134  1.0000  1.0000 
# plot(density(p.alives))
# 
# cbs$bgnbd <- bgnbd.ConditionalExpectedTransactions(params, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)
# 
# 
# #####################
# ### CBG/NBD model ###
# #####################
# 
# params <- cbgnbd.EstimateParameters(cbs)
# round(params, 3)
# # 0.525 6.183 0.891 1.614
# p.alives <- cbgnbd.PAlive(params, cbs$x, cbs$t.x, cbs$T.cal)
# plot(density(p.alives))
# cbs$cbgnbd <- cbgnbd.ConditionalExpectedTransactions(params, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)
# 
# 
# ########################
# ### CBG/CNBD-k model ###
# ########################
# 
# k <- 2
# params <- cbgcnbd.EstimateParameters(cbs, k)
# round(params, 3)
# # 0.580 2.288 1.135 1.585
# p.alives <- cbgcnbd.PAlive(params, k, cbs$x, cbs$t.x, cbs$T.cal)
# plot(density(p.alives))
# cbs$cbgcnbd2 <- cbgcnbd.ConditionalExpectedTransactions(params, k, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)
# 
# k <- 3
# params <- cbgcnbd.EstimateParameters(cbs, k)
# round(params, 3)
# # 0.533 1.258 1.289 1.845
# p.alives <- cbgcnbd.PAlive(params, k, cbs$x, cbs$t.x, cbs$T.cal)
# plot(density(p.alives))
# cbs$cbgcnbd3 <- cbgcnbd.ConditionalExpectedTransactions(params, k, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)
# 
# 
# library(plyr)
# head(ddply(cbs, .(x), summarize, 
#            PNBD=round(mean(pnbd), 2), GGNBD=round(mean(ggnbd), 2), 
#            BGNBD=round(mean(bgnbd), 2), CBGNBD=round(mean(cbgnbd), 2), 
#            CBGcNBD2=round(mean(cbgcnbd2), 2), CBGcNBD3=round(mean(cbgcnbd3), 2), 
#            ACT=round(mean(x.star), 2)))
# # x PNBD GGNBD BGNBD CBGNBD CBGcNBD2 CBGcNBD3  ACT
# # 1 0 0.14  0.14  0.23   0.19     0.15     0.15 0.24
# # 2 1 0.60  0.60  0.52   0.55     0.41     0.37 0.70
# # 3 2 1.20  1.20  1.04   1.05     0.84     0.77 1.39
# # 4 3 1.71  1.71  1.52   1.49     1.21     1.10 1.56
# # 5 4 2.40  2.40  2.16   2.10     1.79     1.66 2.53
# # 6 5 2.91  2.91  2.65   2.56     2.15     1.98 2.95
# 
# est <- data.frame(
#           type = rep(c("P/NBD", "G/G/NBD", "BG/NBD", "CBG/NBD", "CBG/CNBD-2", "CBG/CNBD-3"), each=nrow(cbs)),
#           fcst = c(cbs$pnbd, cbs$ggnbd, cbs$bgnbd, cbs$cbgnbd, cbs$cbgcnbd2, cbs$cbgcnbd3),
#           act = rep(cbs$x.star, 6))
# 
# MSLE <- function(a, f) { return(mean(((log(a+1) - log(f+1)))^2)) }
# MAPE <- function(a, f) { return(sum(abs(a-f)/sum(a))) }
# SUM <- function(a, f) { return((sum(f)-sum(a))/sum(a)) }
# ddply(est, .(type), summarize, 
#       msle=round(MSLE(act, fcst), 3),
#       mape=round(MAPE(act, fcst), 3),
#       sum=round(SUM(act, fcst), 3))
# # type  msle  mape    sum
# # 1     BG/NBD 0.234 0.984 -0.121
# # 2 CBG/CNBD-2 0.238 0.922 -0.310
# # 3 CBG/CNBD-3 0.245 0.916 -0.358
# # 4    CBG/NBD 0.231 0.958 -0.162
# # 5    G/G/NBD 0.238 0.945 -0.115
# # 6      P/NBD 0.238 0.945 -0.115
