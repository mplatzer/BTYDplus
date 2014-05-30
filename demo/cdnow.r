
# Load CDNow event log
data(cdnowElog, package="BTYD")
elog <- data.frame(t(sapply(2:nrow(cdnowElog), function(i) strsplit(as.character(cdnowElog[i,]), split=",")[[1]])))
names(elog) <- c("cust", "sampleid", "date", "cds", "sales")
elog$date <- as.Date(elog$date, "%Y%m%d")

# Transform to CBS (including extra column 'litt')
cbs <- elog2cbs(elog, per="week", T.cal=as.Date("1997-09-30"))

x <- readline("Press Enter to continue...")

# Estimate models via Maximum Likelihood Estimation
(params.nbd <- nbd.EstimateParameters(cbs))
(params.pnbd <- pnbd.EstimateParameters(cbs))
(params.bgnbd <- bgnbd.EstimateParameters(cbs))
(params.cbgnbd <- cbgnbd.EstimateParameters(cbs))
(params.cbgcnbd <- cbgcnbd.EstimateParameters(cbs))
#(params.ggnbd <- ggnbd.EstimateParameters(cbs, trace=10)) # takes several minutes
(params.ggnbd <- c(r=0.552256, alpha=10.568448, b=0.000035, s=0.609494, beta=0.000417))

x <- readline("Press Enter to continue...")

# Compare Log-Likelihoods
rbind("NBD"=nbd.cbs.LL(params.nbd, cbs),
      "Pareto/NBD"=pnbd.cbs.LL(params.pnbd, cbs),
      "GG/NBD"=ggnbd.cbs.LL(params.ggnbd, cbs),
      "BG/NBD"=bgnbd.cbs.LL(params.bgnbd, cbs),
      "CBG/NBD"=cbgnbd.cbs.LL(params.cbgnbd, cbs),
      "CBG/CNBD-k"=cbgcnbd.cbs.LL(params.cbgcnbd, cbs))
# NBD        -9759.985
# Pareto/NBD -9591.412
# GG/NBD     -9591.412
# BG/NBD     -9578.556
# CBG/NBD    -9578.403; CBG/NBD provides best fit according to LL
# CBG/CNBD-k -9578.403; no regularity detected, hence k=1

x <- readline("Press Enter to continue...")

# Estimate Transactions in Holdout period
cbs$nbd <- nbd.ConditionalExpectedTransactions(params.nbd, cbs$T.star, cbs$x, cbs$T.cal)
cbs$pnbd <- pnbd.ConditionalExpectedTransactions(params.pnbd, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)
cbs$ggnbd <- ggnbd.ConditionalExpectedTransactions(params.ggnbd, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)
cbs$bgnbd <- bgnbd.ConditionalExpectedTransactions(params.bgnbd, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)
cbs$cbgnbd <- cbgnbd.ConditionalExpectedTransactions(params.cbgnbd, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)
cbs$cbgcnbd <- cbgcnbd.ConditionalExpectedTransactions(params.cbgcnbd, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)

x <- readline("Press Enter to continue...")

# Estimate P(alive)
cbs$palive.nbd <- 1
cbs$palive.pnbd <- pnbd.PAlive(params=params.pnbd, cbs$x, cbs$t.x, cbs$T.cal)
cbs$palive.ggnbd <- ggnbd.PAlive(params=params.ggnbd, cbs$x, cbs$t.x, cbs$T.cal)
cbs$palive.bgnbd <- bgnbd.PAlive(params=params.bgnbd, cbs$x, cbs$t.x, cbs$T.cal)
cbs$palive.cbgnbd <- cbgnbd.PAlive(params=params.cbgnbd, cbs$x, cbs$t.x, cbs$T.cal)
cbs$palive.cbgcnbd <- cbgcnbd.PAlive(params=params.cbgcnbd, cbs$x, cbs$t.x, cbs$T.cal)

x <- readline("Press Enter to continue...")

# compare forecasting accuracy
MAPE <- function(a, f) { return(sum(abs(a-f)/sum(a))) }
RMSE <- function(a, f) { return(sqrt(mean((a-f)^2))) }
MSLE <- function(a, f) { return(mean(((log(a+1) - log(f+1)))^2)) }

models <- c("nbd", "pnbd", "ggnbd", "bgnbd", "cbgnbd", "cbgcnbd")
acc <- t(sapply(models, function(model) c(MAPE(cbs$x.star, cbs[, model]),
                                          RMSE(cbs$x.star, cbs[, model]),
                                          MSLE(cbs$x.star, cbs[, model]))))
colnames(acc) <- c("MAPE", "RMSE", "MSLE")
round(acc, 3)
#          MAPE  RMSE  MSLE
# nbd     1.303 1.859 0.362
# pnbd    0.945 1.603 0.238
# ggnbd   0.945 1.603 0.238
# bgnbd   0.984 1.608 0.234
# cbgnbd  0.958 1.607 0.231
# cbgcnbd 0.958 1.607 0.231
