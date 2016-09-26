
#' Load transaction records of 2357 CDNow customers.
cdnowElog <- read.csv(system.file("data/cdnowElog.csv", package = "BTYD"),
                      stringsAsFactors = FALSE,
                      col.names = c("cust", "sampleid", "date", "cds", "sales"))
cdnowElog$date <- as.Date(as.character(cdnowElog$date), format = "%Y%m%d")
head(cdnowElog)
range(cdnowElog$date)

#' Convert from event log to customer-by-sufficient-statistic summary.
#' Split into 39 weeks calibration, and 39 weeks holdout period.
cbs <- elog2cbs(cdnowElog, T.cal = "1997-09-30", T.tot = "1998-06-30")


x <- readline("Estimate Pareto/NBD (Abe) without covariates (press Enter)")

#' Estimate with no covariates; see model M1 in Abe (2009)
draws.m1 <- abe.mcmc.DrawParameters(
              cal.cbs = cbs,
              mcmc = 5000, burnin = 5000,
              mc.cores = 1)
round(summary(draws.m1$level_2)$quantiles[, c("2.5%", "50%", "97.5%")], 2)
#' -> Parameter Estimates match Table 3 in Abe (2009);


x <- readline("Estimate Pareto/NBD (Abe) with covariates (press Enter)")

#' Append dollar amount of first purchase to use as covariate
first <- aggregate(sales ~ cust, cdnowElog, function(x) x[1] * 10^-3)
names(first) <- c("cust", "first.sales")
cbs <- merge(cbs, first, by = "cust")

#' Estimate with first purchase spend as covariate; see model M2 in Abe (2009)
draws.m2 <- abe.mcmc.DrawParameters(
              cal.cbs = cbs,
              covariates = c("first.sales"),
              mcmc = 5000, burnin = 5000,
              mc.cores = 1)
round(summary(draws.m2$level_2)$quantiles[, c("2.5%", "50%", "97.5%")], 4)
#' -> Parameter Estimates match Table 3 in Abe (2009), except for
#' `log_lambda_first.sales` and `log_mu_first.sales`; note however, that via
#' simulation we can establish that our implementation is able to re-identify
#' the underlying parameters correctly; see
#' `tests/testthat/test-pareto-nbd-abe.R`


x <- readline("Compare predictive performance of the two models (press Enter)")

#' Predict holdout with M1 and M2 models

#' 1) draw future transaction
xstar.m1.draws <- mcmc.DrawFutureTransactions(cbs, draws.m1, T.star = cbs$T.star)
xstar.m2.draws <- mcmc.DrawFutureTransactions(cbs, draws.m2, T.star = cbs$T.star)

#' 2) calculate mean over future transaction draws for each customer
cbs$xstar.m1 <- apply(xstar.m1.draws, 2, mean)
cbs$xstar.m2 <- apply(xstar.m2.draws, 2, mean)

#' 3) compare mean absolute error at individual level
round(c(`MAE without covariates` = mean(abs(cbs$x.star - cbs$xstar.m1)),
        `MAE with covariates`    = mean(abs(cbs$x.star - cbs$xstar.m2))), 4)
