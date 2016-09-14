#' Load CDNow summary data
cbs <- cdnow.sample()$cbs
set.seed(1234)


x <- readline("Estimate Pareto/NBD Abe without covariates (press Enter)")

#' Estimate with no covariates; see model M1 in Abe (2009)
draws.m1 <- abe.mcmc.DrawParameters(cbs, mcmc = 5000, burnin = 5000, mc.cores = 1)
(params.m1 <- round(summary(draws.m1$level_2)$quantiles[, c("2.5%", "50%", "97.5%")], 2))
#' -> Parameter Estimates match Table 3 in Abe (2009);


x <- readline("Estimate Pareto/NBD Abe with covariates (press Enter)")

#' Append dollar amount of first purchase to use as covariate
first <- aggregate(sales ~ cust, cdnow.sample()$elog, head, n = 1)
names(first) <- c("cust", "first.sales")
cbs <- merge(cbs, first, by = "cust")

#' Estimate with first purchase spend as covariate; see model M2 in Abe (2009)
draws.m2 <- abe.mcmc.DrawParameters(cbs, covariates = c("first.sales"), mcmc = 5000, burnin = 5000, mc.cores = 1)
(params.m2 <- round(summary(draws.m2$level_2)$quantiles[, c("2.5%", "50%", "97.5%")], 4))
#' -> Parameter Estimates match Table 3 in Abe (2009), except for 
#' `log_lambda_first` and `log_mu_first`; note however, that via simulation we
#' can establish that the implementation is able to re-identify the underlying
#' parameters correctly


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
