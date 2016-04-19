require(data.table)

# Load CDNow event log
data(cdnowElog, package = "BTYD")
elog <- data.frame(t(sapply(2:nrow(cdnowElog), function(i) strsplit(as.character(cdnowElog[i, ]), split = ",")[[1]])), 
  stringsAsFactors = FALSE)
names(elog) <- c("cust", "sampleid", "date", "cds", "sales")
elog$date <- as.Date(elog$date, "%Y%m%d")
elog$sales <- as.numeric(elog$sales)

# Transform to CBS (including extra column 'litt')
cbs <- elog2cbs(elog, per = "week", T.cal = as.Date("1997-09-30"))

# Estimate with no covariates, i.e. model M1 in (Abe 2009)
set.seed(1)
draws_m1 <- abe.mcmc.DrawParameters(cbs, mcmc = 5000, burnin = 5000)
plot(draws_m1$level_2, density = FALSE)
params.pnbd_abe.m1 <- round(summary(draws_m1$level_2)$quantiles[, c("2.5%", "50%", "97.5%")], 2)
params.pnbd_abe.m1
# -> Parameter Estimates match (Abe 2009) 2.5% 50% 97.5% log_lambda -3.74 -3.57 -3.37 log_mu -3.98 -3.74 -3.42
# var_log_lambda 1.19 1.42 1.69 cov_log_lambda_log_mu -0.11 0.29 0.76 var_log_mu 1.35 2.66 5.17


# Append dollar amount of first purchase; used as covariate in (Abe 2009)
elog.dt <- data.table(elog)
setkey(elog.dt, cust, date)
cbs <- merge(cbs, elog.dt[, list(first = as.numeric(sales[1]) * 10^-3), by = "cust"], by = "cust")

set.seed(1)
draws_m2 <- abe.mcmc.DrawParameters(cbs, covariates = c("first"), mcmc = 5000, burnin = 5000)
plot(draws_m2$level_2, density = FALSE)
params.pnbd_abe.m2 <- round(summary(draws_m2$level_2)$quantiles[, c("2.5%", "50%", "97.5%")], 4)
params.pnbd_abe.m2
# -> Parameter Estimates match Abe (2009); except for log_lambda_first and log_mu_first; note however, that via
# simulation we can establish that the implementation is able to re-identify the underlying parameters
# correctly 2.5% 50% 97.5% log_lambda_intercept -3.9623 -3.7172 -3.5192 log_mu_intercept -4.0323 -3.6350
# -3.2458 log_lambda_first 2.4690 6.1411 8.9629 log_mu_first -11.0806 1.6959 7.9268 var_log_lambda 1.0881
# 1.3238 1.5536 cov_log_lambda_log_mu -0.2494 0.1835 0.6724 var_log_mu 1.0430 3.0687 5.2276 
