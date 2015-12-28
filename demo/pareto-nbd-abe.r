require(data.table)

# Load CDNow event log
data(cdnowElog, package="BTYD")
elog <- data.frame(t(sapply(2:nrow(cdnowElog), function(i) strsplit(as.character(cdnowElog[i,]), split=",")[[1]])), stringsAsFactors=FALSE)
names(elog) <- c("cust", "sampleid", "date", "cds", "sales")
elog$date <- as.Date(elog$date, "%Y%m%d")
elog$sales <- as.numeric(elog$sales)

# Transform to CBS (including extra column 'litt')
cbs <- elog2cbs(elog, per="week", T.cal=as.Date("1997-09-30"))

# Estimate with no covariates, i.e. model M1 in (Abe 2009)
set.seed(1)
draws_m1 <- abe.mcmc.DrawParameters(cbs, mcmc=5000, burnin=5000)
plot(draws_m1$level_2, density=FALSE)
params.pnbd_abe.m1 <- round(summary(draws_m1$level_2)$quantiles[, c("2.5%", "50%", "97.5%")], 2)
params.pnbd_abe.m1
# -> Parameter Estimates match (Abe 2009)
#                        2.5%   50% 97.5%
# log_lambda            -3.74 -3.57 -3.37
# log_mu                -3.98 -3.74 -3.42
# var_log_lambda         1.19  1.42  1.69
# cov_log_lambda_log_mu -0.11  0.29  0.76
# var_log_mu             1.35  2.66  5.17


# Append dollar amount of first purchase; used as covariate in (Abe 2009)
elog.dt <- data.table(elog)
setkey(elog.dt, cust, date)
cbs <- merge(cbs, elog.dt[, list(first=as.numeric(sales[1])*10^-3), by="cust"], by="cust")

set.seed(1)
draws_m2 <- abe.mcmc.DrawParameters(cbs, covariates=c("first"), mcmc=5000, burnin=5000)
plot(draws_m2$level_2, density=FALSE)
params.pnbd_abe.m2 <- round(summary(draws_m2$level_2)$quantiles[, c("2.5%", "50%", "97.5%")], 4)
params.pnbd_abe.m2
# -> Parameter Estimates match Abe (2009); except for log_lambda_first and
# log_mu_first; note however, that via simulation we can establish that the
# implementation is able to re-identify the underlying parameters correctly
#                           2.5%     50%   97.5%
# log_lambda_intercept   -3.9623 -3.7172 -3.5192
# log_mu_intercept       -4.0323 -3.6350 -3.2458
# log_lambda_first        2.4690  6.1411  8.9629
# log_mu_first          -11.0806  1.6959  7.9268
# var_log_lambda          1.0881  1.3238  1.5536
# cov_log_lambda_log_mu  -0.2494  0.1835  0.6724
# var_log_mu              1.0430  3.0687  5.2276




## attempt to calculate marginal log-likelihood; it should be -1'381 according
## to Table 3 in Abe (2009), but with below code we get -18'292 instead !?

# log-likelihood function of multivariate normal; https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Likelihood_function
ll_mvn <- function(obs, mean, sigma) {
  x1 <- -0.5 * as.numeric(determinant(sigma, logarithm = T)$modulus)
  x2 <- -0.5 * (obs-mean) %*% solve(sigma) %*% (obs-mean) 
  x3 <- -0.5 * length(mean) * log(2*pi)
  x1 + x2 + x3
}

# log-likelihood function of individual customer; equation (5) from Abe (2009)
ll_ind <- function(x, t.x, T.cal, lambda, mu, tau, z) {
  val <- (1-z)*log(mu) - (lambda+mu) * (z*T.cal + (1-z)*tau)
  if (x > 0) val <- val + x*log(lambda) + (x-1)*log(t.x) - lgamma(x)
  val
}

# convert level 1 and level 2 MCMC draws into matrix
samples_1   <- lapply(draws_m1$level_1, function(m) as.matrix(m))
samples_2   <- as.matrix(draws_m1$level_2)

# calculate for each MCMC draw the log-likelihood across all customers
no_of_cust  <- length(samples_1)
no_of_draws <- nrow(samples_2)
LL <- c()
for (j in 1:no_of_draws) {
  cat('calculate log-likelihood for ', j, '-th draw\n')
  LL_j <- c()
  for (i in 1:no_of_cust) {
    cust     <- as.list(cbs[i,])            # the data for i-th customer
    sample_1 <- as.list(samples_1[[i]][j,]) # j-th draw for i-th customer
    sample_2 <- as.list(samples_2[j,])      # j-th draw of heterogeneity parameters
    LL_j[i] <- ll_ind(x = cust$x, t.x = cust$t.x, T.cal = cust$T.cal,
                      lambda = sample_1$lambda, mu = sample_1$mu, tau = sample_1$tau, z = sample_1$z) + 
               ll_mvn(obs   = c(log(sample_1$lambda), log(sample_1$mu)),
                      mean  = c(sample_2$log_lambda, sample_2$log_mu),
                      sigma = matrix(c(sample_2$var_log_lambda, 
                                       sample_2$cov_log_lambda_log_mu, 
                                       sample_2$cov_log_lambda_log_mu, 
                                       sample_2$var_log_mu), nrow=2))
  }
  # sum across customers to get total log-likelihood for the j-th draw
  LL[j] <-  sum(LL_j)
}

# harmonic mean of likelihood; eq (14) from Newton & Raftery (1994) J.Royal Stat. Soc.
mean(LL) - log(mean(exp(mean(LL) - LL)))
# -18'292.42


