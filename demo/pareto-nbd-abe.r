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


# attempt to calculate marginal log-likelihood; it should be -1'381 according to Abe (2009), but is rather -12'952
ll_mvn <- function(obs, mean, sigma) {
  return(0)
  x1 <- -0.5 * determinant(sigma, logarithm = T)$modulus
  x2 <- -0.5 * (obs-mean) %*% solve(sigma) %*% (obs-mean) 
  x3 <- -0.5 * log(2*pi)
  x1 + x2 + x3
}
ll <- function(x, t.x, T.cal, lambda, mu, tau, z) {
  val <- (1-z)*log(mu) - (lambda+mu) * (z*T.cal + (1-z)*tau)
  if (x > 0) val <- val + x*log(lambda)# + (x-1)*log(t.x) - lgamma(x)
  val
}
a <- matrix(NA, nrow=200, ncol=nrow(cbs))
samples_2 <- as.matrix(draws_m1$level_2)
samples_1 <- lapply(draws_m1$level_1, function(m) as.matrix(m))
for (j in 1:200) {
  cat(j, '\n')
  for (i in 1:nrow(cbs)) {
  cust <- as.list(cbs[i,])
  sample_1 <- as.list(samples_1[[i]][j,])
  sample_2 <- as.list(samples_2[j,])
  a[j, i] <- ll(cust$x, cust$t.x, cust$T.cal, 
               sample_1$lambda, sample_1$mu, sample_1$tau, sample_1$z) +
             ll_mvn(obs = c(log(sample_1$lambda), log(sample_1$mu)), 
                    mean = c(sample_2$log_lambda, sample_2$log_mu),
                    sigma = matrix(c(sample_2$var_log_lambda, 
                                   sample_2$cov_log_lambda_log_mu, 
                                   sample_2$cov_log_lambda_log_mu, 
                                   sample_2$var_log_mu), 
                                  nrow=2))
  }
}
dim(a)
ls <- apply(a, 1, sum)
1/(sum(1/ls)/200)
# -12951.93
mean(ls)
# -12956.58



# Append dollar amount of first purchase; used as covariate in (Abe 2009)
elog.dt <- data.table(elog)
setkey(elog.dt, cust, date)
cbs <- merge(cbs, elog.dt[, list(first=as.numeric(sales[1])*10^-3), by="cust"], by="cust")

set.seed(1)
draws_m2 <- abe.mcmc.DrawParameters(cbs, covariates=c("first"), mcmc=5000, burnin=5000)
plot(draws_m2$level_2, density=FALSE)
params.pnbd_abe.m2 <- round(summary(draws_m2$level_2)$quantiles[, c("2.5%", "50%", "97.5%")], 4)
params.pnbd_abe.m2
# -> Parameter Estimates match (Abe 2009)
