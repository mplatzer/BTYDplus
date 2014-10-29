
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
