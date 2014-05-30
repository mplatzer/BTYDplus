
# Load CDNow event log
data(cdnowElog, package="BTYD")
elog <- data.frame(t(sapply(2:nrow(cdnowElog), function(i) strsplit(as.character(cdnowElog[i,]), split=",")[[1]])))
names(elog) <- c("cust", "sampleid", "date", "cds", "sales")
elog$date <- as.Date(elog$date, "%Y%m%d")

# Transform to CBS (including extra column 'litt')
cbs <- elog2cbs(elog, per="week", T.cal=as.Date("1997-09-30"))

# Pareto/NBD MCMC
set.seed(1)
pnbd.draws <- pnbd.mcmc.DrawParameters(cbs, mcmc=1000, burnin=1000, chains=2, thin=10)
summary(pnbd.draws$level_2)
plot(pnbd.draws$level_2)
summary(pnbd.draws$level_2)$quantiles[, "50%"]

# Pareto/CNBD MCMC
set.seed(2)
param_init <- params.pnbd.mcmc
param_init$t <- 10
param_init$gamma <- 10
pcnbd.draws <- pcnbd.mcmc.DrawParameters(cbs, mcmc=100, burnin=100, chains=2, thin=2, param_init=param_init)
# Note: MCMC chain should be run longer; for demo purposes we keep runtime short here

plot(pcnbd.draws$level_2, density=FALSE)
plot(pcnbd.draws$level_2, trace=FALSE)

coda::gelman.diag(pcnbd.draws$level_2)
# -> MCMC chains have not converged yet

summary(pcnbd.draws$level_2)$quantiles[, "50%"]

pcnbd.mcmc.plotRegularityRateHeterogeneity(pcnbd.draws)
# -> very narrow distribution around k=1; 
# CDNow customers purchases indeed seem to follow Poisson process

round(effectiveSize(pcnbd.draws$level_2))
# -> effective sample size are small for such a short chain

# draw future transaction
pnbd.xstar <- pcnbd.mcmc.DrawFutureTransactions(cbs, pnbd.draws, T.star=cbs$T.star)
pcnbd.xstar <- pcnbd.mcmc.DrawFutureTransactions(cbs, pcnbd.draws, T.star=cbs$T.star)

# calculate mean over future transaction draws for each customer
cbs$pnbd.mcmc <- apply(pnbd.xstar, 2, mean)
cbs$pcnbd.mcmc <- apply(pcnbd.xstar, 2, mean)

# calculate P(active)
cbs$pactive.pnbd.mcmc <- apply(pnbd.xstar, 2, function(x) mean(x>0))
cbs$pactive.pcnbd.mcmc <- apply(pcnbd.xstar, 2, function(x) mean(x>0))
