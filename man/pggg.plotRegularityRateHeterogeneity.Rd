% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pareto-ggg-mcmc.R
\name{pggg.plotRegularityRateHeterogeneity}
\alias{pggg.plotRegularityRateHeterogeneity}
\title{Pareto/GGG Plot Regularity Rate Heterogeneity}
\usage{
pggg.plotRegularityRateHeterogeneity(
  draws,
  xmax = NULL,
  fn = NULL,
  title = "Distribution of Regularity Rate k"
)
}
\arguments{
\item{draws}{MCMC draws as returned by \code{\link{pggg.mcmc.DrawParameters}}.}

\item{xmax}{Upper bound for x-scale.}

\item{fn}{Optional function to summarize individual-level draws for k, e.g. 'mean'.}

\item{title}{Plot title.}
}
\description{
Plots and returns the estimated gamma distribution of k (customers'
regularity in interpurchase times).
}
\examples{
data("groceryElog")
cbs <- elog2cbs(groceryElog, T.cal = "2006-12-31")
param.draws <- pggg.mcmc.DrawParameters(cbs,
  mcmc = 20, burnin = 10, thin = 2, chains = 1) # short MCMC to run demo fast
pggg.plotRegularityRateHeterogeneity(param.draws)
}
\references{
Platzer, M., & Reutterer, T. (2016). Ticking away the moments:
  Timing regularity helps to better predict customer activity. Marketing
  Science, 35(5), 779-799. \doi{10.1287/mksc.2015.0963}
}
