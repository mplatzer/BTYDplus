context("mle")

test_that("Pareto/GGG MCMC", {

  # generate artificial Pareto/GGG data 
  params <- list(t=4.5, gamma=1.5, r=0.9, alpha=10, s=0.8, beta=12)
  n <- 5000
  cbs <- pggg.GenerateData(n, 52, 52, params, TRUE)$cbs
  
  # estimate parameters
  draws <- pggg.mcmc.DrawParameters(cbs)
  est <- as.list(summary(draws$level_2)$quantiles[, "50%"])
  
  expect_true(all(c("level_1", "level_2") %in% names(draws)))
  expect_equal(length(draws$level_1), n)
  expect_true(is.mcmc.list(draws$level_1[[1]]))
  expect_true(is.mcmc.list(draws$level_2))
  
  # require less than 10% deviation in estimated parameters
  ape <- function(act, est) abs(act-est)/act
  expect_less_than(ape(params$r, est$r), 0.10)
  expect_less_than(ape(params$alpha, est$alpha), 0.10)
  expect_less_than(ape(params$s, est$s), 0.20) # s hard to estimate
  expect_less_than(ape(params$beta, est$beta), 0.40) # beta hard to estimate
  expect_less_than(ape(params$t, est$t), 0.10)
  expect_less_than(ape(params$gamma, est$gamma), 0.10)
  
  # estimate future transactions & P(alive)
  xstar <- mcmc.DrawFutureTransactions(cbs, draws, T.star=cbs$T.star)
  cbs$x.est <- apply(xstar, 2, mean)
  cbs$pactive <- mcmc.PActive(xstar)
  cbs$palive <- mcmc.PAlive(cbs, draws)

  # require less than 5% deviation
  expect_less_than(ape(sum(cbs$x.star), sum(cbs$x.est)), 0.05)
  expect_less_than(ape(sum(cbs$palive), sum(cbs$alive)), 0.05)
  expect_less_than(ape(sum(cbs$x.star>0), sum(cbs$pactive)), 0.05)

  expect_true(min(cbs$x.star)>=0)
  expect_true(all(cbs$x.star==round(cbs$x.star)))
  expect_true(all(cbs$palive>=0 & cbs$palive<=1))
  
  draws2 <- mcmc.setBurnin(draws, burnin=600)
  expect_equal(start(draws2$level_2), 600)
  expect_equal(start(draws2$level_1[[1]]), 600)
  
  mcmc.plotPActiveDiagnostic(cbs, xstar)
  pggg.mcmc.plotRegularityRateHeterogeneity(draws)
}
