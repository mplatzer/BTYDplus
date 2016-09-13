context("mcmc")

test_that("Pareto/GGG MCMC", {
  cat('test Pareto/GGG\n')
  
  # generate artificial Pareto/GGG data
  params <- list(t = 4.5, gamma = 1.5, r = 0.9, alpha = 10, s = 0.8, beta = 12)
  n <- 100
  expect_silent(pggg.GenerateData(n, 52, c(26, 52), params, TRUE))
  cbs <- pggg.GenerateData(n, 52, 52, params, FALSE)$cbs
  
  # estimate parameters
  draws <- pggg.mcmc.DrawParameters(cbs, mcmc = 100, burnin = 20, thin = 10, chains = 2)
  expect_true(all(c("level_1", "level_2") %in% names(draws)))
  expect_equal(length(draws$level_1), n)
  expect_true(is.mcmc.list(draws$level_1[[1]]))
  expect_true(is.mcmc.list(draws$level_2))
  
  # estimate future transactions
  xstar <- mcmc.DrawFutureTransactions(cbs, draws, T.star = cbs$T.star)
  
  # plot regularity rate
  pggg.plotRegularityRateHeterogeneity(draws)
  
  skip('not run as it takes too long; run this manually')
  skip_on_cran()
  skip_on_travis()
  
  # generate artificial Pareto/GGG data
  params <- list(t = 4.5, gamma = 1.5, r = 0.9, alpha = 10, s = 0.8, beta = 12)
  n <- 5000
  cbs <- pggg.GenerateData(n, 52, 52, params, TRUE)$cbs
  
  # estimate parameters
  draws <- pggg.mcmc.DrawParameters(cbs)
  est <- as.list(summary(draws$level_2)$quantiles[, "50%"])
  
  # require less than 10% deviation in estimated parameters
  ape <- function(act, est) abs(act - est)/act
  expect_lt(ape(params$r, est$r), 0.1)
  expect_lt(ape(params$alpha, est$alpha), 0.1)
  expect_lt(ape(params$s, est$s), 0.2)  # s hard to estimate
  expect_lt(ape(params$beta, est$beta), 0.4)  # beta hard to estimate
  expect_lt(ape(params$t, est$t), 0.1)
  expect_lt(ape(params$gamma, est$gamma), 0.1)
  
  # estimate future transactions & P(alive) & P(active)
  xstar <- mcmc.DrawFutureTransactions(cbs, draws, T.star = cbs$T.star)
  cbs$x.est <- apply(xstar, 2, mean)
  cbs$pactive <- mcmc.PActive(xstar)
  cbs$palive <- mcmc.PAlive(draws)
  
  # require less than 5% deviation
  expect_lt(ape(sum(cbs$x.star), sum(cbs$x.est)), 0.05)
  expect_lt(ape(sum(cbs$palive), sum(cbs$alive)), 0.05)
  expect_lt(ape(sum(cbs$x.star > 0), sum(cbs$pactive)), 0.05)
  
  expect_true(min(cbs$x.star) >= 0)
  expect_true(all(cbs$x.star == round(cbs$x.star)))
  expect_true(all(cbs$palive >= 0 & cbs$palive <= 1))
  
}) 
