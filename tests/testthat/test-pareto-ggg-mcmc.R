context("mcmc")

test_that("Pareto/GGG MCMC", {

  # generate artificial Pareto/GGG data
  params <- list(t = 4.5, gamma = 1.5, r = 0.9, alpha = 10, s = 0.8, beta = 12)
  n <- 100
  expect_silent(pggg.GenerateData(n, 52, c(26, 52), params))
  cbs <- pggg.GenerateData(n, 52, 52, params)$cbs

  # estimate parameters
  draws <- pggg.mcmc.DrawParameters(as.data.table(cbs),
                                    mcmc = 10, burnin = 0, thin = 1, mc.cores = 1,
                                    param_init = list(r = 1, alpha = 1, s = 1, beta = 1, t = 1, gamma = 1))
  draws <- pggg.mcmc.DrawParameters(cbs,
                                    mcmc = 100, burnin = 20, thin = 10, chains = 2, mc.cores = 1,
                                    param_init = params)
  expect_true(all(c("level_1", "level_2") %in% names(draws)))
  expect_equal(length(draws$level_1), n)
  expect_true(coda::is.mcmc.list(draws$level_1[[1]]))
  expect_true(coda::is.mcmc.list(draws$level_2))

  # estimate future transactions
  xstar <- mcmc.DrawFutureTransactions(cbs, draws, T.star = cbs$T.star)

  # plot regularity rate
  pggg.plotRegularityRateHeterogeneity(draws)

  skip("skip long-running test of Pareto/GGG parameter recovery")

  # generate artificial Pareto/GGG data
  set.seed(1)
  params <- list(t = 4.5, gamma = 1.5, r = 0.9, alpha = 5, s = 0.8, beta = 12)
  n <- 5000
  cbs <- pggg.GenerateData(n,
                           round(runif(n, 36, 96) / 12) * 12,
                           36,
                           params)$cbs

  # estimate parameters
  draws <- pggg.mcmc.DrawParameters(cbs, mc.cores = 1)
  est <- as.list(summary(draws$level_2)$quantiles[, "50%"])

  # require less than 10% deviation in estimated parameters
  expect_equal(params, est, tolerance = 0.1)

  # estimate future transactions & P(alive) & P(active)
  xstar <- mcmc.DrawFutureTransactions(cbs, draws, T.star = cbs$T.star)
  cbs$x.est <- apply(xstar, 2, mean)
  cbs$pactive <- mcmc.PActive(xstar)
  cbs$palive <- mcmc.PAlive(draws)

  # require less than 5% deviation
  expect_equal(sum(cbs$x.star), sum(cbs$x.est), tolerance = 0.05)
  expect_equal(sum(cbs$palive), sum(cbs$alive), tolerance = 0.05)
  expect_equal(sum(cbs$x.star > 0), sum(cbs$pactive), tolerance = 0.05)

  expect_true(min(cbs$x.star) >= 0)
  expect_true(all(cbs$x.star == round(cbs$x.star)))
  expect_true(all(cbs$palive >= 0 & cbs$palive <= 1))

})
