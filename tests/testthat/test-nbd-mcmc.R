context("mcmc")

test_that("NBD MCMC", {

  # generate artificial NBD data
  set.seed(1)
  params <- list(r = 0.9, alpha = 5)
  n <- 5000
  sim <- nbd.GenerateData(n,
                          round(runif(n, 36, 96) / 12) * 12,
                          36,
                          params)
  cbs <- sim$cbs

  # test basic parameter estimation
  draws <- nbd.mcmc.DrawParameters(as.data.table(cbs),
                                   mcmc = 10, burnin = 0, thin = 1, mc.cores = 1,
                                   param_init = list(r = 1, alpha = 1))

  # test parameter recovery
  draws <- nbd.mcmc.DrawParameters(cbs, chains = 1, mc.cores = 1, trace = 1000)
  est <- as.list(summary(draws$level_2)$quantiles[, "50%"])

  expect_true(all(c("level_1", "level_2") %in% names(draws)))
  expect_equal(length(draws$level_1), n)
  expect_true(coda::is.mcmc.list(draws$level_1[[1]]))
  expect_true(coda::is.mcmc.list(draws$level_2))

  # require less than 10% deviation in estimated parameters
  expect_equal(params, est, tolerance = 0.10)

  # estimate future transactions & P(alive)
  xstar <- mcmc.DrawFutureTransactions(cbs, draws, T.star = cbs$T.star)
  cbs$x.est <- apply(xstar, 2, mean)
  cbs$pactive <- apply(xstar, 2, function(x) mean(x > 0))
  cbs$palive <- mcmc.PAlive(draws)

  # require less than 5% deviation
  expect_equal(sum(cbs$x.star), sum(cbs$x.est), tolerance = 0.05)

  expect_true(min(cbs$x.star) >= 0)
  expect_true(all(cbs$x.star == round(cbs$x.star)))
  expect_true(all(cbs$palive >= 0 & cbs$palive <= 1))

})
