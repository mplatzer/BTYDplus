context("mcmc")

test_that("MCMC Helpers", {

  mcmc <- 100
  burnin <- 10
  thin <- 10
  chains <- 2

  # generate small Pareto/NBD cohort
  pnbd.params <- list(r = 0.9, alpha = 10, s = 0.8, beta = 12)
  pnbd.cbs <- pnbd.GenerateData(n = 100, 28, 28, pnbd.params)$cbs
  pnbd.draws <- pnbd.mcmc.DrawParameters(pnbd.cbs, mcmc, burnin, thin, chains)

  # generate small Pareto/Abe cohort
  abe.params <- list()
  abe.params$beta <- matrix(c(0.18, -2.5, 0.5, -0.3, -0.2, 0.8), byrow = TRUE, ncol = 2)
  abe.params$gamma <- matrix(c(0.05, 0.1, 0.1, 0.2), ncol = 2)
  abe.cbs <- abe.GenerateData(n = 100, T.cal = 28, T.star = 28, abe.params)$cbs
  abe.draws <- abe.mcmc.DrawParameters(abe.cbs, c("covariate_1", "covariate_2"), mcmc, burnin, thin, chains)

  # generate small Pareto/GGG cohort
  pggg.params <- list(r = 0.9, alpha = 10, s = 0.8, beta = 12, t = 4, gamma = 2)
  pggg.cbs <- pggg.GenerateData(n = 100, 28, 28, pggg.params)$cbs
  pggg.draws <- pggg.mcmc.DrawParameters(pggg.cbs, mcmc / 10, burnin / 10, thin / 10, chains)

  # draw future transactions
  pnbd.xstar.draws <- mcmc.DrawFutureTransactions(pnbd.cbs, pnbd.draws)
  expect_is(pnbd.xstar.draws, "matrix")
  expect_equal(dim(pnbd.xstar.draws), c(mcmc * chains / thin, nrow(pnbd.cbs)))
  size <- 1
  pnbd.xstar.draws1 <- mcmc.DrawFutureTransactions(pnbd.cbs, pnbd.draws, sample_size = size)
  expect_equal(dim(pnbd.xstar.draws1), c(size, nrow(pnbd.cbs)))
  size <- 500
  pnbd.xstar.draws2 <- mcmc.DrawFutureTransactions(pnbd.cbs, pnbd.draws, sample_size = size)
  expect_equal(dim(pnbd.xstar.draws2), c(size, nrow(pnbd.cbs)))
  expect_gt(cor(apply(pnbd.xstar.draws2, 2, mean), apply(pnbd.xstar.draws, 2, mean)), 0.95)
  expect_silent(pnbd.xstar.draws <- mcmc.DrawFutureTransactions(pnbd.cbs, pnbd.draws, T.star = 10))

  # test setBurnin
  burnin2 <- 30
  pnbd.draws2 <- mcmc.setBurnin(pnbd.draws, burnin = burnin2)
  expect_equal(start(pnbd.draws2$level_2), burnin2)
  expect_equal(start(pnbd.draws2$level_1[[1]]), burnin2)
  pnbd.xstar.draws3 <- mcmc.DrawFutureTransactions(pnbd.cbs, pnbd.draws2)
  expect_equal(dim(pnbd.xstar.draws3), c((mcmc + burnin - burnin2) * chains / thin, nrow(pnbd.cbs)))

  # test Palive
  palive <- mcmc.PAlive(pnbd.draws)
  expect_is(palive, "numeric")
  expect_true(all(palive >= 0))
  expect_true(all(palive <= 1))
  expect_is(mcmc.PAlive(pggg.draws), "numeric")
  expect_is(mcmc.PAlive(abe.draws), "numeric")

  # test Pactive
  pactive <- mcmc.PActive(pnbd.xstar.draws)
  expect_is(pactive, "numeric")
  expect_true(all(pactive >= 0))
  expect_true(all(pactive <= 1))

  # test plotPActiveDiagnostic
  expect_silent(mcmc.plotPActiveDiagnostic(pnbd.cbs, pnbd.xstar.draws))

  # test mcmc.pmf
  expect_equal(dim(mcmc.pmf(pggg.draws, c(28, 56), 0:9)), c(10, 2))
  expect_equal(length(mcmc.pmf(pggg.draws, 56, 0:9)), 10)
  expect_equal(length(mcmc.pmf(pggg.draws, 56, 0)), 1)
  expect_equal(sum(mcmc.pmf(pggg.draws, 2, 0:100)), 1)
  expect_equal(sum(mcmc.pmf(pnbd.draws, 2, 0:100)), 1)
  expect_equal(sum(mcmc.pmf(abe.draws, 2, 0:100)), 1)

  # test mcmc.Expectation
  expect_equal(length(mcmc.Expectation(pnbd.draws, c(28, 56))), 2)
  expect_equal(length(mcmc.Expectation(pggg.draws, c(28, 56))), 2)
  expect_equal(length(mcmc.Expectation(abe.draws, 56)), 1)

  # test against P/NBD simulation
  set.seed(1)
  pnbd.params <- list(r = 0.9, alpha = 10, s = 0.8, beta = 12)
  pnbd.sim <- pnbd.GenerateData(n = 1000, 28, 28, pnbd.params)
  pnbd.elog <- pnbd.sim$elog
  pnbd.cbs <- pnbd.sim$cbs
  pnbd.draws <- pnbd.mcmc.DrawParameters(pnbd.cbs, mcmc = 1000, chains = 1)
  expect_equal(mcmc.Expectation(pnbd.draws, 28), mean(pnbd.cbs$x), tolerance = 0.1)
  x <- mcmc.ExpectedCumulativeTransactions(pnbd.draws, T.cal = pnbd.cbs$T.cal, T.tot = 56, n.periods.final = 56)
  expect_equal(length(x), 56)
  expect_true(all(x >= 0))
  expect_equal(x[28], sum(pnbd.cbs$x), tolerance = 0.1)
  expect_equal(x[56], sum(pnbd.cbs$x) +  sum(pnbd.cbs$x.star), tolerance = 0.1)
  pnbd.cum <- elog2cum(pnbd.elog, by = 1)[-1]
  mat <- mcmc.PlotTrackingCum(pnbd.draws, pnbd.cbs$T.cal, T.tot = 56, pnbd.cum)
  expect_equal(mat[1, ], mat[2, ], tolerance = 0.1)
  pnbd.inc <- elog2inc(pnbd.elog, by = 1)
  mat <- mcmc.PlotTrackingInc(pnbd.draws, pnbd.cbs$T.cal, T.tot = 56, pnbd.inc)
  expect_equal(mat[1, ], mat[2, ], tolerance = 0.2)

  expect_equal(mcmc.PlotFrequencyInCalibration(pnbd.draws, pnbd.cbs),
               BTYD::pnbd.PlotFrequencyInCalibration(unlist(pnbd.params), pnbd.cbs, censor = 7),
               tolerance = 0.1)
})
