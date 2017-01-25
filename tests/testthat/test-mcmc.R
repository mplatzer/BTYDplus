context("mcmc")

test_that("MCMC Helpers", {

  mcmc <- 100
  burnin <- 10
  thin <- 10
  chains <- 2

  # generate small Pareto/NBD cohort
  pnbd_params <- list(r = 0.9, alpha = 10, s = 0.8, beta = 12)
  pnbd_cbs <- pnbd.GenerateData(n = 100, 28, 28, pnbd_params)$cbs
  pnbd_draws <- pnbd.mcmc.DrawParameters(pnbd_cbs, mcmc, burnin, thin, chains)

  # generate small Pareto/Abe cohort
  abe_params <- list()
  abe_params$beta <- matrix(c(0.18, -2.5, 0.5, -0.3, -0.2, 0.8), byrow = TRUE, ncol = 2)
  abe_params$gamma <- matrix(c(0.05, 0.1, 0.1, 0.2), ncol = 2)
  abe_cbs <- abe.GenerateData(n = 100, T.cal = 28, T.star = 28, abe_params)$cbs
  abe_draws <- abe.mcmc.DrawParameters(abe_cbs, c("covariate_1", "covariate_2"), mcmc, burnin, thin, chains)

  # generate small Pareto/GGG cohort
  pggg_params <- list(r = 0.9, alpha = 10, s = 0.8, beta = 12, t = 4, gamma = 2)
  pggg_cbs <- pggg.GenerateData(n = 100, 28, 28, pggg_params)$cbs
  pggg_draws <- pggg.mcmc.DrawParameters(pggg_cbs, mcmc / 10, burnin / 10, thin / 10, chains)

  # draw future transactions
  pnbd_xstar_draws <- mcmc.DrawFutureTransactions(pnbd_cbs, pnbd_draws)
  expect_is(pnbd_xstar_draws, "matrix")
  expect_equal(dim(pnbd_xstar_draws), c(mcmc * chains / thin, nrow(pnbd_cbs)))
  size <- 1
  pnbd_xstar_draws1 <- mcmc.DrawFutureTransactions(pnbd_cbs, pnbd_draws, sample_size = size)
  expect_equal(dim(pnbd_xstar_draws1), c(size, nrow(pnbd_cbs)))
  size <- 500
  pnbd_xstar_draws2 <- mcmc.DrawFutureTransactions(pnbd_cbs, pnbd_draws, sample_size = size)
  expect_equal(dim(pnbd_xstar_draws2), c(size, nrow(pnbd_cbs)))
  expect_gt(cor(apply(pnbd_xstar_draws2, 2, mean), apply(pnbd_xstar_draws, 2, mean)), 0.95)
  expect_silent(pnbd_xstar_draws <- mcmc.DrawFutureTransactions(pnbd_cbs, pnbd_draws, T.star = 10))

  # test setBurnin
  burnin2 <- 30
  pnbd_draws2 <- mcmc.setBurnin(pnbd_draws, burnin = burnin2)
  expect_equal(start(pnbd_draws2$level_2), burnin2)
  expect_equal(start(pnbd_draws2$level_1[[1]]), burnin2)
  pnbd_xstar_draws3 <- mcmc.DrawFutureTransactions(pnbd_cbs, pnbd_draws2)
  expect_equal(dim(pnbd_xstar_draws3), c((mcmc + burnin - burnin2) * chains / thin, nrow(pnbd_cbs)))

  # test Palive
  palive <- mcmc.PAlive(pnbd_draws)
  expect_is(palive, "numeric")
  expect_true(all(palive >= 0))
  expect_true(all(palive <= 1))
  expect_is(mcmc.PAlive(pggg_draws), "numeric")
  expect_is(mcmc.PAlive(abe_draws), "numeric")

  # test Pactive
  pactive <- mcmc.PActive(pnbd_xstar_draws)
  expect_is(pactive, "numeric")
  expect_true(all(pactive >= 0))
  expect_true(all(pactive <= 1))

  # test plotPActiveDiagnostic
  expect_silent(mcmc.plotPActiveDiagnostic(pnbd_cbs, pnbd_xstar_draws))

  # test mcmc.pmf
  expect_equal(dim(mcmc.pmf(pggg_draws, c(28, 56), 0:9)), c(10, 2))
  expect_equal(length(mcmc.pmf(pggg_draws, 56, 0:9)), 10)
  expect_equal(length(mcmc.pmf(pggg_draws, 56, 0)), 1)
  expect_equal(sum(mcmc.pmf(pggg_draws, 2, 0:100)), 1)
  expect_equal(sum(mcmc.pmf(pnbd_draws, 2, 0:100)), 1)
  expect_equal(sum(mcmc.pmf(abe_draws, 2, 0:100)), 1)

  # test mcmc.Expectation
  expect_equal(length(mcmc.Expectation(pnbd_draws, c(28, 56))), 2)
  expect_equal(length(mcmc.Expectation(pggg_draws, c(28, 56))), 2)
  expect_equal(length(mcmc.Expectation(abe_draws, 56)), 1)

  # test against P/NBD simulation
  set.seed(1)
  pnbd_params <- list(r = 0.9, alpha = 10, s = 0.8, beta = 12)
  pnbd_sim <- pnbd.GenerateData(n = 1000, 28, 28, pnbd_params)
  pnbd_elog <- pnbd_sim$elog
  pnbd_cbs <- pnbd_sim$cbs
  pnbd_draws <- pnbd.mcmc.DrawParameters(pnbd_cbs, mcmc = 1000, chains = 1)
  expect_equal(mcmc.Expectation(pnbd_draws, 28), mean(pnbd_cbs$x), tolerance = 0.1)
  x <- mcmc.ExpectedCumulativeTransactions(pnbd_draws, T.cal = pnbd_cbs$T.cal, T.tot = 56, n.periods.final = 56)
  expect_equal(length(x), 56)
  expect_true(all(x >= 0))
  expect_equal(x[28], sum(pnbd_cbs$x), tolerance = 0.1)
  expect_equal(x[56], sum(pnbd_cbs$x) +  sum(pnbd_cbs$x.star), tolerance = 0.1)
  pnbd_cum <- elog2cum(pnbd_elog, by = 1)[-1]
  mat <- mcmc.PlotTrackingCum(pnbd_draws, pnbd_cbs$T.cal, T.tot = 56, pnbd_cum)
  expect_equal(mat[1, ], mat[2, ], tolerance = 0.1)
  pnbd_inc <- elog2inc(pnbd_elog, by = 1)
  mat <- mcmc.PlotTrackingInc(pnbd_draws, pnbd_cbs$T.cal, T.tot = 56, pnbd_inc)
  expect_equal(mat[1, ], mat[2, ], tolerance = 0.2)

  expect_equal(mcmc.PlotFrequencyInCalibration(pnbd_draws, pnbd_cbs),
               BTYD::pnbd.PlotFrequencyInCalibration(unlist(pnbd_params), pnbd_cbs, censor = 7),
               tolerance = 0.1)
})
