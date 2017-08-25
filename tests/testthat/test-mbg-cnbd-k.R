context("mle")

test_that("MBG/CNBD-k", {

  # generate artificial MBG/CNBD-k data
  set.seed(1)
  n <- 8000
  params <- c(k = 3, r = 0.85, alpha = 1.45, a = 0.79, b = 2.42)
  sim <- mbgcnbd.GenerateData(n,
                             round(runif(n, 36, 96) / 12) * 12,
                             36,
                             params,
                             "2010-01-01")
  cbs <- sim$cbs
  elog <- sim$elog

  # estimate parameters with fixed k
  est1 <- mbgnbd.EstimateParameters(cbs[, c("x", "t.x", "T.cal", "litt")])
  expect_equal(as.integer(est1[1]), 1L)
  est2 <- mbgcnbd.EstimateParameters(cbs[, c("x", "t.x", "T.cal", "litt")], k = 2)
  expect_equal(as.integer(est2[1]), 2L)

  # estimate parameters with unspecified k
  est <- mbgcnbd.EstimateParameters(cbs[, c("x", "t.x", "T.cal", "litt")])
  expect_equal(est, params, tolerance = 0.1)

  # estimate future transactions & P(alive) with true parameters
  cbs$x.est <- mbgcnbd.ConditionalExpectedTransactions(params, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)
  cbs$palive <- mbgcnbd.PAlive(params, cbs$x, cbs$t.x, cbs$T.cal)

  # require less than 5% deviation
  expect_equal(sum(cbs$x.star), sum(cbs$x.est), tolerance = 0.05)
  expect_equal(sum(cbs$palive), sum(cbs$alive), tolerance = 0.05)

  expect_true(min(cbs$x.star) >= 0)
  expect_true(all(cbs$x.star == round(cbs$x.star)))
  expect_true(all(cbs$palive >= 0 & cbs$palive <= 1))

  expect_equal(dim(mbgcnbd.pmf(params, c(28, 56), 0:9)), c(10, 2))
  expect_equal(length(mbgcnbd.pmf(params, 56, 0:9)), 10)
  expect_equal(length(mbgcnbd.pmf(params, 56, 0)), 1)
  expect_equal(sum(mbgcnbd.pmf(params, 2, 0:100)), 1)
  expect_silent(exp <- mbgcnbd.Expectation(params, 11))
  cum <- mbgcnbd.ExpectedCumulativeTransactions(params, 11, 39, 12)
  expect_true(all(diff(cum) > 0))
  expect_equal(length(cum), 12)
  expect_silent(mbgcnbd.PlotTrackingInc(params, cbs$T.cal, max(cbs$T.cal + cbs$T.star), elog2inc(elog, by = 14)))
  expect_silent(mbgcnbd.PlotTrackingCum(params, cbs$T.cal, max(cbs$T.cal + cbs$T.star), elog2cum(elog, by = 14)))
  mat <- mbgcnbd.PlotFrequencyInCalibration(params, cbs, 7)
  expect_equal(mat[1, ], mat[2, ], tolerance = 0.1)
  mat <- mbgcnbd.PlotFreqVsConditionalExpectedFrequency(params, cbs$T.star, cbs, cbs$x.star, 7)
  expect_equal(mat[1, ], mat[2, ], tolerance = 0.1)

  # check that bias correction does not screw up single estimates
  expect_lt(mbgcnbd.ConditionalExpectedTransactions(params, T.star = 32, x = c(0, 1), t.x = c(0, 12), T.cal = 32)[2], 1)
  expect_true(mbgcnbd.ConditionalExpectedTransactions(params, T.star = 32, x = 3, t.x = 12, T.cal = 32) !=
                mbgcnbd.ConditionalExpectedTransactions(params, T.star = 32, x = 2, t.x = 12, T.cal = 32))

  # check that bias correction actually results in unbiased estimates
  skip("skip long-running test of bias correction")
  set.seed(1)
  bias <- replicate(40, {
    cbs <- mbgcnbd.GenerateData(n = n,
                                T.cal = round(runif(n, 12, 96) / 12) * 12,
                                T.star = 32,
                                params = params)$cbs
    cbs$x.est <- mbgcnbd.ConditionalExpectedTransactions(params, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)
    sum(cbs$x.est) / sum(cbs$x.star)
  })
  expect_equal(mean(bias), 1, 0.01)
})
