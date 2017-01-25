context("mle")

test_that("BG/CNBD-k", {

  # validate against BTYD implementation
  set.seed(1)
  params <- c(1, 0.85, 1.45, 0.79, 2.42)
  n <- 500
  date.zero <- "2010-01-01"
  sim <- bgcnbd.GenerateData(n,
                             round(runif(n, 36, 96) / 12) * 12,
                             36,
                             params,
                             date.zero)
  cbs  <- sim$cbs
  elog <- sim$elog
  expect_is(elog$date, "POSIXct")
  expect_is(cbs$first, "POSIXct")
  expect_equal(min(sim$elog$date), as.POSIXct(date.zero))

  params_est_btyd <- BTYD::bgnbd.EstimateParameters(cbs)
  params_est_btyd_plus <- bgcnbd.EstimateParameters(cbs, k = 1)[-1]
  expect_equal(unname(round(params_est_btyd, 2)),
               unname(round(params_est_btyd_plus, 2)))
  expect_equal(BTYD::bgnbd.PAlive(params[-1], 0, 0, 32),
               bgcnbd.PAlive(params, 0, 0, 32))
  expect_equal(BTYD::bgnbd.PAlive(params[-1], 1, 16, 32),
               bgcnbd.PAlive(params, 1, 16, 32))
  expect_equal(unname(BTYD::bgnbd.Expectation(params[-1], 1:3)),
               unname(bgcnbd.Expectation(params, 1:3)))
  expect_equal(BTYD::bgnbd.ConditionalExpectedTransactions(params[-1], 32, 1, 16, 32),
               bgcnbd.ConditionalExpectedTransactions(params, 32, 1, 16, 32))
  expect_equal(unname(BTYD::bgnbd.pmf(params[-1], 32, 0:2)),
               unname(bgcnbd.pmf(params, 32, 0:2)))
  expect_equal(BTYD::bgnbd.PlotFrequencyInCalibration(params[-1], cbs, 7),
               bgcnbd.PlotFrequencyInCalibration(params, cbs, 7), tolerance = 0.01)
  expect_equal(BTYD::bgnbd.PlotFreqVsConditionalExpectedFrequency(params[-1], T.star = 39, cbs, cbs$x.star, 7),
               bgcnbd.PlotFreqVsConditionalExpectedFrequency(params, T.star = 39, cbs, cbs$x.star, 7), tolerance = 0.01)
  # Note: BTYD::bgnbd.PlotRecVsConditionalExpectedFrequency can't handle missing bins, so we can't compare
  expect_silent(bgcnbd.PlotRecVsConditionalExpectedFrequency(params, cbs, T.star = 39, cbs$x.star))
  inc_tracking <- elog2inc(elog, by = 7)
  expect_equal(BTYD::bgnbd.PlotTrackingInc(params[-1], cbs$T.cal, max(cbs$T.cal) + 32, inc_tracking),
               bgcnbd.PlotTrackingInc(params, cbs$T.cal, max(cbs$T.cal) + 32, inc_tracking),
               tolerance = 0.01)
  expect_silent(bgcnbd.PlotTrackingInc(params, cbs$T.cal, max(cbs$T.cal) + 32, inc_tracking,
                                       xticklab = 1:length(inc_tracking)))
  expect_equal(BTYD::bgnbd.ExpectedCumulativeTransactions(params[-1], 11, 39, 12),
               bgcnbd.ExpectedCumulativeTransactions(params, 11, 39, 12),
               tolerance = 0.01)
  cu_tracking <- cumsum(inc_tracking)
  expect_equal(BTYD::bgnbd.PlotTrackingCum(params[-1], cbs$T.cal, 32 + 32, cu_tracking),
               bgcnbd.PlotTrackingCum(params, cbs$T.cal, 32 + 32, cu_tracking),
               tolerance = 0.01)

  # generate artificial BG/CNBD-k data
  set.seed(1)
  n <- 1000
  params <- c(k = 3, r = 0.85, alpha = 1.45, a = 0.79, b = 2.42)
  sim <- bgcnbd.GenerateData(n,
                              round(runif(n, 36, 96) / 12) * 12,
                              c(32, 64),
                              params = params)
  cbs <- sim$cbs
  elog <- sim$elog

  # estimate regularity & parameters
  k_est <- estimateRegularity(elog)
  est <- bgcnbd.EstimateParameters(cbs[, c("x", "t.x", "T.cal", "litt")])
  est_fixed_k <- bgcnbd.EstimateParameters(cbs[, c("x", "t.x", "T.cal", "litt")], k = params[1])

  expect_equal(params[1], est[1])
  expect_equal(est, est_fixed_k)
  expect_equal(params, est, tolerance = 0.1)

  # estimate future transactions & P(alive) with true parameters
  cbs$x.est32 <- bgcnbd.ConditionalExpectedTransactions(params, 32, cbs$x, cbs$t.x, cbs$T.cal)
  cbs$x.est64 <- bgcnbd.ConditionalExpectedTransactions(params, 64, cbs$x, cbs$t.x, cbs$T.cal)
  cbs$palive <- bgcnbd.PAlive(params, cbs$x, cbs$t.x, cbs$T.cal)

  # require less than 5% deviation
  expect_equal(sum(cbs$x.star32), sum(cbs$x.est32), tolerance = 0.05)
  expect_equal(sum(cbs$x.star64), sum(cbs$x.est64), tolerance = 0.05)
  expect_equal(sum(cbs$palive), sum(cbs$alive), tolerance = 0.05)

  # some basic sanity checks
  expect_true(min(cbs$x.star32) >= 0)
  expect_true(all(cbs$x.star32 == round(cbs$x.star32)))
  expect_true(all(cbs$palive >= 0 & cbs$palive <= 1))

  # test bgcnbd.pmf
  expect_equal(dim(bgcnbd.pmf(params, c(28, 56), 0:9)), c(10, 2))
  expect_equal(length(bgcnbd.pmf(params, 56, 0:9)), 10)
  expect_equal(length(bgcnbd.pmf(params, 56, 0)), 1)
  expect_equal(sum(bgcnbd.pmf(params, 2, 0:100)), 1)

  # test k=6
  set.seed(1)
  n <- 1000
  params <- c(k = 6, r = 0.85, alpha = 1.45, a = 0.79, b = 2.42)
  cbs <- bgcnbd.GenerateData(n,
                             round(runif(n, 36, 96) / 12) * 12,
                             c(32, 64),
                             params = params)$cbs
  est <- bgcnbd.EstimateParameters(cbs[, c("x", "t.x", "T.cal", "litt")])
  expect_equal(est[["k"]], params[["k"]])
})
