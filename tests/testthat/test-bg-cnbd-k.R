context("mle")

test_that("BG/CNBD-k", {
  cat('test BG/CNBD-k\n')
  
  # validate against BTYD implementation
  set.seed(1)
  params <- c(0.85, 1.45, 0.79, 2.42)
  n <- 8000
  data <- bgnbd.GenerateData(n = n, T.cal = round(runif(n, 4, 32)), T.star = 32, params = params, return.elog = TRUE)
  cbs  <- data$cbs
  elog <- data$elog
  params.est.btyd <- BTYD::bgnbd.EstimateParameters(cbs)
  params.est.btyd_plus <- bgcnbd.EstimateParameters(cbs, k = 1)[-1]
  expect_equal(unname(round(params.est.btyd, 2)), 
               unname(round(params.est.btyd_plus, 2)))
  expect_equal(BTYD::bgnbd.PAlive(params, 0, 0, 32), 
               bgcnbd.PAlive(c(1, params), 0, 0, 32))
  expect_equal(BTYD::bgnbd.PAlive(params, 1, 16, 32), 
               bgcnbd.PAlive(c(1, params), 1, 16, 32))
  expect_equal(BTYD::bgnbd.Expectation(params, 1:3), 
               bgcnbd.Expectation(c(1, params), 1:3))
  expect_equal(BTYD::bgnbd.ConditionalExpectedTransactions(params, 32, 1, 16, 32), 
               bgcnbd.ConditionalExpectedTransactions(c(1, params), 32, 1, 16, 32))
  expect_equal(unname(BTYD::bgnbd.pmf(params, 32, 0:2)), 
               unname(bgcnbd.pmf(c(1, params), 32, 0:2)))
  expect_equal(BTYD::bgnbd.PlotFrequencyInCalibration(params, cbs, 7), 
               bgcnbd.PlotFrequencyInCalibration(c(1, params), cbs, 7))
  inc.tracking <- elog2inc(elog, by = 1)
  expect_equal(BTYD::bgnbd.PlotTrackingInc(params, cbs$T.cal, max(cbs$T.cal) + 32, inc.tracking), 
               bgcnbd.PlotTrackingInc(c(1, params), cbs$T.cal, max(cbs$T.cal) + 32, inc.tracking), 
               tolerance = 0.001)
  expect_equal(BTYD::bgnbd.ExpectedCumulativeTransactions(params, 11, 39, 12),
               bgcnbd.ExpectedCumulativeTransactions(c(1, params), 11, 39, 12),
               tolerance = 0.001)
  cu.tracking <- cumsum(inc.tracking)
  expect_equal(BTYD::bgnbd.PlotTrackingCum(params, cbs$T.cal, 32 + 32, cu.tracking), 
               bgcnbd.PlotTrackingCum(c(1, params), cbs$T.cal, 32 + 32, cu.tracking), 
               tolerance = 0.001)
  
  # generate artificial BG/CNBD-k data
  set.seed(1)
  n <- 8000
  params <- c(k = 3, r = 0.85, alpha = 1.45, a = 0.79, b = 2.42)
  data <- bgcnbd.GenerateData(n = n, T.cal = round(runif(n, 12, 96)/4)*4, T.star = c(32, 64), params = params, return.elog = TRUE)
  cbs <- data$cbs
  elog <- data$elog
  
  # estimate regularity & parameters
  k.est <- estimateRegularity(elog)
  est <- bgcnbd.EstimateParameters(cbs[, c("x", "t.x", "T.cal", "litt")])
  est.fixed.k <- bgcnbd.EstimateParameters(cbs[, c("x", "t.x", "T.cal", "litt")], k = params[1])
  
  expect_equal(params[1], est[1])
  expect_equal(est, est.fixed.k)
  
  # require less than 10% deviation in estimated parameters
  ape <- function(act, est) abs(act - est)/act
  expect_true(ape(params[1], k.est) < 0.05)
  expect_true(ape(params[2], est[2]) < 0.1)
  expect_true(ape(params[3], est[3]) < 0.1)
  expect_true(ape(params[4], est[4]) < 0.1)
  expect_true(ape(params[5], est[5]) < 0.1)
  
  # estimate future transactions & P(alive) with true parameters
  cbs$x.est32 <- bgcnbd.ConditionalExpectedTransactions(params, 32, cbs$x, cbs$t.x, cbs$T.cal)
  cbs$x.est64 <- bgcnbd.ConditionalExpectedTransactions(params, 64, cbs$x, cbs$t.x, cbs$T.cal)
  cbs$palive <- bgcnbd.PAlive(params, cbs$x, cbs$t.x, cbs$T.cal)
  
  # require less than 5% deviation
  expect_true(ape(sum(cbs$x.star32), sum(cbs$x.est32)) < 0.05)
  expect_true(ape(sum(cbs$x.star64), sum(cbs$x.est64)) < 0.05)
  expect_true(ape(sum(cbs$palive), sum(cbs$alive)) < 0.05)
  
  # some basic sanity checks
  expect_true(min(cbs$x.star32) >= 0)
  expect_true(all(cbs$x.star32 == round(cbs$x.star32)))
  expect_true(all(cbs$palive >= 0 & cbs$palive <= 1))
  
}) 
