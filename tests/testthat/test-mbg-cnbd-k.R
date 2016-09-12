context("mle")

test_that("MBG/CNBD-k", {
  cat('test MBG/CNBD-k\n')
  
  # generate artificial MBG/CNBD-k data
  set.seed(1)
  n <- 8000
  params <- c(k = 3, r = 0.85, alpha = 1.45, a = 0.79, b = 2.42)
  data <- mbgcnbd.GenerateData(n = n, T.cal = round(runif(n, 12, 96)/4)*4, T.star = 32, params = params, return.elog = TRUE)
  cbs <- data$cbs
  elog <- data$elog
  
  # estimate regularity & parameters
  k.est <- estimateRegularity(elog)
  est <- mbgcnbd.EstimateParameters(cbs[, c("x", "t.x", "T.cal", "litt")])
  est.fixed.k <- mbgcnbd.EstimateParameters(cbs[, c("x", "t.x", "T.cal", "litt")], k = params[1])
  
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
  cbs$x.est <- mbgcnbd.ConditionalExpectedTransactions(params, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)
  cbs$palive <- mbgcnbd.PAlive(params, cbs$x, cbs$t.x, cbs$T.cal)
  
  # require less than 5% deviation
  expect_true(ape(sum(cbs$x.star), sum(cbs$x.est)) < 0.05)
  expect_true(ape(sum(cbs$palive), sum(cbs$alive)) < 0.05)
  
  expect_true(min(cbs$x.star) >= 0)
  expect_true(all(cbs$x.star == round(cbs$x.star)))
  expect_true(all(cbs$palive >= 0 & cbs$palive <= 1))
  
  # test with CDNow
  cdnow <- cdnow.sample()
  cbs <- cdnow$cbs
  elog <- cdnow$elog
  params <- mbgcnbd.EstimateParameters(cbs)
  
  expect_equal(dim(mbgcnbd.pmf(params, c(28, 56), 0:9)), c(10, 2))
  expect_equal(length(mbgcnbd.pmf(params, 56, 0:9)), 10)
  expect_equal(length(mbgcnbd.pmf(params, 56, 0)), 1)
  expect_equal(sum(mbgcnbd.pmf(params, 2, 0:100)), 1)
  expect_silent(exp <- mbgcnbd.Expectation(params, 11))
  cum <- mbgcnbd.ExpectedCumulativeTransactions(params, 11, 39, 12)
  expect_true(all(diff(cum) > 0))
  expect_equal(length(cum), 12)
  expect_silent(mbgcnbd.PlotTrackingInc(params, cbs$T.cal, max(cbs$T.cal+cbs$T.star), elog2inc(elog)))
  expect_silent(mbgcnbd.PlotTrackingCum(params, cbs$T.cal, max(cbs$T.cal+cbs$T.star), elog2cum(elog)))
  mat <- mbgcnbd.PlotFrequencyInCalibration(params, cbs, 7)
  expect_equal(mat[1,], mat[2,], tolerance = 0.1)
}) 
