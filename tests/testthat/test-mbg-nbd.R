context("mle")

test_that("MBG/NBD", {
  cat('test MBG/NBD\n')
  
  # generate artificial MBG/NBD data
  set.seed(1)
  n <- 8000
  params <- c(r = 0.85, alpha = 1.45, a = 0.79, b = 2.42)
  data <- mbgnbd.GenerateData(n = n, T.cal = round(runif(n, 12, 96)/4)*4, T.star = 32, params = params, return.elog = TRUE)
  cbs <- data$cbs
  elog <- data$elog
  
  # estimate regularity & parameters
  est <- mbgnbd.EstimateParameters(cbs[, c("x", "t.x", "T.cal", "litt")])
  expect_equal(params, est, tolerance = 0.1)
  
  # estimate future transactions & P(alive) with true parameters
  cbs$x.est <- mbgnbd.ConditionalExpectedTransactions(params, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)
  cbs$palive <- mbgnbd.PAlive(params, cbs$x, cbs$t.x, cbs$T.cal)
  
  # require less than 5% deviation
  expect_equal(sum(cbs$x.star), sum(cbs$x.est), 0.01)
  expect_equal(sum(cbs$palive), sum(cbs$alive), 0.01)
  
  # test with CDNow
  cdnow <- cdnow.sample()
  cbs <- cdnow$cbs
  elog <- cdnow$elog
  params <- mbgnbd.EstimateParameters(cbs)
  
  expect_equal(dim(mbgnbd.pmf(params, c(28, 56), 0:9)), c(10, 2))
  expect_equal(length(mbgnbd.pmf(params, 56, 0:9)), 10)
  expect_equal(length(mbgnbd.pmf(params, 56, 0)), 1)
  expect_equal(sum(mbgnbd.pmf(params, 2, 0:100)), 1)
  expect_silent(mbgnbd.Expectation(params, 11))
  cum <- mbgnbd.ExpectedCumulativeTransactions(params, 11, 39, 12)
  expect_true(all(diff(cum) > 0))
  expect_equal(length(cum), 12)
  #FIXMEexpect_silent(mbgnbd.PlotTrackingInc(params, cbs$T.cal, max(cbs$T.cal+cbs$T.star), elog2inc(elog)))
  #FIXMEexpect_silent(mbgnbd.PlotTrackingCum(params, cbs$T.cal, max(cbs$T.cal+cbs$T.star), elog2cum(elog)))
  #FIXMEmat <- mbgnbd.PlotFrequencyInCalibration(params, cbs, 7)
  #expect_equal(mat[1,], mat[2,], tolerance = 0.1)
}) 
