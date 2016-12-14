context("mle")

test_that("NBD", {

  # generate artificial NBD data
  set.seed(1)
  params <- c(r = 0.85, alpha = 4.45)
  expect_silent(nbd.GenerateData(100, 32, c(16, 32), params, "2010-01-01"))
  cbs <- nbd.GenerateData(1000, 32, 32, params)$cbs

  # estimate parameters, and compare to true parameters
  est <- nbd.EstimateParameters(cbs[, c("x", "T.cal")])

  # require less than 5% deviation in estimated parameters
  expect_equal(params, est, tolerance = 0.05)

  # estimate future transactions in holdout-period with true params
  cbs$x.est <- nbd.ConditionalExpectedTransactions(params, cbs$T.star, cbs$x, cbs$T.cal)

  # require less than 5% deviation in estimated transactions
  expect_equal(sum(cbs$x.star), sum(cbs$x.est), tolerance = 0.05)

  expect_true(min(cbs$x.star) >= 0)
  expect_true(all(cbs$x.star == round(cbs$x.star)))

})
