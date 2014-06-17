context("mle")

test_that("BG/NBD", {
  
  # generate artificial BG/NBD data 
  set.seed(1)
  params <- c(r=0.85, alpha=4.45, a=0.79, b=2.42)
  n <- 8000
  cbs <- bgnbd.GenerateData(n, runif(n, 12, 96), 32, params)$cbs
  
  # estimate parameters, and compare to true parameters
  est <- bgnbd.EstimateParameters(cbs[, c("x", "t.x", "T.cal")])
  
  # require less than 5% deviation in estimated parameters
  ape <- function(act, est) abs(act-est)/act
  expect_true(ape(params[1], est[1]) < 0.05)
  expect_true(ape(params[2], est[2]) < 0.05)
  expect_true(ape(params[3], est[3]) < 0.05)
  expect_true(ape(params[4], est[4]) < 0.05)
  
  # estimate future transactions & P(alive) with true parameters
  cbs$x.est <- bgnbd.ConditionalExpectedTransactions(params, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)
  cbs$palive <- bgnbd.PAlive(params, cbs$x, cbs$t.x, cbs$T.cal)

  # require less than 5% deviation
  expect_true(ape(sum(cbs$x.star), sum(cbs$x.est)) < 0.05)
  expect_true(ape(sum(cbs$palive), sum(cbs$alive)) < 0.05)

  expect_true(min(cbs$x.star)>=0)
  expect_true(all(cbs$x.star==round(cbs$x.star)))
  expect_true(all(cbs$palive>=0 & cbs$palive<=1))
  
}
