context("mle")

test_that("Gamma/Gompertz/NBD", {
  
  # generate artificial Gamma/Gompertz/NBD data 
  set.seed(1)
  params <- c(r=0.85, alpha=4.45, b=0.12, s=0.75, beta=9.55)
  cbs <- ggnbd.GenerateData(5000, 32, 32, params)$cbs
  
  # estimate parameters, and compare to true parameters
  est <- ggnbd.EstimateParameters(cbs[, c("x", "t.x", "T.cal")], trace=10)
  
  # require less than 5%-15% deviation in estimated parameters
  ape <- function(act, est) abs(act-est)/act
  expect_true(ape(params[1], est[1]) < 0.05)
  expect_true(ape(params[2], est[2]) < 0.05)
  expect_true(ape(params[3], est[3]) < 0.15) # drop-out parameters are hard to identify
  expect_true(ape(params[4], est[4]) < 0.15)
  expect_true(ape(params[5], est[5]) < 0.15)
  
  # estimate future transactions and P(alive) with true parameters
  cbs$x.est <- ggnbd.ConditionalExpectedTransactions(params, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)
  cbs$palive <- ggnbd.PAlive(params, cbs$x, cbs$t.x, cbs$T.cal)
  # require less than 5% deviation
  expect_true(ape(sum(cbs$x.star), sum(cbs$x.est)) < 0.05) # ERROR: ggnbd.ConditionalExpectedTransactions is far too high!?
  expect_true(ape(sum(cbs$palive), sum(cbs$alive)) < 0.05)

  expect_true(min(cbs$x.star)>=0)
  expect_true(all(cbs$x.star==round(cbs$x.star)))
  expect_true(all(cbs$palive>=0 & cbs$palive<=1))

})
