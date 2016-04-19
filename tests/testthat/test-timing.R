context("timing")

test_that("estimateRegularity", {
  cat('test timing\n')
  
  # generate Erlang-3 with various lambdas
  set.seed(1)
  nr_of_customers <- 1000
  k <- 3
  elog <- rbindlist(lapply(1:nr_of_customers, function(i) {
    lambda <- exp(rnorm(1))
    data.table(cust = i, t = cumsum(rgamma(50, k, k * lambda)))
  }))
  
  # Estimate regularity parameter k
  k.est.1 <- estimateRegularity(elog)
  k.est.2 <- estimateRegularity(elog, method = "mle")
  k.est.3 <- estimateRegularity(elog, method = "mle-minka")
  k.est.4 <- estimateRegularity(elog, method = "mle-thom")
  k.est.5 <- estimateRegularity(elog, method = "cv")
  
  # require less than 5% deviation in estimated parameters
  ape <- function(act, est) abs(act - est)/act
  expect_true(ape(k, k.est.1) < 0.05)
  expect_true(ape(k, k.est.2) < 0.05)
  expect_true(ape(k, k.est.3) < 0.05)
  expect_true(ape(k, k.est.4) < 0.05)
  expect_true(ape(k, k.est.5) < 0.05)
  
}) 
