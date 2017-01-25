context("timing")

test_that("estimateRegularity", {

  # generate Erlang-3 with various lambdas
  set.seed(1)
  nr_of_customers <- 1000
  k <- 3
  elog <- rbindlist(lapply(1:nr_of_customers, function(i) {
    lambda <- exp(rnorm(1))
    data.table(cust = i, t = cumsum(rgamma(50, k, k * lambda)))
  }))

  # Estimate regularity parameter k
  k_est_1 <- estimateRegularity(elog)
  k_est_1 <- estimateRegularity(elog, plot = TRUE, min = 2)
  expect_error(k_est_1 <- estimateRegularity(elog, plot = TRUE, min = 60), "sufficient")
  k_est_1 <- estimateRegularity(elog, plot = TRUE, title = "Plot Title")
  k_est_1 <- estimateRegularity(elog, plot = TRUE, title = "")
  k_est_2 <- estimateRegularity(elog, method = "mle", plot = TRUE)
  k_est_3 <- estimateRegularity(elog, method = "mle-minka", plot = TRUE, title = "Plot Title")
  k_est_4 <- estimateRegularity(elog, method = "mle-thom", plot = TRUE, title = "")
  k_est_5 <- estimateRegularity(elog, method = "cv", plot = TRUE)
  expect_equal(k, k_est_1, tolerance = 0.1)
  expect_equal(k, k_est_2, tolerance = 0.1)
  expect_equal(k, k_est_3, tolerance = 0.1)
  expect_equal(k, k_est_4, tolerance = 0.1)
  expect_equal(k, k_est_5, tolerance = 0.1)

})
