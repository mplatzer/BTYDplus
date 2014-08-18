context("mle")

test_that("Pareto/NBD Abe MCMC", {
  
  # generate artificial Pareto/NBD Abe
  set.seed(1)
  params <- list()
  params$beta <- c(log(1.2), log(0.08))
  params$gamma <- matrix(c(0.05, 0, 0, 0.2), ncol=2)
  n <- 1000
  cbs <- abe.GenerateData(n, T.cal=32, T.star=32, params)$cbs
  
  # estimate parameters
  draws <- abe.mcmc.DrawParameters(cbs, mcmc=5000, burnin=5000, thin=100)
  est <- as.list(summary(draws$level_2)$quantiles[, "50%"])
  
  expect_true(all(c("level_1", "level_2") %in% names(draws)))
  expect_equal(length(draws$level_1), n)
  expect_true(is.mcmc.list(draws$level_1[[1]]))
  expect_true(is.mcmc.list(draws$level_2))
  
  # require less than 10% deviation in estimated parameters
  ape <- function(act, est) abs(act-est)/act
  expect_less_than(ape(params$beta[1], est$log_lambda), 0.10)
  expect_less_than(ape(params$beta[2], est$log_mu), 0.10)
  expect_less_than(ape(params$gamma[1,1], est$sigma_lambda_lambda), 0.50)
  expect_less_than(abs(est$sigma_lambda_mu), 0.1)
  ## There are difficulties to estimate heterogeneity in mu accurately
  #expect_less_than(ape(params$gamma[2,2], est$sigma_mu_mu), 0.20) 
  
  # estimate future transactions & P(alive)
  xstar <- abe.mcmc.DrawFutureTransactions(cbs, draws, T.star=cbs$T.star)
  cbs$x.est <- apply(xstar, 2, mean)
  cbs$pactive <- apply(xstar, 2, function(x) mean(x>0))
  cbs$palive <- abe.mcmc.PAlive(cbs, draws)
  
  # require less than 5% deviation
  expect_less_than(ape(sum(cbs$x.star), sum(cbs$x.est)), 0.05)
  expect_less_than(ape(sum(cbs$palive), sum(cbs$alive)), 0.05)
  expect_less_than(ape(sum(cbs$x.star>0), sum(cbs$pactive)), 0.05)
  
  expect_true(min(cbs$x.star)>=0)
  expect_true(all(cbs$x.star==round(cbs$x.star)))
  expect_true(all(cbs$palive>=0 & cbs$palive<=1))
  
}
