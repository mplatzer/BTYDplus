context("mcmc")

test_that("Pareto/NBD Abe MCMC", {
  cat('test Pareto/NBD Abe\n')
  skip_on_cran()

  # generate artificial Pareto/NBD Abe with 2 covariates
  set.seed(1)
  params <- list()
  params$beta <- matrix(c(0.18, -2.5, 0.5, -0.3, -0.2, 0.8), byrow = TRUE, ncol = 2)
  params$gamma <- matrix(c(0.05, 0.1, 0.1, 0.2), ncol = 2)
  expect_silent(abe.GenerateData(n = 100, T.cal = 32, T.star = c(16, 32), params, return.elog = TRUE))
  n <- 5000
  cbs <- abe.GenerateData(n, T.cal = 32, T.star = 32, params, return.elog = FALSE)$cbs
  
  # estimate parameters
  draws <- abe.mcmc.DrawParameters(cbs, covariates = c("covariate_1", "covariate_2"))
  
  expect_true(all(c("level_1", "level_2") %in% names(draws)))
  expect_true(all(c("lambda", "mu", "z", "tau") %in% colnames(as.matrix(draws$level_1[[1]]))))
  expect_true(all(c("log_lambda_intercept", "log_mu_intercept", "log_lambda_covariate_1", "log_mu_covariate_1", 
    "log_lambda_covariate_2", "log_mu_covariate_2", "var_log_lambda", "cov_log_lambda_log_mu", "var_log_mu") %in% 
    colnames(as.matrix(draws$level_2))))
  expect_equal(length(draws$level_1), n)
  expect_true(is.mcmc.list(draws$level_1[[1]]))
  expect_true(is.mcmc.list(draws$level_2))
  
  est <- round(summary(draws$level_2)$quantiles[, "50%"], 3)
  # require less than 20% deviation in estimated location parameter beta
  est.beta <- matrix(est[1:6], ncol = 2, byrow = T)
  expect_lt(max(abs((est.beta - params$beta)/params$beta)), 0.2)
  # variance parameter gamma is difficult to identify, particularly for mu; hence we relax our checks
  expect_lt((est["var_log_lambda"] - params$gamma[1, 1])/params$gamma[1, 1], 0.3)
  expect_lt(abs(est["cov_log_lambda_log_mu"] - params$gamma[1, 2]), 0.05)
  # expect_lt((est['var_log_mu'] - params$gamma[2,2]) / params$gamma[2,2], 0.30)
  
  # estimate future transactions & P(alive)
  xstar <- mcmc.DrawFutureTransactions(cbs, draws, T.star = cbs$T.star)
  cbs$x.est <- apply(xstar, 2, mean)
  cbs$pactive <- apply(xstar, 2, function(x) mean(x > 0))
  cbs$palive <- mcmc.PAlive(draws)
  
  # require less than 10% deviation in aggregated future transactions
  ape <- function(act, est) abs(act - est)/act
  expect_lt(ape(sum(cbs$x.star), sum(cbs$x.est)), 0.1)
  expect_lt(ape(sum(cbs$palive), sum(cbs$alive)), 0.1)
  expect_lt(ape(sum(cbs$x.star > 0), sum(cbs$pactive)), 0.1)
  
  expect_true(min(cbs$x.star) >= 0)
  expect_true(all(cbs$x.star == round(cbs$x.star)))
  expect_true(all(cbs$palive >= 0 & cbs$palive <= 1))
  
}) 
