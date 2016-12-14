context("mcmc")

test_that("Pareto/NBD (Abe) MCMC", {
  skip_on_cran()

  # generate artificial Pareto/NBD (Abe) with 2 covariates
  set.seed(1)
  params <- list()
  params$beta <- matrix(c(0.18, -2.5, 0.5, -0.3, -0.2, 0.8), byrow = TRUE, ncol = 2)
  params$gamma <- matrix(c(0.05, 0.1, 0.1, 0.2), ncol = 2)
  expect_silent(abe.GenerateData(n = 100, T.cal = 32, T.star = c(16, 32), params))
  n <- 5000
  sim <- abe.GenerateData(n,
                          round(runif(n, 36, 96) / 12) * 12,
                          36,
                          params,
                          "2010-01-01")
  cbs <- sim$cbs

  # test basic parameter estimation
  draws <- abe.mcmc.DrawParameters(as.data.table(cbs), covariates = c("covariate_1"),
                                   mcmc = 10, burnin = 0, thin = 1, mc.cores = 1)

  # test parameter recovery
  draws <- abe.mcmc.DrawParameters(cbs, covariates = c("covariate_1", "covariate_2"), mc.cores = 1)

  expect_true(all(c("level_1", "level_2") %in% names(draws)))
  expect_true(all(c("lambda", "mu", "z", "tau") %in% colnames(as.matrix(draws$level_1[[1]]))))
  expect_true(all(c("log_lambda_intercept", "log_mu_intercept", "log_lambda_covariate_1", "log_mu_covariate_1",
    "log_lambda_covariate_2", "log_mu_covariate_2", "var_log_lambda", "cov_log_lambda_log_mu", "var_log_mu") %in%
    colnames(as.matrix(draws$level_2))))
  expect_equal(length(draws$level_1), n)
  expect_true(coda::is.mcmc.list(draws$level_1[[1]]))
  expect_true(coda::is.mcmc.list(draws$level_2))

  est <- round(summary(draws$level_2)$quantiles[, "50%"], 3)
  # require less than 5% deviation in estimated location parameter beta
  expect_equal(matrix(est[1:6], ncol = 2, byrow = T), params$beta, tolerance = 0.05)
  # variance parameter gamma is difficult to identify, particularly for mu;
  # hence we increase tolerance to 10%
  expect_equal(unname(est["var_log_lambda"]), params$gamma[1, 1], tolerance = 0.1)
  expect_equal(unname(est["cov_log_lambda_log_mu"]), params$gamma[1, 2], tolerance = 0.1)
  expect_equal(unname(est["var_log_mu"]), params$gamma[2, 2], tolerance = 0.1)

  # estimate future transactions & P(alive)
  xstar <- mcmc.DrawFutureTransactions(cbs, draws, T.star = cbs$T.star)
  cbs$x.est <- apply(xstar, 2, mean)
  cbs$pactive <- apply(xstar, 2, function(x) mean(x > 0))
  cbs$palive <- mcmc.PAlive(draws)

  # require less than 10% deviation in aggregated future transactions
  expect_equal(sum(cbs$x.star), sum(cbs$x.est), tolerance = 0.1)
  expect_equal(sum(cbs$palive), sum(cbs$alive), tolerance = 0.1)
  expect_equal(sum(cbs$x.star > 0), sum(cbs$pactive), tolerance = 0.1)

  expect_true(min(cbs$x.star) >= 0)
  expect_true(all(cbs$x.star == round(cbs$x.star)))
  expect_true(all(cbs$palive >= 0 & cbs$palive <= 1))

})
