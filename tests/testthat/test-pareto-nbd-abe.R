context("mcmc")

test_that("Pareto/NBD (Abe) MCMC", {

  # test basic data simulation
  n <- 100
  params <- list()
  params$beta <- matrix(c(0.18, -2.5, 0.5, -0.3, -0.2, 0.8), byrow = TRUE, ncol = 2)
  params$gamma <- matrix(c(0.05, 0.1, 0.1, 0.2), ncol = 2)
  expect_silent(abe.GenerateData(n, 36, 36, params,
                                 covariates = matrix(c(rnorm(n), runif(n)), ncol = 2)))
  expect_silent(abe.GenerateData(n, 36, 36, params,
                                 covariates = matrix(c(rnorm(n), runif(n)), ncol = 2,
                                                     dimnames = list(NULL, c("x1", "x2")))))
  expect_error(abe.GenerateData(n, 36, 36, params,
                                covariates = matrix(c(rnorm(n)), ncol = 1,
                                                    dimnames = list(NULL, c("x1")))),
               "covariate columns")
  expect_error(abe.GenerateData(n, 36, 36, params,
                                covariates = matrix(c(rnorm(n), runif(n), runif(n)), ncol = 3,
                                                    dimnames = list(NULL, c("x1", "x2", "x3")))),
               "covariate columns")
  params$beta <- params$beta[1:2, ]
  expect_silent(abe.GenerateData(n, 36, 36, params, covariates = rnorm(n)))
  expect_silent(abe.GenerateData(n, 36, 36, params, covariates = matrix(rnorm(n), ncol = 1,
                                                                  dimnames = list(NULL, "x1"))))

  # generate artificial Pareto/NBD (Abe) with 2 covariates
  set.seed(1)
  params <- list()
  params$beta <- matrix(c(0.18, -2.5, 0.5, -0.3, -0.2, 0.8), byrow = TRUE, ncol = 2)
  params$gamma <- matrix(c(0.05, 0.1, 0.1, 0.2), ncol = 2)
  n <- 5000
  sim <- abe.GenerateData(n,
                          round(runif(n, 36, 96) / 12) * 12,
                          36,
                          params,
                          "2010-01-01")
  # TODO: test param recoverability with covars being provided to abe.GenerateData
  # covars <- matrix(c(runif(n, 0, 10), runif(n, 10, 20)), ncol = 2,
  #                  dimnames = list(NULL, c("x1", "x2")))

  cbs <- sim$cbs
  expect_true(all(c("intercept", "covariate_1", "covariate_2") %in% names(cbs)))

  # test basic parameter estimation
  draws <- abe.mcmc.DrawParameters(as.data.table(cbs),
                                   mcmc = 10, burnin = 0, thin = 1, mc.cores = 1)
  draws <- abe.mcmc.DrawParameters(as.data.table(cbs), covariates = c("covariate_1"),
                                   mcmc = 10, burnin = 0, thin = 1, mc.cores = 1)

  # test parameter recovery
  skip("skip long-running test of MCMC parameter recovery")
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

  # check tracking plots
  inc <- elog2inc(sim$elog)
  mat <- mcmc.PlotTrackingInc(draws,
                              T.cal = cbs$T.cal,
                              T.tot = max(cbs$T.cal + cbs$T.star),
                              actual.inc.tracking.data = inc,
                              covariates = cbs[, c("covariate_1", "covariate_2")])
  mat.sum <- apply(mat, 1, sum)
  expect_equal(mat[1], mat[2], tolerance = 0.1)
  cum <- elog2cum(sim$elog)
  mat <- mcmc.PlotTrackingCum(draws,
                              T.cal = cbs$T.cal,
                              T.tot = max(cbs$T.cal + cbs$T.star),
                              actual.cu.tracking.data = cum,
                              covariates = cbs[, c("covariate_1", "covariate_2")])
  expect_equal(mat[, ncol(mat) - 1], mat[, ncol(mat) - 1], tolerance = 0.1)

})
