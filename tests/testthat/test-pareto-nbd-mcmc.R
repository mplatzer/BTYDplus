context("mcmc")

test_that("Pareto/NBD MCMC", {

  # test Pareto/NBD data generator
  params <- list(r = 0.9, alpha = 10, s = 0.8, beta = 12)
  n <- 5
  T.cal <- c(26, 26, 28.5, 52, 52)
  T.star <- 52
  date.zero <- "2010-01-01"
  sim <- pggg.GenerateData(5, T.cal, T.star, params, date.zero)
  expect_equal(names(sim), c("cbs", "elog"))
  expect_is(sim$cbs, "data.frame")
  expect_is(sim$elog, "data.frame")
  expect_equal(names(sim$elog), c("cust", "t", "date"))
  expect_equal(nrow(sim$cbs), n)
  expect_equal(uniqueN(sim$elog$cust), n)
  expect_is(sim$elog$date, "POSIXct")
  expect_equal(min(sim$elog$date), as.POSIXct(date.zero))
  expect_equal(min(sim$cbs$first), as.POSIXct(date.zero))
  simElog <- sim$elog
  simCBS <- sim$cbs
  # recreate CBS via `elog2cbs`
  date.zero <- min(simElog$date)
  simCBSc <- elog2cbs(simElog,
                      T.cal = date.zero + max(T.cal) * 3600 * 24 * 7,
                      T.tot = date.zero + max(T.cal + T.star) * 3600 * 24 * 7)
  expect_equal(simCBSc, subset(simCBS, select = names(simCBSc)))
  # multiple T.star's
  sim <- pggg.GenerateData(100, 52, c(26, 104), params, date.zero = as.Date("2010-01-01"))
  expect_true(all(c("x.star26", "x.star104") %in% names(sim$cbs)))
  expect_true(all(sim$cbs$x.star104 >= sim$cbs$x.star26))
  expect_true(any(sim$cbs$x.star104 > sim$cbs$x.star26))

  # add test case for issue with varying T.cal's;
  # see https://github.com/mplatzer/BTYDplus/issues/42
  n <- 5000
  params <- list(r = 0.9, alpha = 10, s = 0.8, beta = 12)
  # constant T.cal
  T.cal <- rep(52, n)
  cbs1 <- pnbd.GenerateData(n, T.cal, 52, params)$cbs
  # varying T.cal, with one T.cal being particularly large
  T.cal[n] <- 52
  cbs2 <- pnbd.GenerateData(n, T.cal, 52, params)$cbs
  expect_equal(sum(cbs1[-n]$x.star), sum(cbs2[-n]$x.star), tolerance = 0.2)


  skip_on_cran()

  # generate artificial Pareto/NBD data
  set.seed(1)
  params <- list(r = 0.9, alpha = 5, s = 0.8, beta = 12)
  n <- 5000
  sim <- pnbd.GenerateData(n,
                           round(runif(n, 36, 96) / 12) * 12,
                           36,
                           params)
  cbs <- sim$cbs

  # test basic parameter estimation
  draws <- pnbd.mcmc.DrawParameters(as.data.table(cbs),
                                    mcmc = 10, burnin = 0, thin = 1, mc.cores = 1,
                                    param_init = list(r = 1, alpha = 1, s = 1, beta = 1))

  # test parameter recovery
  draws <- pnbd.mcmc.DrawParameters(cbs, mc.cores = 1, trace = 1000)
  est <- as.list(summary(draws$level_2)$quantiles[, "50%"])

  expect_true(all(c("level_1", "level_2") %in% names(draws)))
  expect_equal(length(draws$level_1), n)
  expect_true(coda::is.mcmc.list(draws$level_1[[1]]))
  expect_true(coda::is.mcmc.list(draws$level_2))

  # require less than 10% deviation in estimated parameters
  expect_equal(params, est, tolerance = 0.10)

  # estimate future transactions & P(alive)
  xstar <- mcmc.DrawFutureTransactions(cbs, draws, T.star = cbs$T.star)
  cbs$x.est <- apply(xstar, 2, mean)
  cbs$pactive <- apply(xstar, 2, function(x) mean(x > 0))
  cbs$palive <- mcmc.PAlive(draws)

  # require less than 5% deviation
  expect_equal(sum(cbs$x.star), sum(cbs$x.est), tolerance = 0.05)
  expect_equal(sum(cbs$palive), sum(cbs$alive), tolerance = 0.05)
  expect_equal(sum(cbs$x.star > 0), sum(cbs$pactive), tolerance = 0.05)

  expect_true(min(cbs$x.star) >= 0)
  expect_true(all(cbs$x.star == round(cbs$x.star)))
  expect_true(all(cbs$palive >= 0 & cbs$palive <= 1))

  # estimate parameters via Ma/Liu
  set.seed(1)
  draws_maliu <- pnbd.mcmc.DrawParameters(cbs,
                                          mcmc = 10, burnin = 0, thin = 2, chains = 2,
                                          use_data_augmentation = FALSE, mc.cores = 1)
  expect_equal(apply(as.matrix(draws_maliu$level_2), 2, mean),
               apply(as.matrix(draws$level_2), 2, mean), tolerance = 0.2)
})
