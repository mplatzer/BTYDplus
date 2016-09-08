context("mcmc")

test_that("MCMC Helpers", {
  
  # generate small Pareto/NBD cohort
  params <- list(k = 1, r = 0.9, alpha = 10, s = 0.8, beta = 12)
  cbs <- pggg.GenerateData(n = 100, 26, 26, params, FALSE)$cbs
  
  # estimate parameters
  mcmc <- 100
  burnin <- 10
  thin <- 10
  chains <- 2
  param.draws <- pnbd.mcmc.DrawParameters(cbs, mc.cores = 1, mcmc, burnin, thin, chains)

  # draw future transactions
  xstar.draws <- mcmc.DrawFutureTransactions(cbs, param.draws)
  expect_is(xstar.draws, "matrix")
  expect_equal(dim(xstar.draws), c(mcmc*chains/thin, nrow(cbs)))
  size <- 1
  xstar.draws <- mcmc.DrawFutureTransactions(cbs, param.draws, size = size)
  expect_equal(dim(xstar.draws), c(size, nrow(cbs)))
  size <- (mcmc*chains / thin) + 1
  xstar.draws <- mcmc.DrawFutureTransactions(cbs, param.draws, size = size)
  expect_equal(dim(xstar.draws), c(size, nrow(cbs)))
  expect_silent(xstar.draws <- mcmc.DrawFutureTransactions(cbs, param.draws, T.star = 10))
  
  # test setBurnin
  burnin2 <- 30
  xstar.draws2 <- mcmc.setBurnin(param.draws, burnin = burnin2)
  expect_equal(start(xstar.draws2$level_2), burnin2)
  expect_equal(start(xstar.draws2$level_1[[1]]), burnin2)
  xstar.draws <- mcmc.DrawFutureTransactions(cbs, xstar.draws2)
  expect_equal(dim(xstar.draws), c((mcmc+burnin-burnin2)*chains/thin, nrow(cbs)))

  # test Palive
  palive <- mcmc.PAlive(cbs, param.draws)
  expect_is(palive, "numeric")
  expect_true(all(palive >= 0))
  expect_true(all(palive <= 1))
  
  # test Pactive
  pactive <- mcmc.PActive(xstar.draws)
  expect_is(pactive, "numeric")
  expect_true(all(pactive >= 0))
  expect_true(all(pactive <= 1))
  
  # test plotPActiveDiagnostic
  expect_silent(mcmc.plotPActiveDiagnostic(cbs, xstar.draws))
  
})