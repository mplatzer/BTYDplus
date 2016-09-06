
context("misc")

test_that("misc", {
  
  params <- c(k=1, r=2, alpha=3, a=4, 6)
  printnames <- c("k", "r", "alpha", "a", "b")
  expect_error(dc.check.model.params.safe(printnames, params, 'foo'))

  params <- c(k=1, r=2, alpha=3, a=4, b=6)
  printnames <- c("k", "r", "alpha", "a", "b")
  expect_silent(dc.check.model.params.safe(printnames, params, 'foo'))

})  

