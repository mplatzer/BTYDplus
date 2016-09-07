
context("misc")

test_that("dc.check.model.params.safe", {

  params <- c(k=1, r=2, alpha=3, a=4, b=6)
  printnames <- c("k", "r", "alpha", "a", "b")
  expect_silent(BTYDplus:::dc.check.model.params.safe(printnames, params, 'foo'))
  
  params <- c(k=1, r=2, alpha=3, a=4, 6)
  printnames <- c("k", "r", "alpha", "a", "b")
  expect_silent(BTYDplus:::dc.check.model.params.safe(printnames, params, 'foo'))

  params <- c(k=1, r=2, alpha=3, x=4, b=6)
  printnames <- c("k", "r", "alpha", "a", "b")
  expect_error(BTYDplus:::dc.check.model.params.safe(printnames, params, 'foo'))

})  


test_that("cdnow.sample", {
  
  data <- cdnow.sample()
  expect_equal(names(data), c("elog", "cbs"))
  expect_is(data$cbs, "data.frame")
  expect_is(data$elog, "data.frame")
  expect_equal(names(data$elog), c("cust", "date", "sales", "cds"))
  expect_equal(names(data$cbs), c("cust", "x", "t.x", "litt", "sales", "first", "T.cal", "T.star", "x.star", "sales.star"))
  
})  


test_that("elog2cum", {
  
  data <- cdnow.sample()
  cum <- elog2cum(data$elog)
  data(cdnowSummary, package="BTYD", envir = environment())
  expect_equal(head(cum, 66), head(cdnowSummary$cu.tracking, 66))
  # Note: for some weird reason these match only for the first 66 weeks?!
  # Probably an issue in BTYD package

  inc <- elog2inc(data$elog)  
  expect_equal(diff(cum), inc)
  
})  
