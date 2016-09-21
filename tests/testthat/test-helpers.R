
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
  expect_equal(nrow(data$cbs), 2357)
  expect_equal(min(data$cbs$T.cal), 27)
})


test_that("elog2cum", {

  data <- cdnow.sample()
  cum <- elog2cum(data$elog)
  utils::data(cdnowSummary, package="BTYD", envir = environment())
  expect_equal(cum, cdnowSummary$cu.tracking)

  inc <- elog2inc(data$elog)
  expect_equal(diff(cum), inc)

  elog <- data.table(cust = c(1, 1, 1, 3, 3), t = c(0, 9, 18, 4, 6))
  expect_equal(elog2cum(elog, by = 7), c(1, 2, 3))
  expect_equal(elog2cum(elog, by = 7, first = TRUE), c(3, 4, 5))
  expect_equal(diff(elog2cum(elog, by = 7, first = TRUE)), elog2inc(elog, by = 7, first = TRUE))
})


test_that("plotSampledTimingPatterns", {

  elog <- cdnow.sample()$elog
  expect_silent(plotSampledTimingPatterns(elog))
  expect_silent(plotSampledTimingPatterns(elog, T.cal = "1997-09-30"))
  expect_silent(plotSampledTimingPatterns(elog, n = 100))
  expect_silent(plotSampledTimingPatterns(head(elog, 5), T.cal = "1997-09-30"))

})
