
context("misc")

test_that("dc.check.model.params.safe", {

  params <- c(k = 1, r = 2, alpha = 3, a = 4, b = 6)
  printnames <- c("k", "r", "alpha", "a", "b")
  expect_silent(BTYDplus:::dc.check.model.params.safe(printnames, params, "foo"))

  params <- c(k = 1, r = 2, alpha = 3, a = 4, 6)
  printnames <- c("k", "r", "alpha", "a", "b")
  expect_silent(BTYDplus:::dc.check.model.params.safe(printnames, params, "foo"))

  params <- c(k = 1, r = 2, alpha = 3, x = 4, b = 6)
  printnames <- c("k", "r", "alpha", "a", "b")
  expect_error(BTYDplus:::dc.check.model.params.safe(printnames, params, "foo"))

})


test_that("elog2cum", {

  cdnowElog <- read.csv(system.file("data/cdnowElog.csv", package = "BTYD"),
                        stringsAsFactors = FALSE,
                        col.names = c("cust", "sampleid", "date", "cds", "sales"))
  cdnowElog$date <- as.Date(as.character(cdnowElog$date), format = "%Y%m%d")
  cum <- elog2cum(cdnowElog)
  utils::data(cdnowSummary, package = "BTYD", envir = environment())
  expect_equal(cum, cdnowSummary$cu.tracking)

  inc <- elog2inc(cdnowElog)
  expect_equal(diff(cum), inc)

  elog <- data.frame(cust = c(1, 1, 1, 3, 3), t = c(0, 9, 20, 4, 6))
  expect_equal(elog2cum(elog, by = 7), c(1, 2, 3))
  expect_equal(elog2cum(elog, by = 7, first = TRUE), c(3, 4, 5))
  expect_equal(diff(elog2cum(elog, by = 7, first = TRUE)), elog2inc(elog, by = 7, first = TRUE))
  expect_equal(tail(elog2cum(elog, by = 1, first = TRUE), 1), nrow(elog))
  expect_equal(tail(elog2cum(elog, by = 1, first = FALSE), 1), nrow(elog) - uniqueN(elog$cust))

  elog <- data.table(cust = c(1, 1, 1, 1, 1, 3, 3), t = c(0, 9, 9, 20, 22, 4, 6))
  expect_equal(elog2cum(elog, by = 7), c(1, 2, 3))
})


test_that("plotTimingPatterns", {

  data("groceryElog")
  expect_silent(plotTimingPatterns(groceryElog))
  expect_silent(plotTimingPatterns(groceryElog, headers = c("X", "Y")))
  expect_silent(plotTimingPatterns(groceryElog, headers = c("X", "Y"), title = ""))
  expect_silent(plotTimingPatterns(groceryElog, T.cal = "2006-12-31", title = NULL))
  expect_silent(plotTimingPatterns(groceryElog, T.cal = "2006-12-31", T.tot = "2007-06-30",
                                   headers = c("Past", "Future")))
  expect_silent(plotTimingPatterns(groceryElog, T.cal = as.Date("2006-12-31")))
  expect_silent(plotTimingPatterns(groceryElog, n = 100))
  expect_silent(plotTimingPatterns(head(groceryElog, 10), T.cal = "2006-12-31", T.tot = "2007-12-30"))

})
