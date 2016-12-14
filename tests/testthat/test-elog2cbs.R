context("elog2cbs")

test_that("elog2cbs", {

  elog <- data.frame(cust = c(1, 1, 1, 1, 1, 2, 3), date = Sys.Date() + c(0, 14, 14, 28, 35, 7, 24))
  elog_time <- data.frame(cust = c(1, 1, 1, 1, 1, 2, 3), date = Sys.time() + c(0, 14, 14, 28, 35, 7, 24))
  elog_dt <- as.data.table(elog)
  elog_dt[, first := min(date), by = "cust"]
  elog_s <- copy(elog_dt)
  elog_s[, sales := .I]
  T.cal <- Sys.Date() + 21
  T.tot <- Sys.Date() + 30

  # check that the return type matches the input type
  expect_is(elog2cbs(elog), "data.frame")
  expect_is(elog2cbs(elog_dt), "data.table")

  # check that we can also properly handle POSIXt
  expect_equal(elog2cbs(elog_time, units = "secs", T.cal = as.character(min(elog_time$date) + 21))[, 1:4],
               elog2cbs(elog, units = "days", T.cal = as.character(Sys.Date() + 21))[, 1:4])

  # check column names
  expect_named(elog2cbs(elog),
               c("cust", "x", "t.x", "litt", "first", "T.cal"))
  expect_named(elog2cbs(elog_dt),
               c("cust", "x", "t.x", "litt", "first", "T.cal"))
  expect_named(elog2cbs(elog_s),
               c("cust", "x", "t.x", "litt", "sales", "first", "T.cal"))
  expect_named(elog2cbs(elog_s, T.cal = T.cal),
               c("cust", "x", "t.x", "litt", "sales", "first", "T.cal", "T.star", "x.star", "sales.star"))
  expect_named(elog2cbs(elog, T.cal = T.cal, T.tot = T.tot),
               c("cust", "x", "t.x", "litt", "first", "T.cal", "T.star", "x.star"))
  expect_named(elog2cbs(elog[, c("cust", "date")], T.cal = T.cal),
               c("cust", "x", "t.x", "litt", "first", "T.cal", "T.star", "x.star"))

  # check number of returned customers
  expect_equal(uniqueN(elog$cust), nrow(elog2cbs(elog)))
  expect_equal(uniqueN(elog[elog$date <= T.cal, "cust"]), nrow(elog2cbs(elog, T.cal = T.cal)))

  # check total number of transactions
  cbs <- elog2cbs(elog_dt)
  expect_equal(uniqueN(elog_dt[, .(cust, date)]), sum(cbs$x) + sum(cbs$x.star) + nrow(cbs))
  expect_true(all(cbs$T.star == 0))
  expect_true(all(cbs$x.star == 0))
  expect_true(all(cbs$sales.star == 0))
  cbs <- elog2cbs(elog_dt, T.cal = T.cal, T.tot = T.tot)
  expect_equal(uniqueN(elog_dt[elog_dt$first <= T.cal & date <= T.tot, .(cust, date)]),
               sum(cbs$x) + sum(cbs$x.star) + nrow(cbs))

  # check sales
  cbs <- elog2cbs(elog_s, T.cal = T.cal, T.tot = T.tot)
  expect_equal(elog_s[first <= T.cal & date <= T.cal, sum(sales)], sum(cbs$sales))
  expect_equal(elog_s[first <= T.cal & date <= T.tot, sum(sales)], sum(cbs$sales) + sum(cbs$sales.star))

  # check that T.cal and T.tot also accept characters
  cbs <- elog2cbs(elog_s, T.cal = as.character(T.cal), T.tot = as.character(T.tot))

  # check that T.tot can be in the future
  cbs <- elog2cbs(elog_s, T.tot = max(elog_s$date) + 70)
  expect_equal(unique(cbs$T.star), 70 / 7)

  # check for `units`
  expect_equal(elog2cbs(elog, units = "days")$T.cal / 7, elog2cbs(elog)$T.cal)
  expect_equal(elog2cbs(elog, units = "hours")$t.x / (7 * 24), elog2cbs(elog)$t.x)
  elog_time <- data.frame(cust = c(1, 1, 1, 1, 1, 2, 3), date = Sys.time() + c(0, 14, 14, 28, 35, 7, 24))
  expect_equal(elog2cbs(elog, units = "days")[, c(1:4, 6)], elog2cbs(elog_time, units = "secs")[, c(1:4, 6)])

  # check for errors
  expect_error(elog2cbs(), "elog")
  expect_error(elog2cbs(data.frame(cust = 1:3)), "cust")
  expect_error(elog2cbs(data.frame(date = Sys.Date())), "date")
})
