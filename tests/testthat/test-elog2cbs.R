context("elog2cbs")

test_that("elog2cbs", {
  cat('test elog2cbs')
  
  elog <- data.frame(cust = c(1, 1, 1, 1, 1, 2, 3), date = Sys.Date() + c(0, 14, 14, 28, 35, 7, 24))
  elog_dt <- as.data.table(elog)
  elog_dt[, first := min(date), by="cust"]
  elog_s <- copy(elog_dt)
  elog_s[, sales := .I]
  T.cal <- Sys.Date() + 21
  T.tot <- Sys.Date() + 30
  
  # check that the return type matches the input type
  expect_is(elog2cbs(elog), "data.frame")
  expect_is(elog2cbs(elog_dt), "data.table")

  # check column names
  cols <- c("cust", "x", "t.x", "litt", "sales", "first", "T.cal", "T.star", "x.star", "sales.star")
  expect_named(elog2cbs(elog), cols)
  expect_named(elog2cbs(elog, T.cal = T.cal), cols)
  expect_named(elog2cbs(elog, T.cal = T.cal, T.tot = T.tot), cols)
  expect_named(elog2cbs(elog_s), cols)
  expect_named(elog2cbs(elog_dt), cols)
  
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
  
  # check for errors
  expect_error(elog2cbs(), "elog")
  expect_error(elog2cbs(data.frame(cust = 1:3)), "cust")
  expect_error(elog2cbs(data.frame(date = Sys.Date())) , "date")
})
