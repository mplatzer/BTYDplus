
#' Faster implementation of BTYD::dc.ElogToCbsCbt that also returns summary statistic for estimating regularity
#'
#' Returns data.frame with
#'      cust:   customer id (unique key)
#'      x:      nr of recurring events in calibration period
#'      t.x:    time between first and last event in calibration period
#'      litt:   sum of logarithmic intertransaction timings durint calibration period 
#'              this is a summary statistic for estimating regularity
#'      sales:  sum of sales in calibration period
#'      first:  date of first transaction in calibration period
#'      T.cal:  time between first event and end of calibration period
#'      T.star: length of holdout period
#'      x.star: nr of events within holdout period
#'      sales.star: sum of sales within holdout period
#'      
#' Customers without any transaction during calibration period are being dropped from the result.
#' Transactions with identical `cust` and `date` field are treated as a single transaction, with `sales` being summed up
#'
#' @param elog data.frame with columns 'cust' and 'date'; optionally with column 'sales'
#' @param per time unit, either 'week', 'day', 'hour', 'min', 'sec'
#' @param T.cal end date of calibration period
#' @param T.tot end date of holdout period
#' @return data.frame
#' @export
elog2cbs <- function(elog, per = "week", T.cal = max(elog$date), T.tot = max(elog$date)) {
  cust <- first <- itt <- T.star <- x.star <- sales <- sales.star <- NULL  # avoid checkUsage warnings
  stopifnot(inherits(elog, "data.frame"))
  stopifnot(all(c("cust", "date") %in% names(elog)))
  stopifnot(any(c("Date", "POSIXt") %in% class(elog$date)))
  
  is.dt <- is.data.table(elog)
  has.sales <- "sales" %in% names(elog)
  # convert to data.table for improved performance
  elog_dt <- data.table(elog)
  setkey(elog_dt, cust, date)
  # check for `sales` column, and populate if missing
  if (!has.sales) {
    elog_dt[, `:=`(sales, 1)]
  } else {
    stopifnot(is.numeric(elog_dt$sales))
  }
  # merge transactions with same dates
  elog_dt <- elog_dt[, list(sales = sum(sales)), by = "cust,date"]
  # determine time since first date for each customer
  elog_dt[, `:=`(first, min(date)), by = "cust"]
  elog_dt[, `:=`(t, as.numeric(difftime(date, first, units = per))), by = "cust"]
  # compute intertransaction times
  elog_dt[, `:=`(itt, c(0, diff(t))), by = "cust"]
  # count events in calibration period
  cbs <- elog_dt[date <= T.cal, list(x = .N - 1, t.x = max(t), litt = sum(log(itt[itt > 0])), sales = sum(sales)), 
    by = "cust,first"]
  cbs[, `:=`(T.cal, as.numeric(difftime(T.cal, first, units = per)))]
  cbs[, `:=`(T.star, as.numeric(difftime(T.tot, first, units = per)) - T.cal)]
  setkey(cbs, cust)
  # count events in validation period
  val <- elog_dt[date > T.cal & date <= T.tot, list(x.star = .N, sales.star = sum(sales)), keyby = "cust"]
  cbs <- merge(cbs, val, all.x = TRUE, by = "cust")
  cbs[is.na(x.star), `:=`(x.star, 0)]
  cbs[is.na(sales.star), `:=`(sales.star, 0)]
  setcolorder(cbs, c("cust", "x", "t.x", "litt", "sales", "first", "T.cal", "T.star", "x.star", "sales.star"))
  # return same object type as was passed
  if (!has.sales) {
    elog_dt[, `:=`(sales, NULL)]
  }
  if (!is.dt) {
    cbs <- data.frame(cbs)
  }
  return(cbs)
}
