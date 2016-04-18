
#' Faster implementation of BTYD::dc.ElogToCbsCbt that also returns summary statistic for estimating regularity
#'
#' Returns data.frame with
#'      cust:   customer id (unique key)
#'      x:      nr of recurring events in calibration period
#'      t.x:    time between first and last event in calibration period
#'      litt:   sum of logarithmic intertransaction timings durint calibration period 
#'              this is a summary statistic for estimating regularity
#'      T.cal:  time between first event and end of calibration period
#'      T.star: length of holdout period
#'      x.star: nr of events within holdout period
#'
#' @param elog data.frame with columns 'cust' and 'date'
#' @param per time unit, either 'week', 'day', 'hour', 'min', 'sec'
#' @param T.cal end date of calibration period
#' @param T.tot end date of holdout period
#' @return data.frame
#' @export
#' @import data.table
elog2cbs <- function(elog, per = "week", T.cal = max(elog$date), T.tot = max(elog$date)) {
  cust <- first <- itt <- T.star <- x.star <- NULL # avoid checkUsage warnings

  is.dt <- is.data.table(elog)
  # convert to data.table for improved performance
  elog_dt <- data.table(elog)
  # merge same dates
  setkey(elog_dt, cust, date)
  elog_dt <- unique(elog_dt)
  # determine time since first date for each customer
  elog_dt[, first := min(date), by="cust"]
  elog_dt[, t := as.numeric(difftime(date, first, units=per)), by="cust"]
  # compute intertransaction times
  elog_dt[, itt := c(0, diff(t)), by="cust"]
  # count events in calibration period
  cbs <- elog_dt[date<=T.cal,
                 list(x=.N-1, t.x=max(t), litt=sum(log(itt[itt>0]))),
                 by="cust,first"]
  cbs[, T.cal := as.numeric(difftime(T.cal, first, units=per))]
  cbs[, T.star := as.numeric(difftime(T.tot, first, units=per)) - T.cal]
  cbs[, first := NULL]
  setkey(cbs, cust)
  # count events in validation period
  val <- elog_dt[date>T.cal & date<=T.tot,
              list(x.star=.N),
              keyby="cust"]
  cbs <- merge(cbs, val, all.x=TRUE, by="cust")
  cbs[is.na(x.star), x.star := 0]
  # return same object type as was passed
  if (!is.dt) cbs <- data.frame(cbs)
  return(cbs)
}
