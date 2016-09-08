
#' Estimate Regularity in Intertransaction Timings
#' 
#' Estimates degree of regularity of intertransaction timings of a customer cohort,
#' 
#' This is done
#'  - by either assuming same regularity across all customers; this then only 
#'      require three transactions per customer
#'      method: wheat
#'  - or by estimating regularity for each customer seperately (as the shape
#'      parameter of a fitted gamma distribution), and then return the 
#'      median across estimates; this requires min. 10 transactions per customer
#'      methods: mle, mle-minka, mle-thom, cv
#' 
#' @param elog data.frame with transaction logs; requires columns customer-id
#'   'cust' and transaction time \code{t} or 'date'
#' @param method \code{wheat}, \code{mle}, \code{mle-minka}, \code{mle-thom}, \code{cv}
#' @param plot if \code{TRUE} then distribution of estimated regularity will be plotted
#' @return estimated real-valued regularity parameter; rounded to an integer,
#'   this can be used as \code{k} for estimating MBG/CNBD-k models
#' @references Wheat, Rita D., and Donald G. Morrison.  'Estimating purchase regularity with two interpurchase times.'
#' @export
#' @examples
#' # generate Erlang-3 with various lambdas
#' n <- 100
#' k <- 3
#' x <- 20
#' elog <- do.call(rbind, lapply(1:n, function(i) {
#'   lambda <- exp(rnorm(1))
#'     data.frame(cust = i, t = cumsum(rgamma(x, k, k * lambda)))
#'     }))
#'     
#' # estimate regularity parameter k
#' estimateRegularity(elog, plot = TRUE, method = 'wheat')
#' estimateRegularity(elog, plot = TRUE, method = 'mle-minka')
#' estimateRegularity(elog, plot = TRUE, method = 'mle-thom')
#' estimateRegularity(elog, plot = TRUE, method = 'cv') 
estimateRegularity <- function(elog, method = "wheat", plot = FALSE) {
  if (!"cust" %in% names(elog)) 
    stop("Error in estimateRegularity: elog must have a column labelled \"cust\"")
  if (!"date" %in% names(elog) & !"t" %in% names(elog)) 
    stop("Error in estimateRegularity: elog must have a column labelled \"t\" or \"date\"")
  if (!"t" %in% names(elog)) 
    elog$t <- as.numeric(elog$date)
  trans <- split(elog, elog$cust)
  if (method == "wheat") {
    # Wheat, Rita D., and Donald G. Morrison.  'Estimating purchase regularity with two interpurchase times.'
    # Journal of Marketing Research (1990): 87-93.
    M <- unlist(lapply(trans, function(df) {
      itt <- diff(sort(unique(df$t)))
      if (length(itt) > 1) {
        # take last two itt's to capture most recent regularity
        itt2 <- rev(itt)[1:2]
        return(sample(itt2)[1]/sum(itt2))
      }
    }))
    r <- (1 - 4 * var(M))/(8 * var(M))
    if (plot) {
      op <- par(mar = c(1, 2, 1, 2))
      plot(density(M), main = "", sub = "", xlab = "", ylab = "", lwd = 2, frame = FALSE, axes = FALSE)
      polygon(density(M), col = "lightgray", border = 1)
      fn1 <- function(x) dbeta(x, 1, 1)
      fnr <- function(x) dbeta(x, round(r), round(r))
      curve(fn1, add = TRUE, lty = 2, lwd = 2)
      curve(fnr, add = TRUE, lty = 2, lwd = 2)
      par(op)
    }
    return(r)
    
  } else {
    if (method == "mle" | method == "mle-minka") {
      # Maximum Likelihood Estimator http://en.wikipedia.org/wiki/Gamma_distribution#Maximum_likelihood_estimation
      # Approximation for MLE by Minka http://research.microsoft.com/en-us/um/people/minka/papers/minka-gamma.pdf
      ks <- unlist(lapply(trans, function(df) {
        itt <- diff(sort(unique(df$t)))
        if (length(itt) >= 9) {
          s <- log(sum(itt)/length(itt)) - sum(log(itt))/length(itt)
          if (method == "mle") {
          fn <- function(v) {
            return((log(v) - digamma(v) - s)^2)
          }
          k <- optimize(fn, lower = 0.1, upper = 50)$min
          } else if (method == "mle-minka") {
          k <- (3 - s + sqrt((s - 3)^2 + 24 * s))/(12 * s)
          }
          return(k)
        }
      }))
      
    } else if (method == "mle-thom") {
      # Approximation for ML estimator Thom (1968); see Dunn, Richard, Steven Reader, and Neil Wrigley.  'An
      # investigation of the assumptions of the NBD model' Applied Statistics (1983): 249-259.
      ks <- unlist(lapply(trans, function(df) {
        itt <- diff(sort(unique(df$t)))
        if (length(itt) >= 9) {
          hm <- function(v) exp(sum(log(v))/length(v))
          mu <- log(mean(itt)/hm(itt))
          d <- (1/(4 * mu)) * (1 + sqrt(1 + 4 * mu/3))
          return(d)
        }
      }))
      
    } else if (method == "cv") {
      # Estimate regularity by analyzing coefficient of variation Wu, Couchen, and H-L. Chen. 'A consumer purchasing
      # model with learning and departure behaviour.'  Journal of the Operational Research Society (2000): 583-591.
      ks <- unlist(lapply(trans, function(df) {
        itt <- diff(sort(unique(df$t)))
        if (length(itt) >= 9) {
          cv <- sd(itt)/mean(itt)
          k <- 1/cv^2
          return(k)
        }
      }))
    }
    if (length(ks) == 0) 
      stop("No customers with 10 or more transactions.")
    
    if (plot) {
      ymax <- median(ks) * 3
      boxplot(ks, horizontal = TRUE, ylim = c(0, ymax), frame = FALSE, axes = FALSE)
      axis(1, at = 0:ymax)
      axis(3, at = 0:ymax, labels = rep("", 1 + ymax))
      abline(v = 1:ymax, lty = "dotted", col = "lightgray")
      boxplot(ks, horizontal = TRUE, add = TRUE, col = "gray", frame = FALSE, axes = FALSE)
    }
    
    return(median(ks))
  }
}


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
  cust <- first <- itt <- T.star <- x.star <- sales <- sales.star <- NULL  # suppress checkUsage warnings
  stopifnot(inherits(elog, "data.frame"))
  if (!all(c("cust", "date") %in% names(elog))) stop("`elog` must have fields `cust` and `date`")
  if (!any(c("Date", "POSIXt") %in% class(elog$date))) stop("`date` field must be of class `Date` or `POSIXt`")
  
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


#' CDNow Sample Data
#' 
#' This is a convenience wrapper for data('cdnowElog', package='BTYD'), with
#' same-day transactions being merged together for each customer.
#' 
#' @references Fader, Peter S. and Bruce G.,S. Hardie, (2001), 'Forecasting
#'   Repeat Sales at CDNOW: A Case Study,' Interfaces, 31 (May-June), Part 2 of
#'   2, S94-S107.
#' @return list with two data.frames: `elog` for the event logs and `cbs` for
#'   the customer by sufficient summary
#' @seealso \code{\link{elog2cbs}} 
#' @export
cdnow.sample <- function() {
  cds <- sales <- NULL  # suppress checkUsage warnings
  elog <- fread(system.file("data/cdnowElog.csv", package = "BTYD"))
  setnames(elog, "masterid", "cust")
  elog <- elog[, list(sales = sum(sales), cds = sum(cds)), by = "cust,date"]
  elog[, `:=`(date, as.Date(as.character(date), "%Y%m%d"))]
  cbs <- elog2cbs(elog, per = "week", T.cal = as.Date("1997-09-30"))
  return(list(elog = as.data.frame(elog), cbs = as.data.frame(cbs)))
}


#' Convert Event Log to (weekly) cumulative number of repeat transactions
#' 
#' @param elog Event Log
#' @param by defaults to 7, which means weekly numbers
#' @return numeric vector
#' @export
#' @seealso \code{\link{elog2inc}}
#' @examples 
#' elog <- cdnow.sample()$elog
#' cum <- elog2cum(elog)
#' plot(cum, typ='l')
elog2cum <- function(elog, by = 7) {
  t0 <- sales <- NULL  # suppress checkUsage warnings
  stopifnot("cust" %in% names(elog))
  is.dt <- is.data.table(elog)
  if (!is.dt) 
    elog <- as.data.table(elog) else elog <- copy(elog)
  if (!"t" %in% names(elog)) {
    stopifnot("date" %in% names(elog))
    cohort_start <- min(as.numeric(elog$date))
    elog[, `:=`(t, as.numeric(date) - cohort_start)]
  }
  elog[, `:=`(t0, min(t)), by = "cust"]
  inc <- elog[t > t0, .N, keyby = list(t = ceiling(t))]$N
  cum <- c(0, cumsum(inc)[seq(by, length(inc) - 1, by = by)], sum(inc))
  return(cum)
}


#' Convert Event Log to (weekly) incremental number of repeat transactions
#' 
#' @param elog Event Log
#' @param by defaults to 7, which means weekly numbers
#' @return numeric vector
#' @export
#' @seealso \code{\link{elog2cum}}
#' @examples
#' elog <- cdnow.sample()$elog
#' inc <- elog2inc(elog)
#' plot(inc, typ='l')
elog2inc <- function(elog, by = 7) {
  cum <- elog2cum(elog = elog, by = by)
  return(diff(cum))
}


#' Wrapper for BTYD::dc.check.model.params with additional check for parameter
#' names if these are present
#'
#' @param printnames Names to print parameter errors
#' @param params model parameters
#' @param func Function calling dc.check.model.params.safe
#' @return stops program if there is something wrong with the parameters
dc.check.model.params.safe <- function(printnames, params, func) {
  # first do basic checks
  dc.check.model.params(printnames, params, func)
  # then check for names, if these are present
  if (!is.null(names(params))) {
    idx <- names(params) != ""
    if (any(printnames[idx] != names(params)[idx])) {
      stop("Error in ", func, ": Parameter names do not match - ", paste0(printnames, collapse = ","), " != ", 
        paste0(names(params), collapse = ","), call. = FALSE)
    }
  }
}
