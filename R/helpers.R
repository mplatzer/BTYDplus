
#' Estimate Regularity in Intertransaction Timings
#' 
#' Estimates degree of regularity of intertransaction timings of a customer
#' cohort. This is done by 1) assuming same regularity across all customers 
#' (\code{method = "wheat"}), or 2) by estimating regularity for each customer 
#' seperately, as the shape parameter of a fitted gamma distribution, and then 
#' return the median across estimates; this requires sufficient (>=10) 
#' transactions per customer
#' 
#' @param elog Event log, a \code{data.frame} with columns \code{cust} and
#'   transaction time \code{t} or \code{date}
#' @param method Either \code{wheat}, \code{mle}, \code{mle-minka}, \code{mle-thom} or
#'   \code{cv}.
#' @param plot If \code{TRUE} then distribution of estimated regularity will be
#'   plotted.
#' @return Estimated real-valued regularity parameter.
#' @references Wheat, Rita D., and Donald G. Morrison.  'Estimating purchase
#'   regularity with two interpurchase times.'
#' @references Dunn, Richard, Steven Reader, and Neil Wrigley. 'An investigation
#'   of the assumptions of the NBD model' Applied Statistics (1983): 249-259.
#' @references Wu, Couchen, and H-L. Chen. 'A consumer purchasing model with
#'   learning and departure behaviour.'  Journal of the Operational Research
#'   Society (2000): 583-591.
#' @references
#'   \url{http://research.microsoft.com/en-us/um/people/minka/papers/minka-gamma.pdf}
#'   
#' @export
#' @examples
#' elog <- cdnow.sample()$elog[, c("cust", "date")]
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


#' Plot timing patterns of sampled customers.
#' 
#' @param elog Event log, a \code{data.frame} with columns \code{cust} and 
#'   transaction time \code{t} or \code{date}.
#' @param T.cal End of calibration period, which is visualized as a vertical line.
#' @param n Number of sampled customers.
#' @param title Plot title.
#' @export
#' @examples
#' elog <- cdnow.sample()$elog
#' plotSampledTimingPatterns(elog, T.cal = "1997-09-30")
plotSampledTimingPatterns <- function(elog, T.cal = NULL, n = 50, title = "Sampled Timing Patterns") {
  
  cust <- first <- t <- V1 <- NULL  # suppress checkUsage warnings
  elog_dt <- setDT(copy(elog))
  custs <- sample(unique(elog_dt$cust), size = min(n, uniqueN(elog_dt$cust)), replace = FALSE)
  n <- length(custs)
  if (!"t" %in% names(elog_dt)) elog_dt[, `:=`(t, as.numeric(date))]
  rg <- range(elog_dt$t)
  elog_dt <- elog_dt[cust %in% custs]
  elog_dt[, `:=`(first, min(t)), by = "cust"]
  if (!is.character(elog_dt$cust)) elog_dt[, `:=`(cust, as.character(cust))]
  custs <- elog_dt[, min(date), by = "cust"][order(V1), cust]
  setkeyv(elog_dt, "cust")
  op <- par(mar = c(0.5,0.5,2.5,0.5))
  ymax <- ifelse(is.null(T.cal), n, ceiling(n * 1.05))
  plot(1, xlim = rg, ylim = c(1, ymax), typ = "n", 
       axes = FALSE, frame = FALSE, 
       xlab = "", ylab = "",
       main = title)
  for (i in 1:n) {
    ts <- unique(elog_dt[custs[i], t])
    segments(min(ts), i, rg[2], i, col = "#efefef", lty = 1, lwd = 1)
    points(min(ts), i, pch = 16, col = "#454545")
    points(ts[ts>min(ts)], rep(i, length(ts)-1), pch = 1, col = "#454545")
  }
  par(op)
  if (!is.null(T.cal)) {
    if ("date" %in% names(elog_dt)) T.cal <- as.Date(T.cal)
    T.cal <- as.numeric(T.cal)
    abline(v = T.cal)
    text("Calibration",  x = T.cal - (T.cal - rg[1]) / 2, y = ymax, col = "#454545")
    text("Holdout",  x = T.cal + (rg[2] - T.cal) / 2, y = ymax, col = "#454545")
  }
  
}


#' Convernt Event Log to customer-level summary statistic.
#' 
#' Takes the event log of a customer cohort, and returns a sufficient summary statistic for applying common BTYD models.
#' 
#' Note: compared to \code{\link[BTYD]{dc.ElogToCbsCbt}} this also adds a summary statistic for estimating regularity.
#'
#' Customers without any transaction during calibration period are being dropped from the result. Transactions with identical \code{cust} and \code{date} field are treated as a single transaction, with `sales` being summed up
#'
#' @param elog Event log, a \code{data.frame} with columns \code{cust} and 
#'   transaction time \code{t} or \code{date}. If column \code{sales} is
#'   present, it will be aggregated as well.
#' @param per Time unit, either 'week', 'day', 'hour', 'min', 'sec'.
#' @param T.cal End date of calibration period.
#' @param T.tot End date of holdout period.
#' @return data.frame with fields
#'  \item{\code{cust}}{customer id (unique key)}
#'  \item{\code{x}}{number of recurring events in calibration period}
#'  \item{\code{t.x}}{time between first and last event in calibration period}
#'  \item{\code{litt}}{sum of logarithmic intertransaction timings durint calibration period}
#'  \item{\code{sales}}{sum of sales in calibration period}
#'  \item{\code{first}}{date of first transaction in calibration period}
#'  \item{\code{T.cal}}{time between first event and end of calibration period}
#'  \item{\code{T.star}}{length of holdout period}
#'  \item{\code{x.star}}{number of events within holdout period}
#'  \item{\code{sales.star}}{sum of sales within holdout period}
#' @export
#' @examples
#' elog <- cdnow.sample()$elog
#' cbs <- elog2cbs(elog, T.cal = "1998-01-01", T.tot = "1998-03-31")
#' head(cbs)
elog2cbs <- function(elog, per = "week", T.cal = max(elog$date), T.tot = max(elog$date)) {
  cust <- first <- itt <- T.star <- x.star <- sales <- sales.star <- NULL  # suppress checkUsage warnings
  stopifnot(inherits(elog, "data.frame"))
  if (!all(c("cust", "date") %in% names(elog))) stop("`elog` must have fields `cust` and `date`")
  if (!any(c("Date", "POSIXt") %in% class(elog$date))) stop("`date` field must be of class `Date` or `POSIXt`")
  if (is.character(T.cal)) T.cal <- as.Date(T.cal)
  if (is.character(T.tot)) T.tot <- as.Date(T.tot)
  stopifnot(T.tot >= T.cal)
  
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
#' This is a convenience wrapper for \code{data('cdnowElog', package='BTYD')},
#' with same-day transactions being merged together for each customer.
#' 
#' @references Fader, Peter S. and Bruce G.,S. Hardie, (2001), 'Forecasting
#'   Repeat Sales at CDNOW: A Case Study,' Interfaces, 31 (May-June), Part 2 of
#'   2, S94-S107.
#' @param T.cal Final date of calibration period.
#' @param T.tot Final date of holdout period.
#' @return List of length 2:
#' \item{\code{cbs}}{A data.frame with a row for each customer and the summary statistic as columns, generated with \code{\link{elog2cbs}}}
#' \item{\code{elog}}{A data.frame with a row for each transaction, and columns \code{cust}, \code{date}, \code{sales} and \code{cds}.}
#' @export
#' @examples
#' cdnow <- cdnow.sample()
#' head(cdnow$cbs)
#' head(cdnow$elog)
cdnow.sample <- function(T.cal = "1997-09-30", T.tot = "1998-06-30") {
  cds <- sales <- NULL  # suppress checkUsage warnings
  elog <- fread(system.file("data/cdnowElog.csv", package = "BTYD"))
  setnames(elog, "masterid", "cust")
  elog <- elog[, list(sales = sum(sales), cds = sum(cds)), by = "cust,date"]
  elog[, `:=`(date, as.Date(as.character(date), "%Y%m%d"))]
  cbs <- elog2cbs(elog, per = "week", T.cal = T.cal, T.tot = T.tot)
  return(list(elog = as.data.frame(elog), cbs = as.data.frame(cbs)))
}


#' Aggregate Event Log to cumulative number of repeat transactions
#' 
#' @param elog Event log, a \code{data.frame} with columns \code{cust} and 
#'   transaction time \code{t} or \code{date}.
#' @param by Only return every \code{}-th number. Defaults to 7, and thus
#'   returns weekly numbers.
#' @return Numeric vector of cumulative repeat transactions.
#' @export
#' @seealso \code{\link{elog2inc}}
#' @examples 
#' elog <- cdnow.sample()$elog
#' cum <- elog2cum(elog)
#' plot(cum, typ="l")
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


#' Aggregate Event Log to incremental number of repeat transactions
#' 
#' @param elog Event log, a \code{data.frame} with columns \code{cust} and 
#'   transaction time \code{t} or \code{date}.
#' @param by Only return every \code{}-th number. Defaults to 7, and thus
#'   returns weekly numbers.
#' @return Numeric vector of repeat transactions.
#' @export
#' @seealso \code{\link{elog2cum}}
#' @examples
#' elog <- cdnow.sample()$elog
#' inc <- elog2inc(elog)
#' plot(inc, typ="l")
elog2inc <- function(elog, by = 7) {
  cum <- elog2cum(elog = elog, by = by)
  return(diff(cum))
}


#' Check Model Parameters
#' 
#' Wrapper for \code{BTYD::dc.check.model.params} with additional check for
#' parameter names if these are present
#'
#' @keywords internal
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


#' Generic Method for Tracking Plots
#' 
#' @keywords internal
dc.PlotTracking <- function(actual, expected, T.cal = NULL,
                            xlab = "", ylab = "", title = "", 
                            xticklab = NULL, ymax = NULL) {

  stopifnot(is.numeric(actual))
  stopifnot(is.numeric(expected))
  stopifnot(all(actual >= 0))
  stopifnot(all(expected >= 0))
  stopifnot(length(actual) == length(expected))
  
  if (is.null(ymax)) ymax <- max(c(actual, expected)) * 1.05
  plot(actual, type = "l", xaxt = "n", xlab = xlab, ylab = ylab, col = 1, ylim = c(0, ymax), main = title)
  lines(expected, lty = 2, col = 2)
  if (is.null(xticklab)) {
    axis(1, at = 1:length(actual), labels = 1:length(actual))
  } else {
    if (length(actual) != length(xticklab)) {
      stop("Plot error, xticklab does not have the correct size")
    }
    axis(1, at = 1:length(actual), labels = xticklab)
  }
  if (!is.null(T.cal)) abline(v = max(T.cal), lty = 2)
  pos <- ifelse(which.max(expected) == length(expected), "bottomright", "topright")
  legend(pos, legend = c("Actual", "Model"), col = 1:2, lty = 1:2, lwd = 1)
  return(rbind(actual, expected))
}
