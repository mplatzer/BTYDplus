
#' Estimate Regularity in Intertransaction Timings
#'
#' The models (M)BG/CNBD-k and Pareto/GGG are capable of leveraging regularity
#' within transaction timings for improving forecast accuracy. This method
#' provides a quick check for the degree of regularity in the event timings. A
#' return value of close to 1 supports the assumption of exponentially
#' distributed intertransaction times, whereas values significantly larger than
#' 1 reveal the presence of regularity.
#'
#' Estimation is either done by 1) assuming the same degree of regularity across
#' all customers (Wheat & Morrison (1990) via \code{method = "wheat"}), or 2) by
#' estimating regularity for each customer separately, as the shape parameter of
#' a fitted gamma distribution, and then return the median across estimates. The
#' latter methods, though, require sufficient (>=\code{min}) transactions per
#' customer.
#'
#' Wheat & Morrison (1990)'s method calculates for each customer a statistic
#' \code{M} based on her last two number of intertransaction times as
#' \code{ITT_1 / (ITT_1 + ITT_2)}. That measure is known to follow a
#' \code{Beta(k, k)} distribution, and \code{k} can be estimated as
#' \code{(1-4*Var(M))/(8*Var(M))}. The corresponding diagnostic plot (\code{plot
#' = TRUE}) shows the actual distribution of \code{M} vs. the theoretical
#' distribution for \code{k = 1} and \code{k = 2}.
#'
#' @param elog Event log, a \code{data.frame} with columns \code{cust} and
#'   transaction time \code{t} or \code{date}
#' @param method Either \code{wheat}, \code{mle}, \code{mle-minka}, \code{mle-thom} or
#'   \code{cv}.
#' @param plot If \code{TRUE} then an additional diagnostic plot is provided.
#' @param title Plot title.
#' @param min Minimum number of intertransaction times per customer. Customers
#'   with less than \code{min} intertransactions are not considered. Defaults to 2
#'   for method `wheat`, and to 10 otherwise.
#' @return Estimated real-valued regularity parameter.
#' @references Wheat, Rita D., and Donald G. Morrison. "Estimating purchase
#'   regularity with two interpurchase times." Journal of Marketing Research
#'   (1990): 87-93.
#' @references Dunn, Richard, Steven Reader, and Neil Wrigley. 'An investigation
#'   of the assumptions of the NBD model' Applied Statistics (1983): 249-259.
#' @references Wu, Couchen, and H-L. Chen. 'A consumer purchasing model with
#'   learning and departure behaviour.'  Journal of the Operational Research
#'   Society (2000): 583-591.
#' @references
#'   \url{https://tminka.github.io/papers/minka-gamma.pdf}
#'
#' @export
#' @examples
#' data("groceryElog")
#' estimateRegularity(groceryElog, plot = TRUE, method = 'wheat')
#' estimateRegularity(groceryElog, plot = TRUE, method = 'mle-minka')
#' estimateRegularity(groceryElog, plot = TRUE, method = 'mle-thom')
#' estimateRegularity(groceryElog, plot = TRUE, method = 'cv')
estimateRegularity <- function(elog, method = "wheat", plot = FALSE, title = "", min = NULL) {
  N <- t <- NULL # suppress checkUsage warnings
  if (!"cust" %in% names(elog))
    stop("elog must have a column labelled \"cust\"")
  if (!"date" %in% names(elog) & !"t" %in% names(elog))
    stop("elog must have a column labelled \"t\" or \"date\"")
  if (!"t" %in% names(elog))
    elog$t <- as.numeric(elog$date)
  elog_dt <- subset(setDT(copy(elog)), select = c("cust", "t"))
  setkey(elog_dt)
  elog_dt <- unique(elog_dt) # nolint
  # discard customers with less than `min` ITTs
  if (is.null(min)) {
    min <- ifelse(method == "wheat", 2, 10)
  } else {
    stopifnot(is.numeric(min))
    stopifnot(min >= 2)
  }
  elog_dt[, N := .N, by = "cust"]
  elog_dt <- elog_dt[N > min]
  if (nrow(elog_dt) == 0) stop("No customers with sufficient number of transactions.")
  # calculate method specific estimate
  if (method == "wheat") {
    # Wheat, Rita D., and Donald G. Morrison.  'Estimating purchase regularity
    # with two interpurchase times.' Journal of Marketing Research (1990): 87-93.
    setkeyv(elog_dt, c("cust", "t"))
    calcM <- function(itts) sample(utils::tail(itts, 2), 1) / sum(utils::tail(itts, 2))
    M <- elog_dt[, calcM(diff(t)), by = "cust"]$V1
    if (var(M) == 0) stop("No customers with sufficient number of transactions.")
    r <- (1 - 4 * var(M)) / (8 * var(M))
    if (plot) {
      mar_top <- ifelse(title != "", 2.5, 1)
      op <- par(mar = c(1, 1, mar_top, 1))
      M_density <- density(M)
      ymax <- max(M_density$y, 1.5)
      plot(M_density, xlim = c(-0.05, 1.05), ylim = c(0, ymax),
           main = title, sub = "", xlab = "", ylab = "",
           lwd = 2, frame = FALSE, axes = FALSE)
      polygon(M_density,
              col = "lightgray", border = 1)
      fn1 <- function(x) dbeta(x, 1, 1)
      fn2 <- function(x) dbeta(x, 2, 2)
      curve(fn1, add = TRUE, lty = 2, lwd = 1, col = "gray12")
      curve(fn2, add = TRUE, lty = 3, lwd = 1, col = "gray12")
      par(op)
    }
    return(r)

  } else {
    if (method == "mle") {
      # Maximum Likelihood Estimator
      # https://en.wikipedia.org/wiki/Gamma_distribution#Maximum_likelihood_estimation
      est_k <- function(itts) {
        s <- log(sum(itts) / length(itts)) - sum(log(itts)) / length(itts)
        fn <- function(v) return( (log(v) - digamma(v) - s) ^ 2)
        k <- optimize(fn, lower = 0.1, upper = 50)$min
        return(k)
      }
    } else if (method == "mle-minka") {
      # Approximation for MLE by Minka
      # https://tminka.github.io/papers/minka-gamma.pdf
      est_k <- function(itts) {
        s <- log(sum(itts) / length(itts)) - sum(log(itts)) / length(itts)
        k <- (3 - s + sqrt( (s - 3) ^ 2 + 24 * s)) / (12 * s)
        return(k)
      }
    } else if (method == "mle-thom") {
      # Approximation for ML estimator Thom (1968); see Dunn, Richard, Steven
      # Reader, and Neil Wrigley. 'An investigation of the assumptions of the
      # NBD model' Applied Statistics (1983): 249-259.
      est_k <- function(itts) {
        hm <- function(v) exp(sum(log(v)) / length(v))
        mu <- log(mean(itts) / hm(itts))
        k <- (1 / (4 * mu)) * (1 + sqrt(1 + 4 * mu / 3))
        return(k)
      }
    } else if (method == "cv") {
      # Estimate regularity by analyzing coefficient of variation Wu, Couchen,
      # and H-L. Chen. 'A consumer purchasing model with learning and departure
      # behaviour.'  Journal of the Operational Research Society (2000): 583-591.
      est_k <- function(itts) {
        cv <- sd(itts) / mean(itts)
        k <- 1 / cv ^ 2
        return (k)
      }
    }
    ks <- elog_dt[, est_k(diff(t)), by = "cust"]$V1
    if (plot) {
      ymax <- median(ks) * 3
      suppressWarnings(boxplot(ks, horizontal = TRUE, ylim = c(0, ymax),
                               frame = FALSE, axes = FALSE, main = title))
      axis(1, at = 0:ymax)
      axis(3, at = 0:ymax, labels = rep("", 1 + ymax))
      abline(v = 1:ymax, lty = "dotted", col = "lightgray")
      suppressWarnings(boxplot(ks, horizontal = TRUE, add = TRUE,
                               col = "gray", frame = FALSE, axes = FALSE))
    }
    return(median(ks))
  }
}


#' Plot timing patterns of sampled customers
#'
#' @param elog Event log, a \code{data.frame} with columns \code{cust} and
#'   transaction time \code{t} or \code{date}.
#' @param n Number of sampled customers.
#' @param T.cal End of calibration period, which is visualized as a vertical line.
#' @param T.tot End of observation period
#' @param title Plot title.
#' @param headers Vector of length 2 for adding headers to the plot, e.g.
#'   \code{c("Calibration", "Holdout")}.
#' @export
#' @examples
#' data("groceryElog")
#' plotTimingPatterns(groceryElog, T.tot = "2008-12-31")
#' plotTimingPatterns(groceryElog, T.cal = "2006-12-31", headers = c("Calibration", "Holdout"))
plotTimingPatterns <- function(elog, n = 40, T.cal = NULL, T.tot = NULL,
                               title = "Sampled Timing Patterns", headers = NULL) {

  cust <- first <- t <- V1 <- NULL  # suppress checkUsage warnings
  elog_dt <- setDT(copy(elog))
  custs <- sample(unique(elog_dt$cust), size = min(n, uniqueN(elog_dt$cust)), replace = FALSE)
  n <- length(custs)
  if (!"t" %in% names(elog_dt)) elog_dt[, `:=`(t, as.numeric(date))]
  T.0 <- min(elog_dt$t)
  if (is.null(T.cal)) {
    T.cal <- max(elog_dt$t)
  } else if (!is.numeric(T.cal)) {
    T.cal <- as.numeric(as.Date(T.cal))
  }
  if (is.null(T.tot)) {
    T.tot <- max(elog_dt$t)
  } else if (!is.numeric(T.tot)) {
    T.tot <- as.numeric(as.Date(T.tot))
  }
  elog_dt <- elog_dt[cust %in% custs & t <= T.tot]
  elog_dt[, `:=`(first, min(t)), by = "cust"]
  if (!is.character(elog_dt$cust)) elog_dt[, `:=`(cust, as.character(cust))]
  custs <- elog_dt[, min(date), by = "cust"][order(V1), cust]
  setkeyv(elog_dt, "cust")
  mar_top <- ifelse(is.null(title) || title == "", 0.5, 2.5)
  op <- par(mar = c(0.5, 0.5, mar_top, 0.5))
  ymax <- ifelse(is.character(headers), ceiling(n * 1.05), n)
  plot(1, xlim = c(T.0, T.tot), ylim = c(1, ymax), typ = "n",
       axes = FALSE, frame = FALSE,
       xlab = "", ylab = "",
       main = title)
  for (i in 1:n) {
    ts <- unique(elog_dt[custs[i], t])
    segments(min(ts), i, T.tot, i, col = "#efefef", lty = 1, lwd = 1)
    points(min(ts), i, pch = 21, col = "#454545", bg = "#454545", cex = 0.8)
    ts.cal <- ts[ts > min(ts) & ts <= as.numeric(T.cal)]
    ts.val <- ts[ts > as.numeric(T.cal)]
    points(ts.cal, rep(i, length(ts.cal)), pch = 21, col = "#454545", bg = "#454545", cex = 0.7)
    points(ts.val, rep(i, length(ts.val)), pch = 21, col = "#454545", bg = "#999999", cex = 0.7)
  }
  par(op)
  if (T.cal < T.tot) abline(v = T.cal, lwd = 1.5, col = "#454545")
  if (is.character(headers)) {
    text(headers[1],  x = T.cal - (T.cal - T.0) / 2, y = ymax, col = "#454545", cex = 0.8)
    if (T.cal < T.tot) text(headers[2], x = T.cal + (T.tot - T.cal) / 2, y = ymax, col = "#454545", cex = 0.8)
  }
  invisible()
}


#' Convert Event Log to customer-level summary statistic
#'
#' Efficient implementation for the conversion of an event log into a
#' customer-by-sufficient-statistic (CBS) \code{data.frame}, with a row for each
#' customer, which is the required data format for estimating model parameters.
#'
#' The time unit for expressing \code{t.x}, \code{T.cal} and \code{litt} are
#' determined via the argument \code{units}, which is passed forward to method
#' \code{difftime}, and defaults to \code{weeks}.
#'
#' Argument \code{T.tot} allows one to specify the end of the observation period,
#' i.e. the last possible date of an event to still be included in the event
#' log. If  \code{T.tot} is not provided, then the date of the last recorded event
#' will be assumed to coincide with the end of the observation period. If
#'  \code{T.tot} is provided, then any event that occurs after that date is discarded.
#'
#' Argument \code{T.cal} allows one to split the summary statistics into a
#' calibration and a holdout period. This can be useful for evaluating
#' forecasting accuracy for a given dataset. If \code{T.cal} is not provided,
#' then the whole observation period is considered, and is then subsequently
#' used for for estimating model parameters. If it is provided, then the
#' returned \code{data.frame} contains two additional fields, with \code{x.star}
#' representing the number of repeat transactions during the holdout period of
#' length \code{T.star}. And only those customers are contained, who have had at
#' least one event during the calibration period.
#'
#' Transactions with identical \code{cust} and \code{date} field are treated as
#' a single transaction, with \code{sales} being summed up.
#'
#' @param elog Event log, a \code{data.frame} with field \code{cust} for the
#'   customer ID and field \code{date} for the date/time of the event, which
#'   should be of type \code{Date} or \code{POSIXt}. If a field \code{sales} is
#'   present, it will be aggregated as well.
#' @param units Time unit, either \code{week}, \code{day}, \code{hour},
#'   \code{min} or \code{sec}. See \code{\link[base]{difftime}}.
#' @param T.cal End date of calibration period. Defaults to
#'   \code{max(elog$date)}.
#' @param T.tot End date of the observation period. Defaults to
#'   \code{max(elog$date)}.
#' @return \code{data.frame} with fields:
#'  \item{\code{cust}}{Customer id (unique key).}
#'  \item{\code{x}}{Number of recurring events in calibration period.}
#'  \item{\code{t.x}}{Time between first and last event in calibration period.}
#'  \item{\code{litt}}{Sum of logarithmic intertransaction timings during calibration period.}
#'  \item{\code{sales}}{Sum of sales in calibration period, incl. initial transaction. Only if \code{elog$sales} is provided.}
#'  \item{\code{sales.x}}{Sum of sales in calibration period, excl. initial transaction. Only if \code{elog$sales} is provided.}
#'  \item{\code{first}}{Date of first transaction in calibration period.}
#'  \item{\code{T.cal}}{Time between first event and end of calibration period.}
#'  \item{\code{T.star}}{Length of holdout period. Only if \code{T.cal} is provided.}
#'  \item{\code{x.star}}{Number of events within holdout period. Only if \code{T.cal} is provided.}
#'  \item{\code{sales.star}}{Sum of sales within holdout period. Only if \code{T.cal} and \code{elog$sales} are provided.}
#' @export
#' @examples
#' data("groceryElog")
#' cbs <- elog2cbs(groceryElog, T.cal = "2006-12-31", T.tot = "2007-12-30")
#' head(cbs)
elog2cbs <- function(elog, units = "week", T.cal = NULL, T.tot = NULL) {
  cust <- first <- itt <- T.star <- x.star <- sales <- sales.star <- sales.x <- NULL  # suppress checkUsage warnings
  stopifnot(inherits(elog, "data.frame"))
  is.dt <- is.data.table(elog)
  if (nrow(elog) == 0) {
    cbs <- data.frame(cust = character(0),
                      x = numeric(0),
                      t.x = numeric(0),
                      litt = numeric(0),
                      first = as.Date(character(0)),
                      T.cal = numeric(0))
    if (is.dt) cbs <- as.data.table(cbs)
    return(cbs)
  }
  if (!all(c("cust", "date") %in% names(elog))) stop("`elog` must have fields `cust` and `date`")
  if (!any(c("Date", "POSIXt") %in% class(elog$date))) stop("`date` field must be of class `Date` or `POSIXt`")
  if ("sales" %in% names(elog) & !is.numeric(elog$sales)) stop("`sales` field must be numeric")
  if (is.data.frame(elog) && nrow(elog) == 0) stop("data.frame supplied to elog2cbs must not be empty")
  if (is.null(T.cal)) T.cal <- max(elog$date)
  if (is.null(T.tot)) T.tot <- max(elog$date)
  if (is.character(T.cal)) T.cal <- if (class(elog$date)[1] == "Date") as.Date(T.cal) else as.POSIXct(T.cal)
  if (is.character(T.tot)) T.tot <- if (class(elog$date)[1] == "Date") as.Date(T.tot) else as.POSIXct(T.tot)
  if (T.tot < T.cal) T.tot <- T.cal
  stopifnot(T.tot >= min(elog$date))

  has.holdout <- T.cal < T.tot
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
  elog_dt[, `:=`(t, as.numeric(difftime(date, first, units = units))), by = "cust"]
  # compute intertransaction times
  elog_dt[, `:=`(itt, c(0, diff(t))), by = "cust"]
  # count events in calibration period
  cbs <- elog_dt[date <= T.cal,
                 list(x = .N - 1,
                      t.x = max(t),
                      litt = sum(log(itt[itt > 0])),
                      sales = sum(sales),
                      sales.x = sum(sales[t > 0])),
                 by = "cust,first"]
  cbs[, `:=`(T.cal, as.numeric(difftime(T.cal, first, units = units)))]
  setkey(cbs, cust)
  # count events in validation period
  if (has.holdout) {
    cbs[, `:=`(T.star, as.numeric(difftime(T.tot, first, units = units)) - T.cal)]
    val <- elog_dt[date > T.cal & date <= T.tot, list(x.star = .N, sales.star = sum(sales)), keyby = "cust"]
    cbs <- merge(cbs, val, all.x = TRUE, by = "cust")
    cbs[is.na(x.star), `:=`(x.star, 0)]
    cbs[is.na(sales.star), `:=`(sales.star, 0)]
    setcolorder(cbs, c("cust", "x", "t.x", "litt", "sales", "sales.x", "first", "T.cal",
                       "T.star", "x.star", "sales.star"))
  } else {
    setcolorder(cbs, c("cust", "x", "t.x", "litt", "sales", "sales.x", "first", "T.cal"))
  }
  # return same object type as was passed
  if (!has.sales) {
    elog_dt[, `:=`(sales, NULL)]
    cbs[, `:=`(sales, NULL)]
    cbs[, `:=`(sales.x, NULL)]
    if (has.holdout) cbs[, `:=`(sales.star, NULL)]
  }
  if (!is.dt) {
    cbs <- data.frame(cbs)
  }
  return(cbs)
}


#' Convert Event Log to Transaction Counts
#'
#' Aggregates an event log to either incremental or cumulative number of
#' transactions. If \code{first=TRUE} then the initial transactions of each
#' customer are included in the count as well.
#'
#' Duplicate transactions with identical \code{cust} and \code{date} (or
#' \code{t}) field are counted only once.
#'
#' @param elog Event log, a \code{data.frame} with columns \code{cust} and
#'   transaction time \code{t} or \code{date}.
#' @param by Only return every \code{by}-th count Defaults to 7, and thus
#'   returns weekly numbers.
#' @param first If TRUE, then the first transaction for each customer is being
#'   counted as well
#' @return Numeric vector of transaction counts.
#' @export
#' @examples
#' data("groceryElog")
#' cum <- elog2cum(groceryElog)
#' plot(cum, typ="l", frame = FALSE)
#' inc <- elog2inc(groceryElog)
#' plot(inc, typ="l", frame = FALSE)
elog2cum <- function(elog, by = 7, first = FALSE) {
  t0 <- N <- cust <- NULL  # suppress checkUsage warnings
  stopifnot("cust" %in% names(elog))
  stopifnot(is.logical(first) & length(first) == 1)
  is.dt <- is.data.table(elog)
  if (!is.dt) {
    elog <- as.data.table(elog)
  } else {
    elog <- copy(elog)
  }
  if (!"t" %in% names(elog)) {
    stopifnot("date" %in% names(elog))
    cohort_start <- min(as.numeric(elog$date))
    elog[, `:=`(t, as.numeric(date) - cohort_start)]
  }
  elog <- unique(elog[, list(cust, t)])
  elog[, `:=`(t0, min(t)), by = "cust"]
  grid <- data.table(t = 0 : ceiling(max(elog$t)))
  grid <- merge(grid, elog[first | t > t0, .N, keyby = list(t = ceiling(t))], all.x = TRUE, by = "t")
  grid <- grid[is.na(N), N := 0L]
  cum <- cumsum(grid$N)
  cum <- cum[seq(by, length(cum), by = by)]
  return(cum)
}


#' @rdname elog2cum
#' @export
elog2inc <- function(elog, by = 7, first = FALSE) {
  cum <- elog2cum(elog = elog, by = by, first = first)
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
  BTYD::dc.check.model.params(printnames, params, func)
  # then check for names, if these are present
  if (!is.null(names(params))) {
    idx <- names(params) != ""
    if (any(printnames[idx] != names(params)[idx])) {
      stop("Parameter names do not match - ", paste0(printnames, collapse = ","), " != ",
        paste0(names(params), collapse = ","), call. = FALSE)
    }
  }
}


#' Generic Method for Tracking Plots
#'
#' @keywords internal
dc.PlotTracking <- function(actual, expected, T.cal = NULL,
                            xlab = "", ylab = "", title = "",
                            xticklab = NULL, ymax = NULL,
                            legend = c("Actual", "Model")) {

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
  if (!is.null(legend) & length(legend) == 2) {
    pos <- ifelse(which.max(expected) == length(expected), "bottomright", "topright")
    legend(pos, legend = legend, col = 1:2, lty = 1:2, lwd = 1)
  }
  return(rbind(actual, expected))
}


#' Generic Method for Plotting Frequency vs. Conditional Expected Frequency
#'
#' @keywords internal
dc.PlotFreqVsConditionalExpectedFrequency <- function(x, actual, expected, censor,
                                                      xlab, ylab, xticklab, title) {

  bin <- bin.size <- transaction.actual <- transaction.expected <- N <- NULL # suppress checkUsage warnings
  if (length(x) != length(actual) | length(x) != length(expected) |
      !is.numeric(x) | !is.numeric(actual) | !is.numeric(expected) |
      any(x < 0) | any(actual < 0) | any(expected < 0))
    stop("x, actual and expected must be non-negative numeric vectors of same length.")

  if (censor > max(x)) censor <- max(x)
  dt <- data.table(x, actual, expected)
  dt[, bin := pmin(x, censor)]
  st <- dt[, list(transaction.actual = mean(actual),
                  transaction.expected = mean(expected),
                  bin.size = .N), keyby = bin]
  st <- merge(data.table(bin = 0:censor), st, all.x = TRUE, by = "bin")
  comparison <- t(st)[-1, ]
  col.names <- paste(rep("freq", length(censor + 1)), (0:censor), sep = ".")
  col.names[censor + 1] <- paste0(col.names[censor + 1], "+")
  colnames(comparison) <- col.names
  if (is.null(xticklab) == FALSE) {
    x.labels <- xticklab
  } else {
    if (censor < ncol(comparison)) {
      x.labels <- 0:(censor)
      x.labels[censor + 1] <- paste0(censor, "+")
    } else {
      x.labels <- 0:(ncol(comparison))
    }
  }
  actual <- comparison[1, ]
  expected <- comparison[2, ]
  ylim <- c(0, ceiling(max(c(actual, expected)) * 1.1))
  plot(actual, type = "l", xaxt = "n", col = 1, ylim = ylim, xlab = xlab, ylab = ylab, main = title)
  lines(expected, lty = 2, col = 2)
  axis(1, at = 1:ncol(comparison), labels = x.labels)
  legend("topleft", legend = c("Actual", "Model"), col = 1:2, lty = 1:2, lwd = 1)
  return(comparison)
}


#' Generic Method for Plotting Frequency vs. Conditional Expected Frequency
#'
#' @keywords internal
dc.PlotRecVsConditionalExpectedFrequency <- function(t.x, actual, expected,
                                                     xlab, ylab, xticklab, title) {
  bin <- bin.size <- N <- NULL # suppress checkUsage warnings
  if (length(t.x) != length(actual) | length(t.x) != length(expected) |
      !is.numeric(t.x) | !is.numeric(actual) | !is.numeric(expected) |
      any(t.x < 0) | any(actual < 0) | any(expected < 0))
    stop("t.x, actual and expected must be non-negative numeric vectors of same length.")

  dt <- data.table(t.x, actual, expected)
  dt[, bin := floor(t.x)]
  st <- dt[, list(actual = mean(actual),
                  expected = mean(expected),
                  bin.size = .N), keyby = bin]
  st <- merge(data.table(bin = 1:floor(max(t.x))), st[bin > 0], all.x = TRUE, by = "bin")
  comparison <- t(st)[-1, ]
  x.labels <- NULL
  if (is.null(xticklab) == FALSE) {
    x.labels <- xticklab
  } else {
    x.labels <- 1:max(st$bin)
  }
  actual <- comparison[1, ]
  expected <- comparison[2, ]
  ylim <- c(0, ceiling(max(c(actual, expected), na.rm = TRUE) * 1.1))
  plot(actual, type = "l", xaxt = "n", col = 1, ylim = ylim, xlab = xlab, ylab = ylab, main = title)
  lines(expected, lty = 2, col = 2)
  axis(1, at = 1:ncol(comparison), labels = x.labels)
  legend("topleft", legend = c("Actual", "Model"), col = 1:2, lty = 1:2, lwd = 1)
  return(comparison)
}
