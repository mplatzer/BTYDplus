
#' Parameter Estimation for the NBD model
#'
#' Estimates parameters for the NBD model via Maximum Likelihood Estimation.
#'
#' @param cal.cbs Calibration period CBS. It must contain columns for frequency
#'   \code{x} and total time observed \code{T.cal}.
#' @param par.start Initial NBD parameters - a vector with \code{r} and \code{alpha} in
#'   that order.
#' @param max.param.value Upper bound on parameters.
#' @return List of estimated parameters.
#' @export
#' @references EHRENBERG, ASC. 'The Pattern of Consumer Purchases.' Quantitative
#'   techniques in marketing analysis: text and readings (1962): 355.
#' @examples
#' data("groceryElog")
#' cbs <- elog2cbs(groceryElog)
#' nbd.EstimateParameters(cbs)
nbd.EstimateParameters <- function(cal.cbs, par.start = c(1, 1), max.param.value = 10000) {
  dc.check.model.params.safe(c("r", "alpha"), par.start, "nbd.EstimateParameters")
  nbd.eLL <- function(params, cal.cbs, max.param.value) {
    params <- exp(params)
    params[params > max.param.value] <- max.param.value
    return(-1 * nbd.cbs.LL(params, cal.cbs))
  }
  logparams <- log(par.start)
  results <- optim(logparams, nbd.eLL, cal.cbs = cal.cbs, max.param.value = max.param.value, method = "L-BFGS-B")
  estimated.params <- exp(results$par)
  estimated.params[estimated.params > max.param.value] <- max.param.value
  names(estimated.params) <- c("r", "alpha")
  return(estimated.params)
}


#' Calculate the log-likelihood of the NBD model
#'
#' @param params NBD parameters - a vector with r and alpha, in that
#'   order.
#' @param cal.cbs Calibration period CBS. It must contain columns for frequency
#'   \code{x} and total time observed \code{T.cal}.
#' @return The total log-likelihood for the provided data.
#' @export
#' @examples
#' data("groceryElog")
#' cbs <- elog2cbs(groceryElog)
#' params <- nbd.EstimateParameters(cbs)
#' nbd.cbs.LL(params, cbs)
nbd.cbs.LL <- function(params, cal.cbs) {
  dc.check.model.params.safe(c("r", "alpha"), params, "nbd.cbs.LL")
  tryCatch(x <- cal.cbs$x,
    error = function(e) stop("cal.cbs must have a frequency column labelled \"x\""))
  tryCatch(T.cal <- cal.cbs$T.cal,
    error = function(e) stop("cal.cbs must have a column for length of time observed labelled \"T.cal\""))
  ll <- nbd.LL(params, x, T.cal)
  return(sum(ll))
}


#' Calculate the log-likelihood of the NBD model
#'
#' @param params NBD parameters - a vector with \code{r} and \code{alpha}, in that
#'   order.
#' @param x Frequency, i.e. number of re-purchases.
#' @param T.cal Total time of observation period.
#' @return A numeric vector of log-likelihoods.
#' @export
#' @seealso \code{\link{nbd.cbs.LL}}
nbd.LL <- function(params, x, T.cal) {
  max.length <- max(length(x), length(T.cal))
  if (max.length %% length(x))
    warning("Maximum vector length not a multiple of the length of x")
  if (max.length %% length(T.cal))
    warning("Maximum vector length not a multiple of the length of T.cal")
  dc.check.model.params.safe(c("r", "alpha"), params, "nbd.LL")
  if (any(x < 0) || !is.numeric(x))
    stop("x must be numeric and may not contain negative numbers.")
  if (any(T.cal < 0) || !is.numeric(T.cal))
    stop("T.cal must be numeric and may not contain negative numbers.")
  x <- rep(x, length.out = max.length)
  T.cal <- rep(T.cal, length.out = max.length)
  r <- params[1]
  alpha <- params[2]
  P1 <- lgamma(r + x) + r * log(alpha)
  P2 <- lgamma(r) + (r + x) * log(alpha + T.cal)
  llh <- P1 - P2
  return(llh)
}


#' NBD Conditional Expected Transactions
#'
#' Uses NBD model parameters and a customer's past transaction behavior to
#' return the number of transactions they are expected to make in a given time
#' period.
#'
#' @param params NBD parameters - a vector with \code{r} and \code{alpha}, in that order.
#' @param T.star Length of time for which we are calculating the expected number
#'   of transactions.
#' @param x Number of repeat transactions in the calibration period \code{T.cal}, or a
#'   vector of calibration period frequencies.
#' @param T.cal Length of calibration period, or a vector of calibration period
#'   lengths.
#' @return Number of transactions a customer is expected to make in a time
#'   period of length t, conditional on their past behavior. If any of the input
#'   parameters has a length greater than 1, this will be a vector of expected
#'   number of transactions.
#' @export
#' @examples
#' data("groceryElog")
#' cbs <- elog2cbs(groceryElog, T.cal = "2006-12-31")
#' params <- nbd.EstimateParameters(cbs)
#' xstar.est <- nbd.ConditionalExpectedTransactions(params, cbs$T.star, cbs$x, cbs$T.cal)
#' sum(xstar.est) # expected total number of transactions during holdout
nbd.ConditionalExpectedTransactions <- function(params, T.star, x, T.cal) {
  max.length <- max(length(T.star), length(x), length(T.cal))
  if (max.length %% length(T.star))
    warning("Maximum vector length not a multiple of the length of T.star")
  if (max.length %% length(x))
    warning("Maximum vector length not a multiple of the length of x")
  if (max.length %% length(T.cal))
    warning("Maximum vector length not a multiple of the length of T.cal")
  dc.check.model.params.safe(c("r", "alpha"), params, "nbd.ConditionalExpectedTransactions")
  if (any(T.star < 0) || !is.numeric(T.star))
    stop("T.star must be numeric and may not contain negative numbers.")
  if (any(x < 0) || !is.numeric(x))
    stop("x must be numeric and may not contain negative numbers.")
  if (any(T.cal < 0) || !is.numeric(T.cal))
    stop("T.cal must be numeric and may not contain negative numbers.")
  T.star <- rep(T.star, length.out = max.length)
  x <- rep(x, length.out = max.length)
  T.cal <- rep(T.cal, length.out = max.length)
  r <- params[1]
  alpha <- params[2]
  return(unname(T.star * (r + x) / (alpha + T.cal)))
}


#' Simulate data according to NBD model assumptions
#'
#' @param n Number of customers.
#' @param T.cal Length of calibration period.
#' @param T.star Length of holdout period. This may be a vector.
#' @param params NBD parameters - a vector with \code{r} and \code{alpha} in that order.
#' @param date.zero Initial date for cohort start. Can be of class character, Date or POSIXt.
#' @return List of length 2:
#' \item{\code{cbs}}{A data.frame with a row for each customer and the summary statistic as columns.}
#' \item{\code{elog}}{A data.frame with a row for each transaction, and columns \code{cust}, \code{date} and \code{t}.}
#' @export
#' @examples
#' n <- 1000  # no. of customers
#' T.cal <- 32  # length of calibration period
#' T.star <- 32  # length of hold-out period
#' params <- c(r = 0.85, alpha = 4.45)  # purchase frequency lambda_i ~ Gamma(r, alpha)
#' data <- nbd.GenerateData(n, T.cal, T.star, params)
#' cbs <- data$cbs  # customer by sufficient summary statistic - one row per customer
#' elog <- data$elog  # Event log - one row per event/purchase
nbd.GenerateData <- function(n, T.cal, T.star, params, date.zero = "2000-01-01") {
  # check model parameters
  dc.check.model.params.safe(c("r", "alpha"), params, "nbd.GenerateData")

  # set start date for each customer, so that they share same T.cal date
  T.cal.fix <- max(T.cal)
  T.cal <- rep(T.cal, length.out = n)
  T.zero <- T.cal.fix - T.cal
  date.zero <- as.POSIXct(date.zero)

  r <- params[1]
  alpha <- params[2]

  # sample intertransaction timings parameter lambda for each customer
  lambdas <- rgamma(n, shape = r, rate = alpha)

  # sample intertransaction timings
  elog_list <- lapply(1:n, function(i) {
    itts <- rexp(10 * (T.cal[i] + max(T.star)) * lambdas[i], rate = lambdas[i])
    ts <- cumsum(c(0, itts))
    ts <- T.zero[i] + ts # shift by T_0
    ts <- ts[ts <= (T.cal.fix + max(T.star))] # trim to observation length
    return(ts)
  })

  # build elog
  elog <- data.table("cust" = rep(1:n, sapply(elog_list, length)), "t" = unlist(elog_list))
  elog[["date"]] <- date.zero + elog[["t"]] * 3600 * 24 * 7

  # build cbs
  date.cal <- date.zero + T.cal.fix * 3600 * 24 * 7
  date.tot <- date.cal + T.star * 3600 * 24 * 7
  cbs <- elog2cbs(elog, T.cal = date.cal)
  if (length(T.star) == 1) set(cbs, j = "T.star", value = T.star[1])
  xstar.cols <- if (length(T.star) == 1) "x.star" else paste0("x.star", T.star)
  for (j in 1:length(date.tot)) {
    set(cbs, j = xstar.cols[j],
        value = sapply(elog_list, function(t) sum(t > T.cal.fix & t <= T.cal.fix + T.star[j])))
  }
  set(cbs, j = "lambda", value = lambdas)

  return(list(cbs = setDF(cbs), elog = setDF(elog)))
}
