
#' Parameter Estimation for the BG/NBD model
#' 
#' Estimates parameters for the BG/NBD model via Maximum Likelihood Estimation.
#' 
#' @param cal.cbs calibration period CBS. It must contain columns for frequency 
#'   \code{x}, for recency \code{t.x} and total time observed \code{T.cal}. Optionally a 
#'   column \code{custs} can be provided, which represents number of customers with a
#'   specific combination of frequency \code{x}, recency \code{t.x} and \code{T.cal}.
#' @param par.start initial BG/NBD parameters - a vector with \code{r}, \code{alpha}, \code{a} and
#'   \code{b} in that order.
#' @param max.param.value the upper bound on parameters
#' @return list of estimated parameters
#' @import BTYD
#' @export
#' @references Fader, Peter S., Bruce GS Hardie, and Ka Lok Lee. "Counting Your Customers the Easy Way: An Alternative to the Pareto/NBD Model." Marketing Science 24.2 (2005): 275-284.
#' @example demo/bg-nbd.r
bgnbd.EstimateParameters <- function(cal.cbs, par.start = c(1, 1, 1, 1), max.param.value = 10000) {
  dc.check.model.params(c("r", "alpha", "a", "b"), par.start, 
    "bgnbd.EstimateParameters")
  bgnbd.eLL <- function(params, cal.cbs, max.param.value) {
    params <- exp(params)
    params[params > max.param.value] <- max.param.value
    return(-1 * bgnbd.cbs.LL(params, cal.cbs))
  }
  logparams <- log(par.start)
  results <- optim(logparams, bgnbd.eLL, cal.cbs = cal.cbs, 
    max.param.value = max.param.value, method = "L-BFGS-B")
  estimated.params <- exp(results$par)
  estimated.params[estimated.params > max.param.value] <- max.param.value
  names(estimated.params) <- c("r", "alpha", "a", "b")
  return(estimated.params)
}


#' Calculate the log-likelihood of the BG/NBD model
#' 
#' @param params BG/NBD parameters - a vector with \code{r}, \code{alpha}, \code{a} and
#'   \code{b} in that order.
#' @param cal.cbs calibration period CBS. It must contain columns for frequency 
#'   \code{x}, for recency \code{t.x} and total time observed \code{T.cal}. Optionally a 
#'   column \code{custs} can be provided, which represents number of customers with a
#'   specific combination of frequency \code{x} and \code{T.cal}.
#' @return the total log-likelihood for the provided data.
#' @export
#' @seealso \code{\link{bgnbd.EstimateParameters}}
bgnbd.cbs.LL <- function(params, cal.cbs) {
  dc.check.model.params(c("r", "alpha", "a", "b"), params, 
   "bgnbd.cbs.LL")  
  tryCatch(x <- cal.cbs[, "x"], error = function(e) stop("Error in bgnbd.cbs.LL: cal.cbs must have a frequency column labelled \"x\""))
  tryCatch(t.x <- cal.cbs[, "t.x"], error = function(e) stop("Error in bgnbd.cbs.LL: cal.cbs must have a recency column labelled \"t.x\""))
  tryCatch(T.cal <- cal.cbs[, "T.cal"], error = function(e) stop("Error in bgnbd.cbs.LL: cal.cbs must have a column for length of time observed labelled \"T.cal\""))
  if ("custs" %in% colnames(cal.cbs)) {
    custs <- cal.cbs[, "custs"]
  } else {
    custs <- rep(1, length(x))
  }
  return(sum(custs * bgnbd.LL(params, x, t.x, T.cal)))
}


#' Calculate the log-likelihood of the BG/NBD model
#' 
#' @param params BG/NBD parameters - a vector with \code{r}, \code{alpha}, \code{a} and
#'   \code{b} in that order.
#' @param x frequency, i.e. number of re-purchases
#' @param t.x recency, i.e. time elapsed from first purchase to last purchase
#' @param T.cal total time of observation period
#' @return a vector of log-likelihoods
#' @export
#' @seealso \code{\link{bgnbd.EstimateParameters}}
bgnbd.LL <- function(params, x, t.x, T.cal) {
  max.length <- max(length(x), length(t.x), length(T.cal))
  if (max.length%%length(x)) 
    warning("Maximum vector length not a multiple of the length of x")
  if (max.length%%length(t.x)) 
    warning("Maximum vector length not a multiple of the length of t.x")
  if (max.length%%length(T.cal)) 
    warning("Maximum vector length not a multiple of the length of T.cal")
  dc.check.model.params(c("r", "alpha", "a", "b"), params, 
   "bgnbd.LL")  
  if (any(x < 0) || !is.numeric(x)) 
    stop("x must be numeric and may not contain negative numbers.")
  if (any(t.x < 0) || !is.numeric(t.x)) 
    stop("t.x must be numeric and may not contain negative numbers.")
  if (any(T.cal < 0) || !is.numeric(T.cal)) 
    stop("T.cal must be numeric and may not contain negative numbers.")
  x <- rep(x, length.out = max.length)
  t.x <- rep(t.x, length.out = max.length)
  T.cal <- rep(T.cal, length.out = max.length)
  r <- params[1]
  alpha <- params[2]
  a <- params[3]
  b <- params[4]
  P1 <- lgamma(r+x) + r * log(alpha) - lgamma(r)
  P2 <- lgamma(a+b) + lgamma(b+x) - lgamma(b) - lgamma(a+b+x)
  P3 <- log(1) # x * log(t.x) - lgamma(x+1)
  P4 <- (1/(alpha + T.cal))^(r+x)
  P5 <- ifelse(x>0, (a/(b+x-1)) * (1/(alpha + t.x))^(r+x), 0)
  llh <- P1 + P2 + P3 + log(P4 + P5)
  return(llh)
}


#' BG/NBD P(alive)
#' 
#' Uses BG/NBD model parameters and a customer's past transaction behavior to
#' return the probability that they are still alive at the end of the
#' calibration period.
#' 
#' @param params BG/NBD parameters - a vector with \code{r}, \code{alpha}, \code{a} and \code{b} in
#'   that order.
#' @param x number of repeat transactions in the calibration period T.cal, or a 
#'   vector of calibration period frequencies.
#' @param t.x recency, i.e. length between first and last transaction during
#'   calibration period.
#' @param T.cal length of calibration period, or a vector of calibration period 
#'   lengths.
#' @return Probability that the customer is still alive at the end of the calibration period.
#' @export
#' @seealso \code{\link{bgnbd.EstimateParameters}}
#' @example demo/bg-nbd.r
bgnbd.PAlive <- function(params, x, t.x, T.cal) {
  max.length <- max(length(x), length(t.x), length(T.cal))
  if (max.length%%length(x)) 
    warning("Maximum vector length not a multiple of the length of x")
  if (max.length%%length(t.x)) 
    warning("Maximum vector length not a multiple of the length of t.x")
  if (max.length%%length(T.cal)) 
    warning("Maximum vector length not a multiple of the length of T.cal")
  dc.check.model.params(c("r", "alpha", "a", "b"), params, "bgnbd.PAlive")
  if (any(x < 0) || !is.numeric(x)) 
    stop("x must be numeric and may not contain negative numbers.")
  if (any(t.x < 0) || !is.numeric(t.x)) 
    stop("t.x must be numeric and may not contain negative numbers.")
  if (any(T.cal < 0) || !is.numeric(T.cal)) 
    stop("T.cal must be numeric and may not contain negative numbers.")
  x <- rep(x, length.out = max.length)
  t.x <- rep(t.x, length.out = max.length)
  T.cal <- rep(T.cal, length.out = max.length)
  r <- params[1]
  alpha <- params[2]
  a <- params[3]
  b <- params[4]  
  P1 <- ifelse(x>0, (a / (b + x - 1)) * ((alpha + T.cal) / (alpha + t.x))^(r + x), 0)
  return (1 / (1 + P1))
}


#' BG/NBD Conditional Expected Transactions
#' 
#' Uses BG/NBD model parameters and a customer's past transaction behavior to 
#' return the number of transactions they are expected to make in a given time 
#' period.
#' 
#' @param params BG/NBD parameters - a vector with \code{r}, \code{alpha}, \code{a} and \code{b} in
#'   that order.
#' @param T.star length of time for which we are calculating the expected number
#'   of transactions.
#' @param x number of repeat transactions in the calibration period T.cal, or a 
#'   vector of calibration period frequencies.
#' @param t.x recency, i.e. length between first and last transaction during
#'   calibration period.
#' @param T.cal length of calibration period, or a vector of calibration period 
#'   lengths.
#' @return Number of transactions a customer is expected to make in a time 
#'   period of length t, conditional on their past behavior. If any of the input
#'   parameters has a length greater than 1, this will be a vector of expected 
#'   number of transactions.
#' @import gsl
#' @export
#' @seealso \code{\link{bgnbd.EstimateParameters}}
#' @example demo/bg-nbd.r
bgnbd.ConditionalExpectedTransactions <- function(params, T.star, x, t.x, T.cal) {
  max.length <- max(length(T.star), length(x), length(t.x), 
    length(T.cal))
  if (max.length%%length(T.star)) 
    warning("Maximum vector length not a multiple of the length of T.star")
  if (max.length%%length(x)) 
    warning("Maximum vector length not a multiple of the length of x")
  if (max.length%%length(t.x)) 
    warning("Maximum vector length not a multiple of the length of t.x")
  if (max.length%%length(T.cal)) 
    warning("Maximum vector length not a multiple of the length of T.cal")
  dc.check.model.params(c("r", "alpha", "a", "b"), params, 
    "bgnbd.ConditionalExpectedTransactions")
  if (any(T.star < 0) || !is.numeric(T.star)) 
    stop("T.star must be numeric and may not contain negative numbers.")
  if (any(x < 0) || !is.numeric(x)) 
    stop("x must be numeric and may not contain negative numbers.")
  if (any(t.x < 0) || !is.numeric(t.x)) 
    stop("t.x must be numeric and may not contain negative numbers.")
  if (any(T.cal < 0) || !is.numeric(T.cal)) 
    stop("T.cal must be numeric and may not contain negative numbers.")
  T.star <- rep(T.star, length.out = max.length)
  x <- rep(x, length.out = max.length)
  t.x <- rep(t.x, length.out = max.length)
  T.cal <- rep(T.cal, length.out = max.length)
  r <- params[1]
  alpha <- params[2]
  a <- params[3]
  b <- params[4]
  P1 <- ((a + b + x - 1)/(a - 1))
  P2 <- 1 - ((alpha + T.cal) / (alpha + T.cal + T.star))^(r + x)*hyperg_2F1(r + x, b + x, a + b + x - 1, T.star/(alpha + T.cal + T.star))
  P3 <- bgnbd.PAlive(params, x, t.x, T.cal)
  return(P1 * P2 * P3)
}


#' Simulate data according to BG/NBD model assumptions
#' 
#' @param n number of customers
#' @param T.cal length of calibration period
#' @param T.star length of holdout period
#' @param params BG/NBD parameters - a vector with \code{r}, \code{alpha}, \code{a} and
#'   \code{b} in that order.
#' @param return.elog boolean - if \code{TRUE} then the event log is returned in
#'   addition to the CBS summary
#' @return list with elements \code{cbs} and \code{elog} containing data.frames
#' @export
#' @seealso \code{\link{bgnbd.EstimateParameters}}
#' @example demo/bg-nbd.r
bgnbd.GenerateData <- function(n, T.cal, T.star, params, return.elog=F) {
  # check model parameters
  dc.check.model.params(c("r", "alpha", "a", "b"), params,
                        "bgnbd.GenerateData")
  
  r <- params[1]
  alpha <- params[2]
  a <- params[3]
  b <- params[4]  
  
  # sample intertransaction timings parameter lambda for each customer
  lambdas <- rgamma(n, shape=r, rate=alpha)

  # sample dropout probability p for each customer
  ps <- rbeta(n, a, b)
  
  # sample intertransaction timings & churn
  cbs <- data.frame()
  elog <- data.frame(cust=numeric(0), t=numeric(0))
  for (i in 1:n) {
    p <- ps[i]
    lambda <- lambdas[i]
    # sample no. of transactions until churn
    churn <- which.max(rbinom(10/p, 1, p))
    # sample transaction times
    times <- cumsum(c(0, rexp(churn, rate=lambda)))
    if (return.elog) elog <- rbind(elog, data.frame(cust=i, t=times[times<(T.cal+T.star)]))
    # determine frequency, recency, etc.
    ts.cal <- times[times<T.cal]
    ts.star <- times[times>=T.cal & times<(T.cal+T.star)]
    cbs[i, "x"] <- length(ts.cal)-1
    cbs[i, "t.x"] <- max(ts.cal)
    cbs[i, "T.cal"] <- T.cal
    cbs[i, "alive"] <- churn>length(ts.cal)
    cbs[i, "x.star"] <- length(ts.star)
    cbs[i, "p"] <- p
    cbs[i, "lambda"] <- lambda
    cbs[i, "churn"] <- churn
    cbs[i, "T.star"] <- T.star
  }
  out <- list(cbs=cbs)
  if (return.elog) out$elog <- elog
  return(out)
}
