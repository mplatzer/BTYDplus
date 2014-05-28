
#' Parameter Estimation for the NBD model
#' 
#' Estimates parameters for the NBD model via Maximum Likelihood Estimation.
#' 
#' @param cal.cbs calibration period CBS. It must contain columns for frequency 
#'   'x' and total time observed 'T.cal'. Optionally a column 'custs' can be 
#'   provided, which represents number of customers with a specific combination 
#'   of frequency 'x' and 'T.cal'.
#' @param par.start initial NBD parameters - a vector with 'r' and 'alpha' in 
#'   that order.
#' @param max.param.value the upper bound on parameters
#' @return list of estimated parameters
#' @import BTYD
#' @export
#' @references EHRENBERG, ASC. "The Pattern of Consumer Purchases." Quantitative
#'   techniques in marketing analysis: text and readings (1962): 355.
#' @examples
#' # generate artificial data 
# ' params <- list(r=0.55, alpha=10.56)
# ' cbs <- nbd.GenerateData(n=1000, T.cal=32, T.star=32, params=params)$cbs
# ' 
# ' # estimate parameters, and compare to true parameters
# ' est <- nbd.EstimateParameters(cbs)
# ' rbind(params, est=round(est, 3))
# ' #           r   alpha 
# ' # params 0.55  10.56 
# ' # est    0.531 10.444
# ' 
# ' # estimate future transactions, and compare to actuals from holdout period
# ' cbs$x.est <- nbd.ConditionalExpectedTransactions(est, cbs$T.star, cbs$x, cbs$T.cal)
# ' sqrt(mean((cbs$x.star-cbs$x.est)^2))
# ' # 1.752681
nbd.EstimateParameters <- function(cal.cbs, par.start = c(1, 1), max.param.value = 10000) {
  dc.check.model.params(c("r", "alpha"), par.start, 
    "nbd.EstimateParameters")
  nbd.eLL <- function(params, cal.cbs, max.param.value) {
    params <- exp(params)
    params[params > max.param.value] <- max.param.value
    return(-1 * nbd.cbs.LL(params, cal.cbs))
  }
  logparams <- log(par.start)
  results <- optim(logparams, nbd.eLL, cal.cbs = cal.cbs, 
    max.param.value = max.param.value, method = "L-BFGS-B")
  estimated.params <- exp(results$par)
  estimated.params[estimated.params > max.param.value] <- max.param.value
  return(estimated.params)
}


#' Calculate the log-likelihood of the NBD model
#' 
#' @param params NBD parameters - a vector with r and alpha, in that
#'   order.
#' @param cal.cbs calibration period CBS. It must contain columns for frequency 
#'   'x' and total time observed 'T.cal'. Optionally a column 'custs' can be 
#'   provided, which represents number of customers with a specific combination 
#'   of frequency 'x' and 'T.cal' data.frame with columns 'x' and 'T.cal' and 
#'   optional 'custs'.
#' @return the total log-likelihood for the provided data.
#' @export
#' @seealso nbd.EstimateParameters
nbd.cbs.LL <- function(params, cal.cbs) {
  dc.check.model.params(c("r", "alpha"), params, 
                        "nbd.cbs.LL")  
  tryCatch(x <- cal.cbs[, "x"], error = function(e) stop("Error in nbd.cbs.LL: cal.cbs must have a frequency column labelled \"x\""))
  tryCatch(T.cal <- cal.cbs[, "T.cal"], error = function(e) stop("Error in nbd.cbs.LL: cal.cbs must have a column for length of time observed labelled \"T.cal\""))
  if ("custs" %in% colnames(cal.cbs)) {
    custs <- cal.cbs[, "custs"]
  } else {
    custs <- rep(1, length(x))
  }
  return(sum(custs * nbd.LL(params, x, T.cal)))
}


#' Calculate the log-likelihood of the NBD model
#' 
#' @param params NBD parameters - a vector with r and alpha, in that
#'   order.
#' @param cal.cbs calibration period CBS. It must contain columns for frequency 
#'   'x' and total time observed 'T.cal'. Optionally a column 'custs' can be 
#'   provided, which represents number of customers with a specific combination 
#'   of frequency 'x' and 'T.cal' data.frame with columns 'x' and 'T.cal' and 
#'   optional 'custs'.
#' @return a vector of log-likelihoods
#' @export
#' @seealso nbd.EstimateParameters
nbd.LL <- function(params, x, T.cal) {
  max.length <- max(length(x), length(T.cal))
  if (max.length%%length(x)) 
    warning("Maximum vector length not a multiple of the length of x")
  if (max.length%%length(T.cal)) 
    warning("Maximum vector length not a multiple of the length of T.cal")
  dc.check.model.params(c("r", "alpha"), params, 
                        "nbd.LL")  
  if (any(x < 0) || !is.numeric(x)) 
    stop("x must be numeric and may not contain negative numbers.")
  if (any(T.cal < 0) || !is.numeric(T.cal)) 
    stop("T.cal must be numeric and may not contain negative numbers.")
  x <- rep(x, length.out = max.length)
  T.cal <- rep(T.cal, length.out = max.length)
  r <- params[1]
  alpha <- params[2]
  P1 <- lgamma(r+x) + r*log(alpha)
  P2 <- lgamma(r) + (r+x)*log(alpha+T.cal)
  llh <- P1 - P2
  return(llh)
}


#' NBD Conditional Expected Transactions
#' 
#' Uses NBD model parameters and a customer's past transaction behavior to
#' return the number of transactions they are expected to make in a given time
#' period.
#' 
#' @param params NBD parameters - a vector with r and alpha, in that order.
#' @param T.star length of time for which we are calculating the expected number
#'   of transactions.
#' @param x number of repeat transactions in the calibration period T.cal, or a
#'   vector of calibration period frequencies.
#' @param T.cal length of calibration period, or a vector of calibration period
#'   lengths.
#' @return Number of transactions a customer is expected to make in a time
#'   period of length t, conditional on their past behavior. If any of the input
#'   parameters has a length greater than 1, this will be a vector of expected
#'   number of transactions.
#' @export
#' @seealso nbd.EstimateParameters
nbd.ConditionalExpectedTransactions <- function(params, T.star, x, T.cal) {
  max.length <- max(length(T.star), length(x), length(T.cal))
  if (max.length%%length(T.star)) 
    warning("Maximum vector length not a multiple of the length of T.star")
  if (max.length%%length(x)) 
    warning("Maximum vector length not a multiple of the length of x")
  if (max.length%%length(T.cal)) 
    warning("Maximum vector length not a multiple of the length of T.cal")
  dc.check.model.params(c("r", "alpha"), params, 
    "nbd.ConditionalExpectedTransactions")
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
  return (T.star * (r+x) / (alpha+T.cal))
}


#' Simulate data according to NBD model assumptions
#' 
#' @param n number of customers
#' @param T.cal length of calibration period
#' @param T.star length of holdout period
#' @param params NBD parameters - a list with 'r' and 'alpha'
#' @param return.elog boolean - if TRUE then the event log is returned in
#'   addition to the CBS summary
#' @return list with elements 'cbs' and 'elog' containing data.frames
#' @export
#' @seealso nbd.EstimateParameters
nbd.GenerateData <- function(n, T.cal, T.star, params, return.elog=F) {
  # sample intertransaction timings parameter lambda for each customer
  lambdas <- rgamma(n, shape=params$r, rate=params$alpha)
  
  # sample intertransaction timings, and generate CBS data frame
  cbs <- data.frame()
  for (i in 1:n) {
    lambda <- lambdas[i]
    # sample transaction times
    times <- cumsum(c(0, rexp(10*(T.cal+T.star)*lambda, rate=lambda)))
    if (return.elog) elog <- rbind(elog, data.frame(cust=i, t=times[times<(T.cal+T.star)]))
    # determine frequency, recency, etc.
    ts.cal <- times[times<T.cal]
    ts.star <- times[times>=T.cal & times<(T.cal+T.star)]
    cbs[i, "x"] <- length(ts.cal)-1
    cbs[i, "t.x"] <- max(ts.cal)
    cbs[i, "T.cal"] <- T.cal
    cbs[i, "x.star"] <- length(ts.star)
    cbs[i, "T.star"] <- T.star
    cbs[i, "lambda"] <- lambda
  }
  out <- list(cbs=cbs)
  if (return.elog) {
    out$elog <- transform(elog, date=as.Date("2001/01/01") + t)
  }
  return(out)
}
