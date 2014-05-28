
#' Estimates parameters for the Gamma/Gompertz/NBD model (Bemmaor, Glady, 2012) 
#' via Maximum Likelihood Estimation.
#' 
#' @param cal.cbs calibration period CBS. It must contain columns for frequency 
#'   'x', for recency 't.x.' and total time observed 'T.cal'. Optionally a 
#'   column 'custs' can be provided, which represents number of customers with a
#'   specific combination of frequency 'x', recency 't.x' and 'T.cal'.
#' @param par.start initial Gamma/Gompertz/NBD parameters - a vector with 'r',
#'   'alpha', 'b', 's' and 'beta' in that order.
#' @param max.param.value the upper bound on parameters
#' @return list of estimated parameters
#' @export
ggnbd.EstimateParameters <- function(cal.cbs, par.start = c(1, 1, .5, 1, 1), min.param.value=1e-6, max.param.value=1e+4, trace=0) {
  dc.check.model.params(c("r", "alpha", "b", "s", "beta"), par.start, 
    "ggnbd.EstimateParameters")
  count <- 0
  ggnbd.eLL <- function(params, cal.cbs) {
    params <- exp(params)
    loglik <- ggnbd.cbs.LL(params, cal.cbs)
    count <<- count + 1
    if (trace>0 & count%%trace==0)
      cat("iter", count, ":", sprintf("%12.2f", loglik), ":", sprintf("%10.6f", params), "\n")
    return(-1 * loglik)
  }
  logparams <- log(par.start)
  results <- optim(par = logparams, 
                   fn = ggnbd.eLL, 
                   cal.cbs = cal.cbs, 
                   method = "L-BFGS-B",
                   lower = log(min.param.value), 
                   upper = log(max.param.value))
  estimated.params <- exp(results$par)
  return(estimated.params)
}


#' Calculates the log-likelihood of the Gamma/Gompertz/NBD model.
#' 
#' @param params Gamma/Gompertz/NBD parameters - a vector with 'r', 'alpha',
#'   'b', 's' and 'beta' in that order.
#' @param cal.cbs calibration period CBS. It must contain columns for frequency 
#'   'x', for recency 't.x.' and total time observed 'T.cal'. Optionally a 
#'   column 'custs' can be provided, which represents number of customers with a
#'   specific combination of frequency 'x' and 'T.cal'.
#' @return the total log-likelihood for the provided data.
#' @export
ggnbd.cbs.LL <- function(params, cal.cbs) {
  dc.check.model.params(c("r", "alpha", "b", "s", "beta"), params, 
   "ggnbd.cbs.LL")  
  tryCatch(x <- cal.cbs[, "x"], error = function(e) stop("Error in ggnbd.cbs.LL: cal.cbs must have a frequency column labelled \"x\""))
  tryCatch(t.x <- cal.cbs[, "t.x"], error = function(e) stop("Error in ggnbd.cbs.LL: cal.cbs must have a recency column labelled \"t.x\""))
  tryCatch(T.cal <- cal.cbs[, "T.cal"], error = function(e) stop("Error in ggnbd.cbs.LL: cal.cbs must have a column for length of time observed labelled \"T.cal\""))
  if ("custs" %in% colnames(cal.cbs)) {
    custs <- cal.cbs[, "custs"]
  } else {
    custs <- rep(1, length(x))
  }
  return(sum(custs * ggnbd.LL(params, x, t.x, T.cal)))
}


#' Calculates the log-likelihood of the Gamma/Gompertz/NBD model.
#' 
#' @param params Gamma/Gompertz/NBD parameters - a vector with 'r', 'alpha', 
#'   'b', 's' and 'beta' in that order.
#' @param cal.cbs calibration period CBS. It must contain columns for frequency 
#'   'x', for recency 't.x.' and total time observed 'T.cal'. Optionally a 
#'   column 'custs' can be provided, which represents number of customers with a
#'   specific combination of frequency 'x' and 'T.cal'.
#' @return a vector of log-likelihoods
#' @export
ggnbd.LL <- function(params, x, t.x, T.cal) {
  max.length <- max(length(x), length(t.x), length(T.cal))
  if (max.length%%length(x)) 
    warning("Maximum vector length not a multiple of the length of x")
  if (max.length%%length(t.x)) 
    warning("Maximum vector length not a multiple of the length of t.x")
  if (max.length%%length(T.cal)) 
    warning("Maximum vector length not a multiple of the length of T.cal")
  dc.check.model.params(c("r", "alpha", "b", "s", "beta"), params, 
   "ggnbd.LL")
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
  b <- params[3]
  s <- params[4]
  beta <- params[5]
  
  intl <- numeric(0)
  for (i in 1:max.length) {
    intl[i] <- integrate(function(y) (y+alpha)^-(r+x[i]) * (beta+exp(b*y)-1)^-(s+1) * exp(b*y), 
                      lower=t.x[i], upper=T.cal[i])$value
  }
  
  L1 <- lgamma(r+x) - lgamma(r) + r*(log(alpha)-log(alpha+T.cal)) - x*log(alpha+T.cal) + s*(log(beta)-log(beta-1+exp(b*T.cal)))
  L2 <- lgamma(r+x) - lgamma(r) + log(b) + r*log(alpha) + log(s) + s*log(beta) + log(intl)
  
  llh <- log(exp(L1) + exp(L2))
  
  return(llh)
}


#' Uses Gamma/Gompertz/NBD model parameters and a customer's past transaction
#' behavior to return the probability that they are still alive at the end of
#' the calibration period.
#' 
#' @param params Gamma/Gompertz/NBD parameters - a vector with 'r', 'alpha',
#'   'b', 's' and 'beta' in that order.
#' @param x number of repeat transactions in the calibration period T.cal, or a 
#'   vector of calibration period frequencies.
#' @param t.x recency, i.e. length between first and last transaction during 
#'   calibration period.
#' @param T.cal length of calibration period, or a vector of calibration period 
#'   lengths.
#' @return Probability that the customer is still alive at the end of the 
#'   calibration period.
#' @export
ggnbd.PAlive <- function(params, x, t.x, T.cal) {
  max.length <- max(length(x), length(t.x), length(T.cal))
  if (max.length%%length(x)) 
    warning("Maximum vector length not a multiple of the length of x")
  if (max.length%%length(t.x)) 
    warning("Maximum vector length not a multiple of the length of t.x")
  if (max.length%%length(T.cal)) 
    warning("Maximum vector length not a multiple of the length of T.cal")
  dc.check.model.params(c("r", "alpha", "b", "s", "beta"), params, "ggnbd.PAlive")
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
  b <- params[3]
  s <- params[4]
  beta <- params[5]
  P1 <- lgamma(r+x) - lgamma(r)
  P2 <- r * log(alpha/(alpha+T.cal)) + x * log(1/(alpha+T.cal)) + s * log(beta/(beta-1+exp(b*T.cal)))
  P3 <- ggnbd.LL(params, x, t.x, T.cal)
  return(exp(P1 + P2 - P3))
}


#' Uses Gamma/Gompertz/NBD model parameters and a customer's past transaction
#' behavior to return the number of transactions they are expected to make in a
#' given time period.
#' 
#' @param params Gamma/Gompertz/NBD parameters - a vector with 'r', 'alpha',
#'   'b', 's' and 'beta' in that order.
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
#' @export
ggnbd.ConditionalExpectedTransactions <- function(params, T.star, x, t.x, T.cal) {
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
  dc.check.model.params(c("r", "alpha", "b", "s", "beta"), params, 
    "ggnbd.ConditionalExpectedTransactions")
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
  b <- params[3]
  s <- params[4]
  beta <- params[5]
  P1 <- ((r+x)/(alpha+T.cal))
  beta.star <- beta + exp(b*T.cal) - 1
  P2 <- (beta.star/(beta.star+exp(b*T.star)-1))^s * T.star
  P3a <- b * s * beta.star^s
  P3b <- numeric()
  for (i in 1:max.length) {
    P3b[i] <- integrate(function(tau) tau * exp(b*tau) * (beta.star[i]+exp(b*tau)-1)^(-s-1),
                        lower=0, upper=T.star[i])$value
  }
  PAlive <- ggnbd.PAlive(params, x, t.x, T.cal)
  return(P1 * (P2+P3a*P3b) * PAlive)
}


#' Create simulated data set according to Gamma/Gompertz/NBD model assumptions.
#' 
#' @param n number of customers
#' @param T.cal length of calibration period
#' @param T.star length of holdout period
#' @param params Gamma/Gompertz/NBD parameters - a vector with 'r', 'alpha',
#'   'b', 's' and 'beta' in that order.
#' @param return.elog boolean - if TRUE then the event log is returned in 
#'   addition to the CBS summary
#' @return list with elements 'cbs' and 'elog' containing data.frames
#' @export
ggnbd.GenerateData <- function(n, T.cal, T.star, params, return.elog=F) {
  # sample intertransaction timings parameter lambda for each customer
  lambdas <- rgamma(n, shape=params$r, rate=params$alpha)
  
  # sample lifetimes for each customer  
  taus <- (1/params$b)*log(1-params$beta+params$beta/(1-runif(n))^(1/params$s))
  #etas <- rgamma(n, shape=params$s, rate=params$beta)
  #taus <- flexsurv::rgompertz(n, shape=etas, rate=b)
  
  # sample intertransaction timings & churn
  cbs <- data.frame()
  elog <- data.frame(cust=numeric(0), t=numeric(0))
  for (i in 1:n) {
    lambda <- lambdas[i]
    #eta <- etas[i]
    tau <- taus[i]
    # sample 'sufficiently' large amount of inter-transaction times
    minT <- min(T.cal+T.star, tau)
    draws <- max(10, minT * lambda)
    itts <- rexp(draws * 2, rate=lambda)
    if (sum(itts)<minT) itts <- c(itts, rexp(draws * 4, rate=lambda))
    if (sum(itts)<minT) itts <- c(itts, rexp(draws * 8, rate=lambda))
    if (sum(itts)<minT)
      stop("not enough inter-transaction times sampled: ", sum(itts), " < ", tau)
    times <- cumsum(c(0, itts))
    times <- times[times<tau]
    if (return.elog) elog <- rbind(elog, data.frame(cust=i, t=times[times<(T.cal+T.star)]))
    # determine frequency, recency, etc.
    ts.cal <- times[times<T.cal]
    ts.star <- times[times>=T.cal & times<(T.cal+T.star)]
    cbs[i, "x"] <- length(ts.cal)-1
    cbs[i, "t.x"] <- max(ts.cal)
    cbs[i, "T.cal"] <- T.cal
    cbs[i, "alive"] <- tau>T.cal
    cbs[i, "x.star"] <- length(ts.star)
    cbs[i, "lambda"] <- lambda
    #cbs[i, "eta"] <- eta
    cbs[i, "tau"] <- tau
  }
  elog$date <- as.Date("2001/01/01") + elog$t
  out <- list(cbs=cbs)
  if (return.elog) {
    out$elog <- transform(elog, date=as.Date("2001/01/01") + t)
  }
  return(out)
}
