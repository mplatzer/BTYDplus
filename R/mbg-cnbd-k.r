
#' Parameter Estimation for the MBG/CNBD-k model
#' 
#' Estimates parameters for the MBG/CNBD-k model via Maximum Likelihood 
#' Estimation.
#' 
#' @param cal.cbs calibration period CBS. It must contain columns for frequency 
#'   \code{x}, for recency \code{t.x} and total time observed \code{T.cal}. Optionally a 
#'   column \code{custs} can be provided, which represents number of customers with a
#'   specific combination of frequency \code{x}, recency \code{t.x} and \code{T.cal}.
#' @param k specified degree of regularity for Erlang-k distributed 
#'   interpurchase times; needs to be integer-value; if this is not specified, 
#'   then grid search from 1 to 12 is performed; this however requires column
#'   \code{litt} to be present in cal.cbs, which represents sum of logarithmic
#'   interpurchase times during calibration period;
#' @param par.start initial MBG/CNBD-k parameters - a vector with \code{r}, \code{alpha}, 
#'   \code{a} and \code{b} in that order.
#' @param max.param.value the upper bound on parameters
#' @param trace print logging step every \code{trace} iteration
#' @return list of estimated parameters
#' @import BTYD
#' @export
#' @seealso \code{\link{elog2cbs}}
#' @references Platzer, Michael. "Stochastic models of noncontractual consumer 
#'   relationships." Master of Science in Business Administration thesis, Vienna
#'   University of Economics and Business Administration, Austria (2008). 
#'   \url{https://sites.google.com/site/michaelplatzer/stochastic-models-of-noncontractual-consumer-relationships}
#' @example demo/mbg-cnbd-k.r
mbgcnbd.EstimateParameters <- function(cal.cbs, k=NULL, par.start=c(1, 1, 1, 1), max.param.value=10000, trace=0) {
  
  dc.check.model.params(c("r", "alpha", "a", "b"), par.start, 
    "mbgcnbd.EstimateParameters")
  
  # Either 'k' or \code{litt} need to be present
  if (is.null(k) & !"litt" %in% colnames(cal.cbs))
    stop("Either regularity parameter k need to be specified, or a column with logarithmic interpurchase times litt need to be present in cal.cbs")
  
  # if k is not specified we do grid search for k
  if (is.null(k)) {
    params <- list()
    LL <- c()
    for (k in 1:12) {
      params[[k]] <- tryCatch(mbgcnbd.EstimateParameters(cal.cbs, par.start, max.param.value, k=k, trace=trace),
                              error = function(e) { e })
      if (inherits(params[[k]], "error")) {
        params[[k]] <- NULL
        break # stop if parameters could not be estimated, e.g. if mbgcnbd.LL returns Inf
      }
      LL[k] <- mbgcnbd.cbs.LL(params[[k]], cal.cbs)
      if (k > 4 && LL[k] < LL[k-1] && LL[k-1] < LL[k-2]) 
        break # stop if LL gets worse for increasing k
    }
    k <- which.max(LL)
    return(params[[k]])
  }
  
  # if \code{litt} is missing, we set it to zero, so that mbgcnbd.cbs.LL does not
  # complain; however this makes LL values for different k values not comparable
  if (!"litt" %in% colnames(cal.cbs))
    cal.cbs[, "litt"] <- 0 
  
  count <- 0  
  mbgcnbd.eLL <- function(params, k, cal.cbs, max.param.value) {
    params <- exp(params)
    params[params > max.param.value] <- max.param.value
    params <- c(k, params)
    loglik <- mbgcnbd.cbs.LL(params, cal.cbs)
    count <<- count + 1
    if (trace>0 & count%%trace==0)
      cat("mbgcnbd.EstimateParameters - iter", count, ":", sprintf("%12.2f", loglik), ":", sprintf("%10.6f", params), "\n")
    return(-1 * loglik)
  }
  
  logparams <- log(par.start)
  results <- optim(logparams, mbgcnbd.eLL, cal.cbs = cal.cbs, k = k,
    max.param.value = max.param.value, method = "L-BFGS-B")
  estimated.params <- exp(results$par)
  estimated.params[estimated.params > max.param.value] <- max.param.value
  estimated.params <- c(k, estimated.params)
  names(estimated.params) <- c("k", "r", "alpha", "a", "b")
  return(estimated.params)
}


#' Calculate the log-likelihood of the MBG/CNBD-k model
#' 
#' @param params MBG/CNBD-k parameters - a vector with 'k', \code{r}, \code{alpha}, \code{a} 
#'   and \code{b} in that order.
#' @param cal.cbs calibration period CBS. It must contain columns for frequency 
#'   \code{x}, for recency \code{t.x}, for sum of logarithmic interpurchase times \code{litt}
#'   and total time observed \code{T.cal}. Optionally a column \code{custs} can be
#'   provided, which represents number of customers with a specific combination
#'   of frequency \code{x} and \code{T.cal}.
#' @return the total log-likelihood for the provided data.
#' @export
#' @seealso \code{\link{mbgcnbd.EstimateParameters}} 
mbgcnbd.cbs.LL <- function(params, cal.cbs) {
  dc.check.model.params(c("k", "r", "alpha", "a", "b"), params, 
   "mbgcnbd.cbs.LL")
  tryCatch(x <- cal.cbs$x, error = function(e) stop("Error in mbgcnbd.cbs.LL: cal.cbs must have a frequency column labelled \"x\""))
  tryCatch(t.x <- cal.cbs$t.x, error = function(e) stop("Error in mbgcnbd.cbs.LL: cal.cbs must have a recency column labelled \"t.x\""))
  tryCatch(T.cal <- cal.cbs$T.cal, error = function(e) stop("Error in mbgcnbd.cbs.LL: cal.cbs must have a column for length of time observed labelled \"T.cal\""))
  tryCatch(litt <- cal.cbs$litt, error = function(e) stop("Error in mbgcnbd.cbs.LL: cal.cbs must have a column for sum over logarithmic inter-transaction-times labelled \"litt\""))
  if ("custs" %in% colnames(cal.cbs)) {
    custs <- cal.cbs$custs
  } else {
    custs <- rep(1, length(x))
  }
  return(sum(custs * mbgcnbd.LL(params, x, t.x, T.cal, litt)))
}


#' Calculate the log-likelihood of the MBG/CNBD-k model
#' 
#' @param params MBG/CNBD-k parameters - a vector with 'k', \code{r}, \code{alpha}, \code{a}
#'   and \code{b} in that order.
#' @param x frequency, i.e. number of re-purchases
#' @param t.x recency, i.e. time elapsed from first purchase to last purchase
#' @param T.cal total time of observation period
#' @param litt sum of logarithmic interpurchase times
#' @return a vector of log-likelihoods
#' @export
#' @seealso \code{\link{mbgcnbd.EstimateParameters}}
mbgcnbd.LL <- function(params, x, t.x, T.cal, litt) {
  max.length <- max(length(x), length(t.x), length(T.cal))
  if (max.length%%length(x)) 
    warning("Maximum vector length not a multiple of the length of x")
  if (max.length%%length(t.x)) 
    warning("Maximum vector length not a multiple of the length of t.x")
  if (max.length%%length(T.cal)) 
    warning("Maximum vector length not a multiple of the length of T.cal")
  if (max.length%%length(litt)) 
    warning("Maximum vector length not a multiple of the length of litt")
  dc.check.model.params(c("k", "r", "alpha", "a", "b"), params, 
   "mbgcnbd.LL")
  if (params[1] != floor(params[1]) | params[1] < 1)
    stop("k must be integer being greater or equal to 1.")
  if (any(x < 0) || !is.numeric(x)) 
    stop("x must be numeric and may not contain negative numbers.")
  if (any(t.x < 0) || !is.numeric(t.x)) 
    stop("t.x must be numeric and may not contain negative numbers.")
  if (any(T.cal < 0) || !is.numeric(T.cal)) 
    stop("T.cal must be numeric and may not contain negative numbers.")
  x <- rep(x, length.out = max.length)
  t.x <- rep(t.x, length.out = max.length)
  T.cal <- rep(T.cal, length.out = max.length)
  litt <- rep(litt, length.out = max.length)
  k <- params[1]
  r <- params[2]
  alpha <- params[3]
  a <- params[4]
  b <- params[5]
  P0 <- (k-1) * litt - x*log(factorial(k-1))
  P1 <- lgamma(a+b) + lgamma(b+x+1) - lgamma(b) - lgamma(a+b+x+1)
  P2 <- lgamma(r+k*x) + r * log(alpha) - lgamma(r)
  P3 <- -1 * (r+k*x) * log(alpha+T.cal)
  P4 <- (a/(b+x))*((alpha+T.cal)/(alpha+t.x))^(r+k*x)
  P5 <- 1
  if (k > 1) {
    for (j in 1:(k-1)) {
      P5a <- 1
      for (i in 0:(j-1)) P5a <- P5a * (r + k * x + i)
      P5 <- P5 + (P5a * (T.cal - t.x)^j) / (factorial(j) * (alpha + T.cal)^j)
    }
  }
  return(P0 + P1 + P2 + P3 + log(P4 + P5))
}


#' MBG/CNBD-k P(alive)
#' 
#' Uses MBG/CNBD-k model parameters and a customer's past transaction behavior 
#' to return the probability that they are still alive at the end of the 
#' calibration period.
#' 
#' @param params MBG/CNBD-k parameters - a vector with 'k', \code{r}, \code{alpha}, \code{a}
#'   and \code{b} in that order.
#' @param x number of repeat transactions in the calibration period T.cal, or a 
#'   vector of calibration period frequencies.
#' @param t.x recency, i.e. length between first and last transaction during 
#'   calibration period.
#' @param T.cal length of calibration period, or a vector of calibration period 
#'   lengths.
#' @return Probability that the customer is still alive at the end of the 
#'   calibration period.
#' @export
#' @example demo/mbg-cnbd-k.r
#' @seealso \code{\link{mbgcnbd.EstimateParameters}}
mbgcnbd.PAlive <- function(params, x, t.x, T.cal) {
  max.length <- max(length(x), length(t.x), length(T.cal))
  if (max.length%%length(x)) 
    warning("Maximum vector length not a multiple of the length of x")
  if (max.length%%length(t.x)) 
    warning("Maximum vector length not a multiple of the length of t.x")
  if (max.length%%length(T.cal)) 
    warning("Maximum vector length not a multiple of the length of T.cal")
  dc.check.model.params(c("k", "r", "alpha", "a", "b"), params, 
   "mbgcnbd.PAlive")  
  if (params[1] != floor(params[1]) | params[1] < 1)
    stop("k must be integer being greater or equal to 1.")
  if (any(x < 0) || !is.numeric(x)) 
    stop("x must be numeric and may not contain negative numbers.")
  if (any(t.x < 0) || !is.numeric(t.x)) 
    stop("t.x must be numeric and may not contain negative numbers.")
  if (any(T.cal < 0) || !is.numeric(T.cal)) 
    stop("T.cal must be numeric and may not contain negative numbers.")
  x <- rep(x, length.out = max.length)
  t.x <- rep(t.x, length.out = max.length)
  T.cal <- rep(T.cal, length.out = max.length)
  k <- params[1]
  r <- params[2]
  alpha <- params[3]
  a <- params[4]
  b <- params[5]
  P1 <- (a / (b + x)) * ((alpha + T.cal) / (alpha + t.x))^(r + k * x)
  P2 <- 1
  if (k>1) {
    for (j in 1:(k-1)) {
      P2a <- 1
      for (i in 0:(j-1)) P2a <- P2a * (r+k*x+i)
      P2 <- P2 + ((P2a * (T.cal - t.x)^j) / (factorial(j) * (alpha + T.cal)^j))
    }
  }
  return (1 / (1 + P1 / P2))
}


#' MBG/CNBD-k Conditional Expected Transactions
#' 
#' Uses MBG/CNBD-k model parameters and a customer's past transaction behavior 
#' to return the number of transactions they are expected to make in a given 
#' time period.
#' 
#' @param params MBG/CNBD-k parameters - a vector with 'k', \code{r}, \code{alpha}, \code{a}
#'   and \code{b} in that order.
#' @param T.star length of time for which we are calculating the expected number
#'   of transactions.
#' @param x number of repeat transactions in the calibration period T.cal, or a 
#'   vector of calibration period frequencies.
#' @param t.x recency, i.e. length between first and last transaction during 
#'   calibration period.
#' @param T.cal length of calibration period, or a vector of calibration period 
#'   lengths.
#' @param method `exact` or `approx`; method `exact` computes `P(Y(s, s+t)=j)`
#'   for sufficiently large j and sums over these; method `approx` is only
#'   approximative and uses closed-form results from MBG/NBD and down-scales
#'   them by factor k; the latter method is significantly faster;
#' @return Number of transactions a customer is expected to make in a time 
#'   period of length t, conditional on their past behavior. If any of the input
#'   parameters has a length greater than 1, this will be a vector of expected 
#'   number of transactions.
#' @import gsl
#' @export
#' @example demo/mbg-cnbd-k.r
#' @seealso \code{\link{mbgcnbd.EstimateParameters}}
mbgcnbd.ConditionalExpectedTransactions <- function(params, T.star, x, t.x, T.cal, method="exact") {
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
  dc.check.model.params(c("k", "r", "alpha", "a", "b"), params,
    "mbgcnbd.ConditionalExpectedTransactions")
  if (params[1] != floor(params[1]) | params[1] < 1)
    stop("k must be integer being greater or equal to 1.")
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
  k <- params[1]
  r <- params[2]
  alpha <- params[3]
  a <- params[4]
  b <- params[5]
  # calculate probabilities for encountering 0 to k-1 censored events within T.cal-t.x
  probs <- matrix(NA_real_, ncol=k, nrow=length(x))
  for (i in 0:(k-1)) {
    probs[, i+1] <- ((T.cal-t.x)^i / factorial(i)) * (exp(lgamma(r+x+i) - lgamma(r+x))) * (alpha+1)^(r+x) / (alpha+1+T.cal-t.x)^(r+x+i)
  }
  probs <- probs / rowSums(probs)
  # implementation of two methods for calculating conditional expected number of transactions
  if (method=='approx') {
    if (k>1) warning('method `approx` will result in approximative results for k>1')
    # iterate over cases of encountering 0 to k-1 censored events within T.cal-t.x
    e_x <- numeric(length(x))
    for (i in 0:(k-1)) {
      a_ <- a
      b_ <- b+k*x+i
      r_ <- r+k*x+i
      alpha_ <- alpha+T.cal
      e_x <- e_x + probs[,i+1] * (1/k) * ((a_+b_)/(a_-1)) * 
        (1-(alpha_/(alpha_+T.star))^r_ * gsl::hyperg_2F1(r_, b_+1, a_+b_, T.star/(alpha_+T.star)))
    }
  } else if (method=='exact1') {
    warning('this is not working :(')
    nbd_pdf <- function(t, j, r, alpha) exp(j * log(t) - lfactorial(j) + r * log(alpha) + lgamma(r+j) - lgamma(r) - (r+j) * log(alpha+t))
    nbd_cdf <- function(t, j, r, alpha) apply(sapply(0:j, function(j) nbd_pdf(t, j, r, alpha)), 1, sum)
    max_x  <- min(500, round(max(cbs$T.star * cbs$x / cbs$T.cal))*2) #set max number of x for which we compute P(X(t)=x)
    seq_x  <- 1:max_x
    a_     <- a
    b_     <- b+x
    r_     <- r+x
    alpha_ <- alpha+T.cal
    e_x <- numeric(length(x))
    for (j in seq_x) e_x <- e_x + j * (
      nbd_pdf(cbs$T.star, j, r_, alpha_) * beta(a_, b_+j+1) / beta(a_, b_) + 
      (1-nbd_cdf(cbs$T.star, j-1, r_, alpha_)) * beta(a_+1, b_+j) / beta(a_, b_)
    )
  } else if (method=='exact-base') {
    nbd_pdf <- function(t, j, r, alpha) exp(j * log(t) - lfactorial(j) + r * log(alpha) + lgamma(r+j) - lgamma(r) - (r+j) * log(alpha+t))
    nbd_cdf <- function(t, j, r, alpha) sum(sapply(0:j, function(j) nbd_pdf(t, j, r, alpha)))
    max_x  <- min(500, round(max(cbs$T.star * cbs$x / cbs$T.cal))*2) #set max number of x for which we compute P(X(t)=x)
    seq_x  <- 1:max_x
    e_x <- 0
    for (j in seq_x) e_x <- e_x + j * (
      nbd_pdf(T.star, j, r, alpha) * beta(a, b+j+1) / beta(a, b) + 
      (1-nbd_cdf(T.star, j-1, r, alpha)) * beta(a+1, b+j) / beta(a, b)
    )
  } else if (method=='approx-base') {
    (b/(a-1)) * (1-(alpha/(alpha+T.star))^r * gsl::hyperg_2F1(r, b+1, a+b, T.star/(alpha+T.star)))
  } else if (method=='exact') {
    warning('this is not working :(')
    # sum over sufficiently large series of P(X(t)=x) terms
    nbd_pdf <- function(t, j, r, alpha) exp(j * log(t) - lfactorial(j) + r * log(alpha) + lgamma(r+j) - lgamma(r) - (r+j) * log(alpha+t))
    e_x <- numeric(length(x))
    for (i in 0:(k-1)) {
      max_x  <- min(500, round(max(T.star * cbs$x / cbs$T.cal))*2) #set max number of x for which we compute P(X(t)=x)
      seq_x  <- matrix(1:max_x,     nrow = nrow(cbs), ncol = max_x, byrow = T)
      a_     <- matrix(a,           nrow = nrow(cbs), ncol = max_x)
      b_     <- matrix(b+k*x+i,     nrow = nrow(cbs), ncol = max_x)
      r_     <- matrix(r+k*x+i,     nrow = nrow(cbs), ncol = max_x)
      alpha_ <- matrix(alpha+T.cal, nrow = nrow(cbs), ncol = max_x)
      seq_j  <- matrix(0:((max_x+1)*k-1), nrow=nrow(cbs), ncol=(max_x+1)*k, byrow=T)
      mat_pdf <- nbd_pdf(T.star, seq_j, r_[,1], alpha_[,1])
      mat_cdf <- t(apply(mat_pdf, 1, cumsum))
      F1 <- exp(lbeta(a_, b_+seq_x+1) - lbeta(a_, b_))
      T1 <- mat_pdf %*% matrix(c(rep(c(rep(1, k), rep(0, (max_x+1)*k)), max_x), rep(1, k)), nrow=(max_x+1)*k, ncol=max_x+1)[,-1]
      F2 <- exp(lbeta(a_+1, b_+seq_x) - lbeta(a_, b_))
      T2 <- (1-mat_cdf[, k*(1:max_x)])
      e_x <- e_x + probs[,i+1] * rowSums(seq_x * (F1*T1 + F2*T2))
    }
  }
  p_alive <- mbgcnbd.PAlive(params, x, t.x, T.cal)
  return(e_x * p_alive)    
}



#' MBG/CNBD-k Unconditional Probability Distribution of Transactions
#' 
#' Uses MBG/CNBD-k model parameters to return the probability distribution of
#' purchase frequencies for a random customer in a given time period, i.e.
#' P(X(t)=x|r,alpha,a,b)
#' 
#' @param params MBG/CNBD-k parameters - a vector with 'k', \code{r}, \code{alpha}, \code{a} and \code{b}
#'   in that order.
#' @param t length of time for which we are calculating the expected number of 
#'   transactions.
#' @param x number of transactions for which probability is calculated.
#' @return P(X(t)=x|r,alpha,a,b). If any of the input parameters has a length
#'   greater than 1, this will be a vector of expected number of transactions.
#' @export
#' @seealso \code{\link{mbgcnbd.EstimateParameters}}
mbgcnbd.Px <- function(params, t, x) {
  dc.check.model.params(c("k", "r", "alpha", "a", "b"), params, 
                        "mbgcnbd.Px")
  if (params[1] != floor(params[1]) | params[1] < 1)
    stop("k must be integer being greater or equal to 1.")  
  
  k <- params[1]
  r <- params[2]
  alpha <- params[3]
  a <- params[4]
  b <- params[5]
  
  nbd.Px <- function(params, t, x) {
    return ((gamma(r+x)*alpha^r*t^x)/(factorial(x)*gamma(r)*(alpha+t)^(r+x)))
  }

  P1 <- (gamma(b+x+1) * gamma(a+b)) / (gamma(b) * gamma(a+b+x+1))
  P2a <- sum(nbd.Px(params, t, (k*x):(k*x+k-1)))
  P2b <- a/(b+x)
  if (x>0) P2b <- P2b * (1-sum(nbd.Px(params, t, 0:(k*x-1))))
  return (P1 * (P2a + P2b))
}


#' Simulate data according to MBG/CNBD-k model assumptions
#' 
#' @param n number of customers
#' @param T.cal length of calibration period
#' @param T.star length of holdout period
#' @param params MBG/CNBD-k parameters - a vector with 'k', \code{r}, \code{alpha}, \code{a} and \code{b}
#'   in that order.
#' @param return.elog boolean - if \code{TRUE} then the event log is returned in 
#'   addition to the CBS summary
#' @return list with elements \code{cbs} and \code{elog} containing data.frames
#' @export
mbgcnbd.GenerateData <- function(n, T.cal, T.star=T.cal, params, return.elog=FALSE) {
  # check model parameters
  dc.check.model.params(c("k", "r", "alpha", "a", "b"), params,
                        "mbgcnbd.GenerateData")
  if (params[1] != floor(params[1]) | params[1] < 1)
    stop("k must be integer being greater or equal to 1.")  
  
  k <- params[1]
  r <- params[2]
  alpha <- params[3]
  a <- params[4]
  b <- params[5]  

  if (length(T.cal)==1) T.cal <- rep(T.cal, n)  
  if (length(T.star)==1) T.star <- rep(T.star, n)  
  
  # sample intertransaction timings parameter lambda for each customer
  lambdas <- rgamma(n, shape=r, rate=alpha)
  
  # sample churn-probability p for each customer
  ps <- rbeta(n, a, b)
  
  # sample intertransaction timings & churn
  cbs_list <- list()
  elog_list <- list()
  for (i in 1:n) {
    p <- ps[i]
    lambda <- lambdas[i]
    # sample no. of transactions until churn
    churn <- which.max(rbinom(min(10000, 10/p), 1, p)) - 1
    # sample transaction times
    times <- cumsum(c(0, rgamma(churn, shape=k, rate=lambda)))
    if (return.elog)
      elog_list[[i]] <- data.frame(cust=i, t=times[times<(T.cal[i]+T.star[i])])
    # determine frequency, recency, etc.
    ts.cal  <- times[times<T.cal[i]]
    ts.star <- times[times>=T.cal[i] & times<(T.cal[i]+T.star[i])]
    cbs_list[[i]] <- list(cust   = i,
                          x      = length(ts.cal)-1,
                          t.x    = max(ts.cal),
                          litt   = sum(log(diff(ts.cal))),
                          churn  = churn,
                          alive  = churn > (length(ts.cal)-1),
                          x.star = length(ts.star))
  }
  cbs <- do.call(rbind.data.frame, cbs_list)
  cbs$lambda <- lambdas
  cbs$p      <- ps
  cbs$k      <- k
  cbs$T.cal  <- T.cal
  cbs$T.star <- T.star
  rownames(cbs) <- NULL
  out <- list(cbs=cbs)
  if (return.elog) {
    elog <- do.call(rbind.data.frame, elog_list)
    out$elog <- elog
  }
  return(out)
}
