
#' Parameter Estimation for the BG/CNBD-k model
#' 
#' Estimates parameters for the BG/CNBD-k via Maximum Likelihood Estimation.
#' 
#' @param cal.cbs calibration period CBS. It must contain columns for frequency 
#'   \code{x}, for recency \code{t.x} and total time observed \code{T.cal}.
#'   Optionally a column \code{custs} can be provided, which represents number
#'   of customers with a specific combination of frequency \code{x}, recency
#'   \code{t.x} and \code{T.cal}.
#' @param k specified degree of regularity for Erlang-k distributed 
#'   interpurchase times; needs to be integer-value; if this is not specified, 
#'   then grid search from 1 to 12 is performed; this however requires column 
#'   \code{litt} to be present in cal.cbs, which represents sum of logarithmic 
#'   interpurchase times during calibration period;
#' @param par.start initial BG/CNBD-k parameters - a vector with \code{r},
#'   \code{alpha}, \code{a} and \code{b} in that order.
#' @param max.param.value the upper bound on parameters
#' @param trace print logging step every \code{trace} iteration
#' @param dropout_at_zero Boolean; the mbg-methods are simple wrapper methods,
#'   which set this parameter to TRUE
#' @return list of estimated parameters
#' @import BTYD
#' @export
#' @seealso \code{\link{elog2cbs}}
#' @references Platzer Michael, and Thomas Reutterer (forthcoming)
#' @example demo/bg-cnbd-k.r
bgcnbd.EstimateParameters <- function(cal.cbs, 
                                      k = NULL, 
                                      par.start = c(1, 3, 1, 3), 
                                      max.param.value = 10000, 
                                      trace = 0,
                                      dropout_at_zero = FALSE) {
  
  dc.check.model.params(c("r", "alpha", "a", "b"), par.start, 
    "bgcnbd.EstimateParameters")
  
  # Either \code{k} or \code{litt} need to be present
  if (is.null(k) & !"litt" %in% colnames(cal.cbs))
    stop("Either regularity parameter k need to be specified, or a column with logarithmic interpurchase times litt need to be present in cal.cbs")
  
  # if k is not specified we do grid search for k
  if (is.null(k)) {
    params <- list()
    LL <- c()
    for (k in 1:12) {
      params[[k]] <- tryCatch(bgcnbd.EstimateParameters(cal.cbs, par.start, max.param.value, k=k, trace=trace),
                              error = function(e) { e })
      if (inherits(params[[k]], "error")) {
        params[[k]] <- NULL
        break # stop if parameters could not be estimated, e.g. if bgcnbd.LL returns Inf
      }
      LL[k] <- bgcnbd.cbs.LL(params[[k]], cal.cbs)
      if (k > 4 && LL[k] < LL[k-1] && LL[k-1] < LL[k-2]) 
        break # stop if LL gets worse for increasing k
    }
    k <- which.max(LL)
    return(params[[k]])
  }
  
  # if \code{litt} is missing, we set it to zero, so that bgcnbd.cbs.LL does not
  # complain; however this makes LL values for different k values not comparable
  if (!"litt" %in% colnames(cal.cbs))
    cal.cbs[, "litt"] <- 0 
  
  count <- 0
  bgcnbd.eLL <- function(params, k, cal.cbs, max.param.value) {
    params <- exp(params)
    params[params > max.param.value] <- max.param.value
    params <- c(k, params)
    loglik <- bgcnbd.cbs.LL(params, cal.cbs, dropout_at_zero)
    count <<- count + 1
    if (trace>0 & count%%trace==0)
      cat("bgcnbd.EstimateParameters - iter", count, ":", sprintf("%12.2f", loglik), ":", sprintf("%10.6f", params), "\n")
    return(-1 * loglik)
  }
  
  logparams <- log(par.start)
  results <- optim(logparams, bgcnbd.eLL, cal.cbs = cal.cbs, k = k,
    max.param.value = max.param.value, method = "L-BFGS-B")
  estimated.params <- exp(results$par)
  estimated.params[estimated.params > max.param.value] <- max.param.value
  estimated.params <- c(k, estimated.params)
  names(estimated.params) <- c("k", "r", "alpha", "a", "b")
  return(estimated.params)
}


#' Calculate the log-likelihood of the BG/CNBD-k model
#' 
#' @param params BG/CNBD-k parameters - a vector with \code{k}, \code{r},
#'   \code{alpha}, \code{a} and \code{b} in that order.
#' @param cal.cbs calibration period CBS. It must contain columns for frequency 
#'   \code{x}, for recency \code{t.x}, for sum of logarithmic interpurchase
#'   times \code{litt} and total time observed \code{T.cal}. Optionally a column
#'   \code{custs} can be provided, which represents number of customers with a
#'   specific combination of frequency \code{x} and \code{T.cal}.
#' @param dropout_at_zero Boolean; the mbg-methods are simple wrapper methods,
#'   which set this parameter to TRUE
#' @return the total log-likelihood for the provided data.
#' @export
#' @seealso \code{\link{bgcnbd.EstimateParameters}}
bgcnbd.cbs.LL <- function(params, cal.cbs, dropout_at_zero = FALSE) {
  dc.check.model.params(c("k", "r", "alpha", "a", "b"), params, 
   "bgcnbd.cbs.LL")
  tryCatch(x <- cal.cbs$x, error = function(e) stop("Error in bgcnbd.cbs.LL: cal.cbs must have a frequency column labelled \"x\""))
  tryCatch(t.x <- cal.cbs$t.x, error = function(e) stop("Error in bgcnbd.cbs.LL: cal.cbs must have a recency column labelled \"t.x\""))
  tryCatch(T.cal <- cal.cbs$T.cal, error = function(e) stop("Error in bgcnbd.cbs.LL: cal.cbs must have a column for length of time observed labelled \"T.cal\""))
  tryCatch(litt <- cal.cbs$litt, error = function(e) stop("Error in bgcnbd.cbs.LL: cal.cbs must have a column for sum over logarithmic inter-transaction-times labelled \"litt\""))
  if ("custs" %in% colnames(cal.cbs)) {
    custs <- cal.cbs$custs
  } else {
    custs <- rep(1, length(x))
  }
  return(sum(custs * bgcnbd.LL(params, x, t.x, T.cal, litt, dropout_at_zero)))
}


#' Calculate the log-likelihood of the BG/CNBD-k model
#' 
#' @param params BG/CNBD-k parameters - a vector with \code{k}, \code{r}, \code{alpha}, \code{a}
#'   and \code{b} in that order.
#' @param x frequency, i.e. number of re-purchases
#' @param t.x recency, i.e. time elapsed from first purchase to last purchase
#' @param T.cal total time of observation period
#' @param litt sum of logarithmic interpurchase times
#' @param dropout_at_zero Boolean; the mbg-methods are simple wrapper methods,
#'   which set this parameter to TRUE
#' @return a vector of log-likelihoods
#' @export
#' @seealso \code{\link{bgcnbd.EstimateParameters}}
bgcnbd.LL <- function(params, x, t.x, T.cal, litt, dropout_at_zero = FALSE) {
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
   "bgcnbd.LL")
  if (params[1] != floor(params[1]) | params[1] < 1)
    stop("k must be integer being greater or equal to 1.")
  if (any(x < 0) || !is.numeric(x)) 
    stop("x must be numeric and may not contain negative numbers.")
  if (any(t.x < 0) || !is.numeric(t.x)) 
    stop("t.x must be numeric and may not contain negative numbers.")
  if (any(T.cal < 0) || !is.numeric(T.cal)) 
    stop("T.cal must be numeric and may not contain negative numbers.")
  x     <- rep(x,     length.out = max.length)
  t.x   <- rep(t.x,   length.out = max.length)
  T.cal <- rep(T.cal, length.out = max.length)
  litt  <- rep(litt,  length.out = max.length)
  k     <- params[1]
  r     <- params[2]
  alpha <- params[3]
  a     <- params[4]
  b     <- params[5]
  P1 <- (k-1) * litt - x*log(factorial(k-1))
  P2 <- lbeta(a, b + x + ifelse(dropout_at_zero, 1, 0)) - lbeta(a, b)
  P3 <- lgamma(r+k*x) - lgamma(r) + r*log(alpha)
  P4 <- -1 * (r+k*x) * log(alpha+T.cal)
  S1 <- as.numeric(x>0) * a / (b + x - 1 + ifelse(dropout_at_zero, 1, 0)) * ((alpha+T.cal)/(alpha+t.x))^(r+k*x)
  S2 <- 1
  if (k > 1) {
    for (j in 1:(k-1)) {
      S2a <- 1
      for (i in 0:(j-1)) S2a <- S2a * (r + k * x + i)
      S2 <- S2 + (S2a * (T.cal - t.x)^j) / (factorial(j) * (alpha + T.cal)^j)
    }
  }
  return(P1 + P2 + P3 + P4 + log(S1 + S2))
}


#' BG/CNBD-k P(alive)
#' 
#' Uses BG/CNBD-k model parameters and a customer's past transaction behavior 
#' to return the probability that they are still alive at the end of the 
#' calibration period.
#' 
#' @param params BG/CNBD-k parameters - a vector with \code{k}, \code{r}, \code{alpha}, \code{a}
#'   and \code{b} in that order.
#' @param x number of repeat transactions in the calibration period T.cal, or a 
#'   vector of calibration period frequencies.
#' @param t.x recency, i.e. length between first and last transaction during 
#'   calibration period.
#' @param T.cal length of calibration period, or a vector of calibration period 
#'   lengths.
#' @param dropout_at_zero Boolean; the mbg-methods are simple wrapper methods,
#'   which set this parameter to TRUE
#' @return Probability that the customer is still alive at the end of the 
#'   calibration period.
#' @export
#' @seealso \code{\link{bgcnbd.EstimateParameters}}
bgcnbd.PAlive <- function(params, x, t.x, T.cal, dropout_at_zero = FALSE) {
  max.length <- max(length(x), length(t.x), length(T.cal))
  if (max.length%%length(x)) 
    warning("Maximum vector length not a multiple of the length of x")
  if (max.length%%length(t.x)) 
    warning("Maximum vector length not a multiple of the length of t.x")
  if (max.length%%length(T.cal)) 
    warning("Maximum vector length not a multiple of the length of T.cal")
  dc.check.model.params(c("k", "r", "alpha", "a", "b"), params, 
   "bgcnbd.PAlive")  
  if (params[1] != floor(params[1]) | params[1] < 1)
    stop("k must be integer being greater or equal to 1.")
  if (any(x < 0) || !is.numeric(x)) 
    stop("x must be numeric and may not contain negative numbers.")
  if (any(t.x < 0) || !is.numeric(t.x)) 
    stop("t.x must be numeric and may not contain negative numbers.")
  if (any(T.cal < 0) || !is.numeric(T.cal)) 
    stop("T.cal must be numeric and may not contain negative numbers.")
  x     <- rep(x,     length.out = max.length)
  t.x   <- rep(t.x,   length.out = max.length)
  T.cal <- rep(T.cal, length.out = max.length)
  k     <- params[1]
  r     <- params[2]
  alpha <- params[3]
  a     <- params[4]
  b     <- params[5]
  P1 <- (a / (b + x - 1 + ifelse(dropout_at_zero, 1, 0))) * ((alpha + T.cal) / (alpha + t.x))^(r + k * x)
  P2 <- 1
  if (k>1) {
    for (j in 1:(k-1)) {
      P2a <- 1
      for (i in 0:(j-1)) P2a <- P2a * (r+k*x+i)
      P2 <- P2 + ((P2a * (T.cal - t.x)^j) / (factorial(j) * (alpha + T.cal)^j))
    }
  }
  palive <- (1 / (1 + P1 / P2))
  if (dropout_at_zero==FALSE) palive[x==0] <- 1
  return (palive)
}


#' BG/CNBD-k Conditional Expected Transactions
#' 
#' Uses BG/CNBD-k model parameters and a customer's past transaction behavior 
#' to return the number of transactions they are expected to make in a given 
#' time period.
#' 
#' @param params BG/CNBD-k parameters - a vector with \code{k}, \code{r}, \code{alpha}, \code{a}
#'   and \code{b} in that order.
#' @param T.star length of time for which we are calculating the expected number
#'   of transactions.
#' @param x number of repeat transactions in the calibration period T.cal, or a 
#'   vector of calibration period frequencies.
#' @param t.x recency, i.e. length between first and last transaction during 
#'   calibration period.
#' @param T.cal length of calibration period, or a vector of calibration period 
#'   lengths.
#' @param dropout_at_zero Boolean; the mbg-methods are simple wrapper methods,
#'   which set this parameter to TRUE
#' @return Number of transactions a customer is expected to make in a time 
#'   period of length t, conditional on their past behavior. If any of the input
#'   parameters has a length greater than 1, this will be a vector of expected 
#'   number of transactions.
#' @import gsl
#' @export
#' @seealso \code{\link{bgcnbd.EstimateParameters}}
bgcnbd.ConditionalExpectedTransactions <- function(params, T.star, x, t.x, T.cal, dropout_at_zero = FALSE) {
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
    "bgcnbd.ConditionalExpectedTransactions")
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
  x      <- rep(x,     length.out = max.length)
  t.x    <- rep(t.x,   length.out = max.length)
  T.cal  <- rep(T.cal, length.out = max.length)
  T.star <- rep(T.star, length.out = max.length)
  k      <- params[1]
  r      <- params[2]
  alpha  <- params[3]
  a      <- params[4]
  b      <- params[5]
  if (k>1) cat("note: conditional expected transactions can only be approximated for k>1\n")
  G <- function(r, alpha, a, b) 1 - (alpha/(alpha+T.star))^r * gsl::hyperg_2F1(r, b+1, a+b, T.star/(alpha+T.star))
  P1 <- (a+b+x-1+ifelse(dropout_at_zero, 1, 0)) / (a-1)
  P2 <- G(r+x, k*alpha+T.cal, a, b+x-1+ifelse(dropout_at_zero, 1, 0))
  P3 <- bgcnbd.PAlive(params, x, t.x, T.cal, dropout_at_zero)
  return (P1 * P2 * P3)
}


#' BG/CNBD-k Probability Mass Function
#' 
#' Uses BG/CNBD-k model parameters to return the probability distribution of
#' purchase frequencies for a random customer in a given time period, i.e.
#' P(X(t)=x|r,alpha,a,b)
#' 
#' @param params BG/CNBD-k parameters - a vector with \code{k}, \code{r}, \code{alpha}, \code{a} and \code{b}
#'   in that order.
#' @param t length of time for which we are calculating the expected number of 
#'   transactions.
#' @param x number of transactions for which probability is calculated.
#' @param dropout_at_zero Boolean; the mbg-methods are simple wrapper methods,
#'   which set this parameter to TRUE
#' @return P(X(t)=x|r,alpha,a,b). If any of the input parameters has a length
#'   greater than 1, this will be a vector of expected number of transactions.
#' @export
#' @seealso \code{\link{bgcnbd.EstimateParameters}}
bgcnbd.pmf <- function(params, t, x, dropout_at_zero = FALSE) {
  if (length(t)>1 & length(x)>1 & length(t)!=length(x)) {
    stop("parameters t and x must be of same length")
  } else if (length(t)>1 | length(x)>1) {
    max.length <- max(length(t), length(x))
    t <- rep(t, length.out = max.length)
    x <- rep(x, length.out = max.length)
    return (sapply(1:max.length, function(i) bgcnbd.pmf(params, t[i], x[i], dropout_at_zero)))
  }
  dc.check.model.params(c("k", "r", "alpha", "a", "b"), params, 
                        "bgcnbd.pmf")
  if (params[1] != floor(params[1]) | params[1] < 1)
    stop("k must be integer being greater or equal to 1.")
  k      <- params[1]
  r      <- params[2]
  alpha  <- params[3]
  a      <- params[4]
  b      <- params[5]
  
  nbd.pmf <- function(params, t, x) {
    return ((gamma(r+x)*alpha^r*t^x)/(factorial(x)*gamma(r)*(alpha+t)^(r+x)))
  }

  P1 <- beta(a, b+x+ifelse(dropout_at_zero, 1, 0)) / beta(a, b)
  P2a <- sum(nbd.pmf(params, t, (k*x):(k*x+k-1)))
  if (dropout_at_zero==FALSE & x==0) {
    P2b <- 0
  } else {
    P2b <- a/(b+x-1+ifelse(dropout_at_zero, 1, 0))
    if (x>0) P2b <- P2b * (1-sum(nbd.pmf(params, t, 0:(k*x-1))))
  }
  return (P1 * (P2a + P2b))
}


#' Simulate data according to BG/CNBD-k model assumptions
#' 
#' @param n number of customers
#' @param T.cal length of calibration period
#' @param T.star length of holdout period
#' @param params BG/CNBD-k parameters - a vector with \code{k}, \code{r}, \code{alpha}, \code{a} and \code{b}
#'   in that order.
#' @param dropout_at_zero Boolean; the mbg-methods are simple wrapper methods,
#'   which set this parameter to TRUE
#' @param return.elog boolean - if \code{TRUE} then the event log is returned in 
#'   addition to the CBS summary
#' @return list with elements \code{cbs} and \code{elog} containing data.frames
#' @export
bgcnbd.GenerateData <- function(n, T.cal, T.star=T.cal, params, 
                                return.elog = FALSE, dropout_at_zero = FALSE) {
  # check model parameters
  dc.check.model.params(c("k", "r", "alpha", "a", "b"), params,
                        "bgcnbd.GenerateData")
  if (params[1] != floor(params[1]) | params[1] < 1)
    stop("k must be integer being greater or equal to 1.")  
  
  k      <- params[1]
  r      <- params[2]
  alpha  <- params[3]
  a      <- params[4]
  b      <- params[5]

  if (length(T.cal)==1)  T.cal  <- rep(T.cal, n)  
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
    coins <- rbinom(min(10000, 10/p), 1, p)
    churn <- ifelse(any(coins==1), min(which(coins==1)), 10000)
    if (dropout_at_zero) churn <- churn - 1
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
