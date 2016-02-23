
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
  
  # either `k` or `litt` need to be present
  if (is.null(k) & !"litt" %in% colnames(cal.cbs))
    stop("Either regularity parameter k need to be specified, or a column with logarithmic interpurchase times litt need to be present in cal.cbs")
  
  # if `k` is not specified we do grid search for `k`
  if (is.null(k)) {
    params <- list()
    LL <- c()
    for (k in 1:12) {
      params[[k]] <- tryCatch(bgcnbd.EstimateParameters(cal.cbs = cal.cbs, k = k, par.start = par.start,
                                                        max.param.value = max.param.value, trace = trace,
                                                        dropout_at_zero = dropout_at_zero),
                              error = function(e) { e })
      if (inherits(params[[k]], "error")) {
        params[[k]] <- NULL
        break # stop if parameters could not be estimated, e.g. if bgcnbd.LL returns Inf
      }
      LL[k] <- bgcnbd.cbs.LL(params = params[[k]], cal.cbs = cal.cbs, 
                             dropout_at_zero = dropout_at_zero)
      if (k > 4 && LL[k] < LL[k-1] && LL[k-1] < LL[k-2]) 
        break # stop if LL gets worse for increasing k
    }
    k <- which.max(LL)
    return(params[[k]])
  }
  
  # if `litt` is missing, we set it to zero, so that bgcnbd.cbs.LL does not
  # complain; however this makes LL values for different k values not comparable
  if (!"litt" %in% colnames(cal.cbs))
    cal.cbs[, "litt"] <- 0 
  
  count <- 0
  bgcnbd.eLL <- function(params, k, cal.cbs, max.param.value, dropout_at_zero) {
    params <- exp(params)
    params[params > max.param.value] <- max.param.value
    params <- c(k, params)
    loglik <- bgcnbd.cbs.LL(params = params, cal.cbs = cal.cbs, dropout_at_zero = dropout_at_zero)
    count <<- count + 1
    if (trace>0 & count%%trace==0)
      cat("bgcnbd.EstimateParameters - iter", count, ":", sprintf("%12.2f", loglik), ":", sprintf("%10.6f", params), "\n")
    return(-1 * loglik)
  }
  
  logparams <- log(par.start)
  results <- optim(logparams, bgcnbd.eLL, cal.cbs = cal.cbs, k = k,
    max.param.value = max.param.value, dropout_at_zero = dropout_at_zero, 
    method = "L-BFGS-B")
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
#' @references Platzer Michael, and Thomas Reutterer (forthcoming)
bgcnbd.cbs.LL <- function(params, cal.cbs, dropout_at_zero = FALSE) {
  dc.check.model.params(c("k", "r", "alpha", "a", "b"), params, "bgcnbd.cbs.LL")
  tryCatch(x <- cal.cbs$x, error = function(e) stop("Error in bgcnbd.cbs.LL: cal.cbs must have a frequency column labelled \"x\""))
  tryCatch(t.x <- cal.cbs$t.x, error = function(e) stop("Error in bgcnbd.cbs.LL: cal.cbs must have a recency column labelled \"t.x\""))
  tryCatch(T.cal <- cal.cbs$T.cal, error = function(e) stop("Error in bgcnbd.cbs.LL: cal.cbs must have a column for length of time observed labelled \"T.cal\""))
  tryCatch(litt <- cal.cbs$litt, error = function(e) stop("Error in bgcnbd.cbs.LL: cal.cbs must have a column for sum over logarithmic inter-transaction-times labelled \"litt\""))
  if ("custs" %in% colnames(cal.cbs)) {
    custs <- cal.cbs$custs
  } else {
    custs <- rep(1, length(x))
  }
  return(sum(custs * bgcnbd.LL(params = params, 
                               x = x, t.x = t.x, T.cal = T.cal, litt = litt, 
                               dropout_at_zero = dropout_at_zero)))
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
#' @references Platzer Michael, and Thomas Reutterer (forthcoming)
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
  dc.check.model.params(c("k", "r", "alpha", "a", "b"), params, "bgcnbd.LL")
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
  S1 <- as.numeric(dropout_at_zero | x>0) * a / (b + x - 1 + ifelse(dropout_at_zero, 1, 0)) * ((alpha+T.cal)/(alpha+t.x))^(r+k*x)
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
#' @references Platzer Michael, and Thomas Reutterer (forthcoming)
bgcnbd.pmf <- function(params, t, x, dropout_at_zero = FALSE) {

  dc.check.model.params(c("k", "r", "alpha", "a", "b"), params, "bgcnbd.pmf")
  if (params[1] != floor(params[1]) | params[1] < 1)
    stop("k must be integer being greater or equal to 1.")
  if (length(t) > 1 & length(x) > 1 & length(t) != length(x))
    stop("parameters t and x must be of same length")

  max.length <- max(length(t), length(x))
  t <- rep(t, length.out = max.length)
  x <- rep(x, length.out = max.length)
  res <- sapply(1:max.length, function(i) {
    # call C++ implementation
    bgcnbd_pmf_cpp(params, t[i], x[i], dropout_at_zero)
  })
  return(res)
}


#' BG/CNBD-k Expectation
#' 
#' Returns the number of repeat transactions that a randomly chosen customer
#' (for whom we have no prior information) is expected to make in a given time
#' period.
#' 
#' E(X(t) | k, r, alpha, a, b)
#' 
#' @param params BG/CNBD-k parameters - a vector with \code{k}, \code{r}, \code{alpha}, \code{a} and \code{b}
#'   in that order.
#' @param t length of time for which we are calculating the expected number of repeat transactions.
#' @param dropout_at_zero Boolean; the mbg-methods are simple wrapper methods,
#'   which set this parameter to TRUE
#' @return Number of repeat transactions a customer is expected to make in a time period of length t.
#' @references Platzer Michael, and Thomas Reutterer (forthcoming)
#' @seealso \code{\link{bgcnbd.EstimateParameters}}
bgcnbd.Expectation <- function(params, t, dropout_at_zero = FALSE) {
  
  dc.check.model.params(c("k", "r", "alpha", "a", "b"), params, "bgcnbd.Expectation")
  if (any(t < 0) || !is.numeric(t)) 
    stop("t must be numeric and may not contain negative numbers.")
  
  # to save computation time, we collapse vector `t` on to its unique values
  ts <- unique(t); names(ts) <- ts
  ts_map <- sapply(ts, function(t) {
    # call C++ implementation
    bgcnbd_exp_cpp(params, t, dropout_at_zero)
  })
  res <- unname(ts_map[as.character(t)])
  return(res)
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
#' @references Platzer Michael, and Thomas Reutterer (forthcoming)
bgcnbd.PAlive <- function(params, x, t.x, T.cal, dropout_at_zero = FALSE) {
  max.length <- max(length(x), length(t.x), length(T.cal))
  if (max.length%%length(x)) 
    warning("Maximum vector length not a multiple of the length of x")
  if (max.length%%length(t.x)) 
    warning("Maximum vector length not a multiple of the length of t.x")
  if (max.length%%length(T.cal)) 
    warning("Maximum vector length not a multiple of the length of T.cal")
  dc.check.model.params(c("k", "r", "alpha", "a", "b"), params, "bgcnbd.PAlive")  
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
  return(palive)
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
#' @references Platzer Michael, and Thomas Reutterer (forthcoming)
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
  dc.check.model.params(c("k", "r", "alpha", "a", "b"), params, "bgcnbd.ConditionalExpectedTransactions")
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
  if (round(a, 2)==1) a <- a + 0.01 # P1 not defined for a=1, so we add slight noise in such rare cases
  if (k>1) cat("note: conditional expected transactions can only be approximated for k>1\n")
  # approximate via expression for conditional expected transactions for BG/NBD model, but adjust scale parameter by k
  G <- function(r, alpha, a, b) 1 - (alpha/(alpha+T.star))^r * gsl::hyperg_2F1(r, b+1, a+b, T.star/(alpha+T.star))
  P1 <- (a+b+x-1+ifelse(dropout_at_zero, 1, 0)) / (a-1)
  P2 <- G(r+x, k*alpha+T.cal, a, b+x-1+ifelse(dropout_at_zero, 1, 0))
  P3 <- bgcnbd.PAlive(params = params, x = x, t.x = t.x, T.cal = T.cal, dropout_at_zero = dropout_at_zero)
  exp <- P1 * P2 * P3
  # adjust bias BG/NBD-based approximation by scaling via the Unconditional Expectations (for wich we have exact expression)
  if (k>1) {
    sum.cal <- sum(bgcnbd.Expectation(params = params, t = T.cal,        dropout_at_zero = dropout_at_zero))
    sum.tot <- sum(bgcnbd.Expectation(params = params, t = T.cal+T.star, dropout_at_zero = dropout_at_zero))
    exp <- exp * (sum.tot-sum.cal) / sum(exp)
  }
  return(exp)
}


#' Simulate data according to BG/CNBD-k model assumptions
#' 
#' @param n number of customers
#' @param T.cal length of calibration period; if vector then it is assumed that
#'   customers have different "birth" dates, i.e. max(T.cal)-T.cal
#' @param T.star length(s) of holdout period(s); assumed to be same for all
#'   customers
#' @param params BG/CNBD-k parameters - a vector with \code{k}, \code{r},
#'   \code{alpha}, \code{a} and \code{b} in that order.
#' @param dropout_at_zero Boolean; the mbg-methods are simple wrapper methods,
#'   which set this parameter to TRUE
#' @param return.elog boolean - if \code{TRUE} then the event log is returned in 
#'   addition to the CBS summary
#' @return list with elements \code{cbs} and \code{elog} containing data.frames
#' @export
#' @seealso \code{\link{bgcnbd.EstimateParameters}}
#' @references Platzer Michael, and Thomas Reutterer (forthcoming)
bgcnbd.GenerateData <- function(n, T.cal, T.star=NULL, params, 
                                return.elog = FALSE, dropout_at_zero = FALSE) {
  dc.check.model.params(c("k", "r", "alpha", "a", "b"), params, "bgcnbd.GenerateData")
  if (params[1] != floor(params[1]) | params[1] < 1)
    stop("k must be integer being greater or equal to 1.")  
  
  k      <- params[1]
  r      <- params[2]
  alpha  <- params[3]
  a      <- params[4]
  b      <- params[5]

  if (length(T.cal)==1)  T.cal  <- rep(T.cal, n)  
  
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
    times <- cumsum(c(max(T.cal)-T.cal[i], rgamma(churn, shape=k, rate=lambda)))
    if (return.elog)
      elog_list[[i]] <- data.frame(cust=i, t=times[times<(max(T.cal)+max(T.star))])
    # determine frequency, recency, etc.
    ts.cal  <- times[times < max(T.cal)]
    ts.star <- times[times >= max(T.cal) & times<(max(T.cal)+T.star[i])]
    cbs_list[[i]] <- list(cust   = i,
                          x      = length(ts.cal)-1,
                          t.x    = max(ts.cal) - (max(T.cal)-T.cal[i]),
                          litt   = sum(log(diff(ts.cal))),
                          churn  = churn,
                          alive  = churn > (length(ts.cal)-1))
    for (tstar in T.star) {
      colname <- paste0('x.star', ifelse(length(T.star)>1, tstar, ''))
      cbs_list[[i]][[colname]] <- length(times[times >= max(T.cal) & times<(max(T.cal)+tstar)])
    }
  }
  cbs <- do.call(rbind.data.frame, cbs_list)
  cbs$lambda <- lambdas
  cbs$p      <- ps
  cbs$k      <- k
  cbs$T.cal  <- T.cal
  if (length(T.star)==1) cbs$T.star <- T.star
  rownames(cbs) <- NULL
  out <- list(cbs=cbs)
  if (return.elog) {
    elog <- do.call(rbind.data.frame, elog_list)
    out$elog <- elog
  }
  return(out)
}


#' BG/CNBD-k Plot Frequency in Calibration Period
#' 
#' Plots a histogram and returns a matrix comparing the actual and expected
#' number of customers who made a certain number of repeat transactions in the
#' calibration period, binned according to calibration period frequencies.
#' 
#' @param params BG/CNBD-k parameters - a vector with \code{k}, \code{r}, \code{alpha}, \code{a} and \code{b}
#'   in that order.
#' @param cal.cbs calibration period CBS (customer by sufficient statistic). It
#'   must contain columns for frequency ("x") and total time observed ("T.cal").
#' @param censor integer used to censor the data.
#' @param plotZero If FALSE, the histogram will exclude the zero bin.
#' @param xlab descriptive label for the x axis.
#' @param ylab descriptive label for the y axis.
#' @param title title placed on the top-center of the plot.
#' @param dropout_at_zero Boolean; the mbg-methods are simple wrapper methods,
#'   which set this parameter to TRUE
#' @return Calibration period repeat transaction frequency comparison matrix (actual vs. expected).
#'
#' @seealso \code{\link{bgcnbd.EstimateParameters}}
#' @references Platzer Michael, and Thomas Reutterer (forthcoming)
#' @examples 
#' set.seed(1)
#' params <- c(k=3, r=0.85, alpha=1.45, a=0.79, b=2.42)
#' cal.cbs <- bgcnbd.GenerateData(n=5000, 32, 0, params)$cbs
#' bgcnbd.PlotFrequencyInCalibration(params, cal.cbs, censor = 7)
bgcnbd.PlotFrequencyInCalibration <- function(params, cal.cbs, 
                                              censor = NULL, 
                                              plotZero = TRUE, 
                                              xlab = "Calibration period transactions",
                                              ylab = "Customers", 
                                              title = "Frequency of Repeat Transactions",
                                              dropout_at_zero = FALSE) {
  tryCatch(x <- cal.cbs$x, error = function(e) stop("Error in bgcnbd.PlotFrequencyInCalibration: cal.cbs must have a frequency column labelled \"x\""))
  tryCatch(T.cal <- cal.cbs$T.cal, error = function(e) stop("Error in bgcnbd.PlotFrequencyInCalibration: cal.cbs must have a column for length of time observed labelled \"T.cal\""))
  dc.check.model.params(c("k", "r", "alpha", "a", "b"), params, "bgcnbd.PlotFrequencyInCalibration")
  if (is.null(censor) || censor > max(x)) censor <- max(x)
  n.x    <- rep(0, max(x) + 1)
  for (xx in unique(x)) {
      n.x[xx + 1] <- sum(xx == x)
  }
  n.x.censor <- sum(n.x[(censor + 1):length(n.x)])
  n.x.actual <- c(n.x[1:censor], n.x.censor)
  T.value.counts <- table(T.cal)
  T.values <- as.numeric(names(T.value.counts))
  n.T.values <- length(T.values)
  n.x.expected <- rep(0, length(n.x.actual))
  n.x.expected.all <- rep(0, max(x) + 1)
  for (xx in 0:max(x)) {
    this.x.expected <- 0
    for (T.idx in 1:n.T.values) {
      t <- T.values[T.idx]
      if (t == 0) 
          next
      n.T <- T.value.counts[T.idx]
      prob.of.this.x.for.this.T <- bgcnbd.pmf(params, t, xx, dropout_at_zero=dropout_at_zero)
      expected.given.x.and.T <- n.T * prob.of.this.x.for.this.T
      this.x.expected <- this.x.expected + expected.given.x.and.T
    }
    n.x.expected.all[xx + 1] <- this.x.expected
  }
  n.x.expected[1:censor] <- n.x.expected.all[1:censor]
  n.x.expected[censor + 1] <- sum(n.x.expected.all[(censor + 1):(max(x) + 1)])
  col.names <- paste(rep("freq", length(censor + 1)), (0:censor), sep = ".")
  col.names[censor + 1] <- paste(col.names[censor + 1], "+", sep = "")
  censored.freq.comparison <- rbind(n.x.actual, n.x.expected)
  colnames(censored.freq.comparison) <- col.names
  cfc.plot <- censored.freq.comparison
  if (plotZero == FALSE) {
    cfc.plot <- cfc.plot[, -1]
  }
  n.ticks <- ncol(cfc.plot)
  if (plotZero == TRUE) {
    x.labels <- 0:(n.ticks - 1)
    x.labels[n.ticks] <- paste(n.ticks - 1, "+", sep = "")
  }
  ylim <- c(0, ceiling(max(cfc.plot) * 1.1))
  barplot(cfc.plot, names.arg = x.labels, beside = TRUE, ylim = ylim, 
      main = title, xlab = xlab, ylab = ylab, col = 1:2)
  model <- paste0(ifelse(dropout_at_zero, "MBG", "BG"), "/",
                  ifelse(params[1]>1, paste0("CNBD-", params[1]), "NBD"))
  legend("topright", legend = c("Actual", model), col = 1:2, 
      lwd = 2)
  return(censored.freq.comparison)
}


#' BG/CNBD-k Expected Cumulative Transactions
#' 
#' Calculates the expected cumulative total repeat transactions by all customers
#' for the calibration and holdout periods.
#' 
#' @param params BG/CNBD-k parameters - a vector with \code{k}, \code{r}, \code{alpha}, \code{a} and \code{b}
#'   in that order.
#' @param T.cal a vector to represent customers' calibration period lengths (in
#'   other words, the "T.cal" column from a customer-by-sufficient-statistic
#'   matrix).
#' @param T.tot end of holdout period. Must be a single value, not a vector.
#' @param n.periods.final number of time periods in the calibration and holdout periods.
#' @param dropout_at_zero Boolean; the mbg-methods are simple wrapper methods,
#'   which set this parameter to TRUE
#' @return Vector of expected cumulative total repeat transactions by all customers.
#' 
#' @seealso \code{\link{bgnbd.ExpectedCumulativeTransactions}}
bgcnbd.ExpectedCumulativeTransactions <- function(params, T.cal, T.tot, n.periods.final, 
                                                  dropout_at_zero = FALSE) {
  dc.check.model.params(c("k", "r", "alpha", "s", "beta"), params, "bgcnbd.ExpectedCumulativeTransactions")
  if (any(T.cal < 0) || !is.numeric(T.cal))
    stop("T.cal must be numeric and may not contain negative numbers.")
  if (length(T.tot) > 1 || T.tot < 0 || !is.numeric(T.tot))
    stop("T.cal must be a single numeric value and may not be negative.")
  if (length(n.periods.final) > 1 || n.periods.final < 0 || !is.numeric(n.periods.final)) 
    stop("n.periods.final must be a single numeric value and may not be negative.")
  
  intervals <- seq(T.tot/n.periods.final, T.tot, length.out = n.periods.final)
  cust.birth.periods <- max(T.cal) - T.cal
  expected.transactions <- sapply(intervals, function(interval) {
    if (interval <= min(cust.birth.periods)) 
      return(0)
    sum(bgcnbd.Expectation(params = params, 
                           t = interval - cust.birth.periods[cust.birth.periods < interval], 
                           dropout_at_zero = dropout_at_zero))
  })
  return(expected.transactions)
}


#' BG/CNBD-k Tracking Cumulative Transactions Plot
#' 
#' Plots the actual and expected cumulative total repeat transactions by all
#' customers for the calibration and holdout periods, and returns this
#' comparison in a matrix.
#' 
#' @param params BG/CNBD-k parameters - a vector with \code{k}, \code{r}, \code{alpha}, \code{a} and \code{b}
#'   in that order.
#' @param T.cal a vector to represent customers' calibration period lengths (in
#'   other words, the "T.cal" column from a customer-by-sufficient-statistic
#'   matrix).
#' @param T.tot end of holdout period. Must be a single value, not a vector.
#' @param actual.cu.tracking.data vector containing the cumulative number of
#'   repeat transactions made by customers for each period in the total time
#'   period (both calibration and holdout periods).
#' @param xlab descriptive label for the x axis.
#' @param ylab descriptive label for the y axis.
#' @param xticklab vector containing a label for each tick mark on the x axis.
#' @param title title placed on the top-center of the plot.
#' @param ymax upper boundary for y axis.
#' @param dropout_at_zero Boolean; the mbg-methods are simple wrapper methods,
#'   which set this parameter to TRUE
#' @return Matrix containing actual and expected cumulative repeat transactions.
#' 
#' @seealso \code{\link{bgnbd.PlotTrackingCum}}
bgcnbd.PlotTrackingCum <- function(params, T.cal, T.tot, actual.cu.tracking.data, 
    xlab = "Week", ylab = "Cumulative Transactions", 
    xticklab = NULL, title = "Tracking Cumulative Transactions", 
    ymax = NULL, dropout_at_zero = FALSE) {
  
  dc.check.model.params(c("k", "r", "alpha", "a", "b"), params, "bgcnbd.PlotTrackingCum")
  if (any(T.cal < 0) || !is.numeric(T.cal)) 
    stop("T.cal must be numeric and may not contain negative numbers.")
  if (any(actual.cu.tracking.data < 0) || !is.numeric(actual.cu.tracking.data)) 
    stop("actual.cu.tracking.data must be numeric and may not contain negative numbers.")
  if (length(T.tot) > 1 || T.tot < 0 || !is.numeric(T.tot)) 
    stop("T.cal must be a single numeric value and may not be negative.")
  
  actual <- actual.cu.tracking.data
  expected <- bgcnbd.ExpectedCumulativeTransactions(params, T.cal, T.tot, length(actual), dropout_at_zero = dropout_at_zero)
  cu.tracking.comparison <- rbind(actual, expected)
  if (is.null(ymax)) ymax <- max(c(actual, expected)) * 1.05
  plot(actual, type = "l", xaxt = "n", xlab = xlab, ylab = ylab, 
       col = 1, ylim = c(0, ymax), main = title)
  lines(expected, lty = 2, col = 2)
  if (is.null(xticklab) == FALSE) {
    if (ncol(cu.tracking.comparison) != length(xticklab)) {
      stop("Plot error, xticklab does not have the correct size")
    }
    axis(1, at = 1:ncol(cu.tracking.comparison), labels = xticklab)
  } else {
    axis(1, at = 1:length(actual), labels = 1:length(actual))
  }
  abline(v = max(T.cal), lty = 2)
  model <- paste0(ifelse(dropout_at_zero, "MBG", "BG"), "/",
                  ifelse(params[1]>1, paste0("CNBD-", params[1]), "NBD"))
  legend("bottomright", legend = c("Actual", model), col = 1:2, lty = 1:2, lwd = 1)
  return(cu.tracking.comparison)
}


#' BG/CNBD-k Tracking Incremental Transactions Comparison
#' 
#' Plots the actual and expected incremental total repeat transactions by all
#' customers for the calibration and holdout periods, and returns this
#' comparison in a matrix.
#' 
#' @param params BG/CNBD-k parameters - a vector with \code{k}, \code{r}, \code{alpha}, \code{a} and \code{b}
#'   in that order.
#' @param T.cal a vector to represent customers' calibration period lengths (in
#'   other words, the "T.cal" column from a customer-by-sufficient-statistic
#'   matrix).
#' @param T.tot end of holdout period. Must be a single value, not a vector.
#' @param actual.cu.tracking.data vector containing the cumulative number of
#'   repeat transactions made by customers for each period in the total time
#'   period (both calibration and holdout periods).
#' @param xlab descriptive label for the x axis.
#' @param ylab descriptive label for the y axis.
#' @param xticklab vector containing a label for each tick mark on the x axis.
#' @param title title placed on the top-center of the plot.
#' @param ymax upper boundary for y axis.
#' @param dropout_at_zero Boolean; the mbg-methods are simple wrapper methods,
#'   which set this parameter to TRUE
#' @return Matrix containing actual and expected incremental repeat transactions.
#' 
#' @seealso \code{\link{bgnbd.PlotTrackingInc}}
bgcnbd.PlotTrackingInc <- function(params, T.cal, T.tot, actual.inc.tracking.data, 
    xlab = "Week", ylab = "Transactions",
    xticklab = NULL, title = "Tracking Weekly Transactions", 
    ymax = NULL, dropout_at_zero = FALSE) {
  
  dc.check.model.params(c("k", "r", "alpha", "a", "b"), params, "bgcnbd.PlotTrackingInc")
  if (any(T.cal < 0) || !is.numeric(T.cal)) 
   stop("T.cal must be numeric and may not contain negative numbers.")
  if (any(actual.inc.tracking.data < 0) || !is.numeric(actual.inc.tracking.data)) 
    stop("actual.inc.tracking.data must be numeric and may not contain negative numbers.")
  if (length(T.tot) > 1 || T.tot < 0 || !is.numeric(T.tot)) 
    stop("T.cal must be a single numeric value and may not be negative.")
  
  actual <- actual.inc.tracking.data
  expected <- BTYD::dc.CumulativeToIncremental(bgcnbd.ExpectedCumulativeTransactions(params, T.cal, T.tot, length(actual), dropout_at_zero=dropout_at_zero))
  
  if (is.null(ymax)) ymax <- max(c(actual, expected)) * 1.05
  plot(actual, type = "l", xaxt = "n", xlab = xlab, ylab = ylab, 
       col = 1, ylim = c(0, ymax), main = title)
  lines(expected, lty = 2, col = 2)
  if (is.null(xticklab) == FALSE) {
    if (length(actual) != length(xticklab)) {
      stop("Plot error, xticklab does not have the correct size")
    }
    axis(1, at = 1:length(actual), labels = xticklab)
  } else {
    axis(1, at = 1:length(actual), labels = 1:length(actual))
  }
  abline(v = max(T.cal), lty = 2)
  
  model <- paste0(ifelse(dropout_at_zero, "MBG", "BG"), "/",
                  ifelse(params[1]>1, paste0("CNBD-", params[1]), "NBD"))
  legend("topright", legend = c("Actual", model), col = 1:2, lty = 1:2, lwd = 1)
  return(rbind(actual, expected))
}
