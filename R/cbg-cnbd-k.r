
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
#'   'cust' and transaction time 't' or 'date'
#' @param method 'wheat', 'mle', 'mle-minka', 'mle-thom', 'cv'
#' @param plot if TRUE then distribution of estimated regularity will be plotted
#' @return estimated real-valued regularity parameter; rounded to an integer,
#'   this can be used as 'k' for estimating CBG/CNBD-k models
#' @export
#' @seealso cbgcnbd.EstimateParameters
#' @example demo/cbg-cnbd-k.r
cbgcnbd.EstimateRegularity <- function(elog, method="wheat", plot=F) {
  if (!"cust" %in% names(elog)) 
    stop("Error in cbgcnbd.EstimateRegularity: elog must have a column labelled \"cust\"")
  if (!"date" %in% names(elog) & !"t" %in% names(elog)) 
    stop("Error in cbgcnbd.EstimateRegularity: elog must have a column labelled \"t\" or \"date\"")
  if (!"t" %in% names(elog)) elog$t <- as.numeric(elog$date)
  trans <- split(elog, elog$cust)
  if (method=="wheat")
  {
    # Wheat, Rita D., and Donald G. Morrison. 
    #   "Estimating purchase regularity with two interpurchase times." 
    #   Journal of Marketing Research (1990): 87-93.
    M <- unlist(lapply(trans, function(df) { 
      itt <- diff(sort(unique(df$t)))
      if (length(itt)>1) {
        # take last two itt's to capture most recent regularity
        itt2 <- rev(itt)[1:2]
        return(sample(itt2)[1]/sum(itt2))
      }
    }))
    r <- (1-4*var(M)) / (8*var(M))
    if (plot) {
      op <- par(mar=c(1, 2, 1, 2))
      plot(density(M), main="", sub="", xlab="", ylab="", lwd=2, frame=F, axes=F)
      polygon(density(M), col="lightgray", border=1)
      curve(dbeta(x, 1, 1), add=T, lty=2, lwd=2)
      curve(dbeta(x, round(r), round(r)), add=T, lty=2, lwd=2)
      par(op)
    }
    return(r)
  } 
  else 
  {
    if (method=="mle" | method=="mle-minka")
    {
      # Maximum Likelihood Estimator
      #  http://en.wikipedia.org/wiki/Gamma_distribution#Maximum_likelihood_estimation
      # Approximation for MLE by Minka
      #  http://research.microsoft.com/en-us/um/people/minka/papers/minka-gamma.pdf
      ks <- unlist(lapply(trans, function(df) {
        itt <- diff(sort(unique(df$t)))
        if (length(itt)>=9) {
          s <- log(sum(itt)/length(itt)) - sum(log(itt))/length(itt)
          if (method=="mle") {
            fn <- function(v) {
              return ((log(v)-digamma(v)-s)^2)
            }
            k <- optimize(fn, lower=.1, upper=50)$min
          } else if (method=="mle-minka") {
            k <- (3-s+sqrt((s-3)^2+24*s))/(12*s)
          }
          return(k)
        }
      }))
    }
    else if (method=="mle-thom")
    {
      # Approximation for ML estimator
      #   Thom (1968); see Dunn, Richard, Steven Reader, and Neil Wrigley. 
      #   "An investigation of the assumptions of the NBD model" Applied Statistics (1983): 249-259.
      ks <- unlist(lapply(trans, function(df) {
        itt <- diff(sort(unique(df$t)))
        if (length(itt)>=9) {
          hm <- function(v) exp(sum(log(v))/length(v))
          mu <- log(mean(itt) / hm(itt))
          d <- (1/(4*mu))*(1+sqrt(1+4*mu/3))
          return(d)
        }
      }))
    }
    else if (method=="cv")
    {
      # Estimate regularity by analyzing coefficient of variation
      #   Wu, Couchen, and H-L. Chen. "A consumer purchasing model with learning and departure behaviour." 
      #   Journal of the Operational Research Society (2000): 583-591.
      ks <- unlist(lapply(trans, function(df) {
        itt <- diff(sort(unique(df$t)))    
        if (length(itt)>=9) {
          cv <- sd(itt) / mean(itt)
          k <- 1 / cv^2
          return(k)
        }
      }))
    }
    if (length(ks)==0)
      stop("No customers with 10 or more transactions.")
    if (plot) {
      ymax <- median(ks)*3
      boxplot(ks, horizontal=T, ylim=c(0, ymax), frame=F, axes=F)
      axis(1, at=0:ymax)
      axis(3, at=0:ymax, labels=rep("", 1+ymax))
      abline(v=1:ymax, lty="dotted", col="lightgray")
      boxplot(ks, horizontal=T, add=T, col="gray", frame=F, axes=F)
    }
    return(median(ks))
  }
}

#' Parameter Estimation for the CBG/CNBD-k model
#' 
#' Estimates parameters for the CBG/CNBD-k model via Maximum Likelihood 
#' Estimation.
#' 
#' @param cal.cbs calibration period CBS. It must contain columns for frequency 
#'   'x', for recency 't.x.' and total time observed 'T.cal'. Optionally a 
#'   column 'custs' can be provided, which represents number of customers with a
#'   specific combination of frequency 'x', recency 't.x' and 'T.cal'.
#' @param par.start initial CBG/CNBD-k parameters - a vector with 'r', 'alpha', 
#'   'a' and 'b' in that order.
#' @param max.param.value the upper bound on parameters
#' @param k specified degree of regularity for Erlang-k distributed 
#'   interpurchase times; needs to be integer-value; if this is not specified, 
#'   then grid search from 1 to 12 is performed; this however requires column
#'   'litt' to be present in cal.cbs, which represents sum of logarithmic
#'   interpurchase times during calibration period;
#' @return list of estimated parameters
#' @import BTYD
#' @export
#' @seealso elog2cbs
#' @references Platzer, Michael. "Stochastic models of noncontractual consumer 
#'   relationships." Master of Science in Business Administration thesis, Vienna
#'   University of Economics and Business Administration, Austria (2008). 
#'   \url{https://sites.google.com/site/michaelplatzer/stochastic-models-of-noncontractual-consumer-relationships}
#' @example demo/cbg-cnbd-k.r
cbgcnbd.EstimateParameters <- function(cal.cbs, par.start = c(1, 1, 1, 1), max.param.value = 10000, k = NULL) {
  
  dc.check.model.params(c("r", "alpha", "a", "b"), par.start, 
    "cbgcnbd.EstimateParameters")
  
  # Either 'k' or 'litt' need to be present
  if (is.null(k) & !"litt" %in% colnames(cal.cbs))
    stop("Either regularity parameter k need to be specified, or a column with logarithmic interpurchase times 'litt' need to be present in cal.cbs")
  
  # if k is not specified we do grid search for k
  if (is.null(k)) {
    params <- list()
    LL <- c()
    for (k in 1:12) {
      params[[k]] <- cbgcnbd.EstimateParameters(cal.cbs, par.start, max.param.value, k=k)
      LL[k] <- cbgcnbd.cbs.LL(params[[k]], cal.cbs)
      if (k > 4 && LL[k] < LL[k-1] && LL[k-1] < LL[k-2]) break; # stop if LL gets worse for increasing k
    }
    k <- which.max(LL)
    return(params[[k]])
  }
  
  # if 'litt' is missing, we set it to zero, so that cbgcnbd.cbs.LL does not
  # complain; however this makes LL values for different k values not comparable
  if (!"litt" %in% colnames(cal.cbs))
    cal.cbs[, "litt"] <- 0 
  
  cbgcnbd.eLL <- function(params, k, cal.cbs, max.param.value) {
    params <- exp(params)
    params[params > max.param.value] <- max.param.value
    params <- c(k, params)
    return(-1 * cbgcnbd.cbs.LL(params, cal.cbs))
  }
  
  logparams <- log(par.start)
  results <- optim(logparams, cbgcnbd.eLL, cal.cbs = cal.cbs, k = k,
    max.param.value = max.param.value, method = "L-BFGS-B")
  estimated.params <- exp(results$par)
  estimated.params[estimated.params > max.param.value] <- max.param.value
  estimated.params <- c(k, estimated.params)
  names(estimated.params) <- c("k", "r", "alpha", "a", "b")
  return(estimated.params)
}


#' Calculate the log-likelihood of the CBG/CNBD-k model
#' 
#' @param params CBG/CNBD-k parameters - a vector with 'k', 'r', 'alpha', 'a' 
#'   and 'b' in that order.
#' @param cal.cbs calibration period CBS. It must contain columns for frequency 
#'   'x', for recency 't.x.', for sum of logarithmic interpurchase times 'litt'
#'   and total time observed 'T.cal'. Optionally a column 'custs' can be
#'   provided, which represents number of customers with a specific combination
#'   of frequency 'x' and 'T.cal'.
#' @return the total log-likelihood for the provided data.
#' @export
#' @seealso cbgcnbd.EstimateParameters 
cbgcnbd.cbs.LL <- function(params, cal.cbs) {
  dc.check.model.params(c("k", "r", "alpha", "a", "b"), params, 
   "cbgcnbd.cbs.LL")
  tryCatch(x <- cal.cbs[, "x"], error = function(e) stop("Error in cbgcnbd.cbs.LL: cal.cbs must have a frequency column labelled \"x\""))
  tryCatch(t.x <- cal.cbs[, "t.x"], error = function(e) stop("Error in cbgcnbd.cbs.LL: cal.cbs must have a recency column labelled \"t.x\""))
  tryCatch(T.cal <- cal.cbs[, "T.cal"], error = function(e) stop("Error in cbgcnbd.cbs.LL: cal.cbs must have a column for length of time observed labelled \"T.cal\""))
  tryCatch(litt <- cal.cbs[, "litt"], error = function(e) stop("Error in cbgcnbd.cbs.LL: cal.cbs must have a column for sum over logarithmic inter-transaction-times labelled \"litt\""))
  if ("custs" %in% colnames(cal.cbs)) {
    custs <- cal.cbs[, "custs"]
  } else {
    custs <- rep(1, length(x))
  }
  return(sum(custs * cbgcnbd.LL(params, x, t.x, T.cal, litt)))
}


#' Calculate the log-likelihood of the CBG/CNBD-k model
#' 
#' @param params CBG/CNBD-k parameters - a vector with 'k', 'r', 'alpha', 'a'
#'   and 'b' in that order.
#' @param x frequency, i.e. number of re-purchases
#' @param t.x recency, i.e. time elapsed from first purchase to last purchase
#' @param T.cal total time of observation period
#' @param litt sum of logarithmic interpurchase times
#' @return a vector of log-likelihoods
#' @export
#' @seealso cbgcnbd.EstimateParameters
cbgcnbd.LL <- function(params, x, t.x, T.cal, litt) {
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
   "cbgcnbd.LL")
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


#' CBG/CNBD-k P(alive)
#' 
#' Uses CBG/CNBD-k model parameters and a customer's past transaction behavior 
#' to return the probability that they are still alive at the end of the 
#' calibration period.
#' 
#' @param params CBG/CNBD-k parameters - a vector with 'k', 'r', 'alpha', 'a'
#'   and 'b' in that order.
#' @param x number of repeat transactions in the calibration period T.cal, or a 
#'   vector of calibration period frequencies.
#' @param t.x recency, i.e. length between first and last transaction during 
#'   calibration period.
#' @param T.cal length of calibration period, or a vector of calibration period 
#'   lengths.
#' @return Probability that the customer is still alive at the end of the 
#'   calibration period.
#' @export
#' @example demo/cbg-cnbd-k.r
#' @seealso cbgcnbd.EstimateParameters
cbgcnbd.PAlive <- function(params, x, t.x, T.cal) {
  max.length <- max(length(x), length(t.x), length(T.cal))
  if (max.length%%length(x)) 
    warning("Maximum vector length not a multiple of the length of x")
  if (max.length%%length(t.x)) 
    warning("Maximum vector length not a multiple of the length of t.x")
  if (max.length%%length(T.cal)) 
    warning("Maximum vector length not a multiple of the length of T.cal")
  dc.check.model.params(c("k", "r", "alpha", "a", "b"), params, 
   "cbgcnbd.PAlive")  
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


#' CBD/CNBD-k Conditional Expected Transactions
#' 
#' Uses CBG/CNBD-k model parameters and a customer's past transaction behavior
#' to return the number of transactions they are expected to make in a given
#' time period.
#' 
#' @param params CBG/CNBD-k parameters - a vector with 'k', 'r', 'alpha', 'a'
#'   and 'b' in that order.
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
#' @example demo/cbg-cnbd-k.r
#' @seealso cbgcnbd.EstimateParameters
cbgcnbd.ConditionalExpectedTransactions <- function(params, T.star, x, t.x, T.cal) {
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
    "cbgcnbd.ConditionalExpectedTransactions")
  if (params[1] != floor(params[1]) | params[1] < 1)
    stop("k must be integer being greater or equal to 1.")
  if (params[1]>1)
    message("Results for k>1 are approximative")
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
  # for k>1 we simply scale rate parameter alpha accordingly to approximate results
  alpha <- alpha * k
  P1 <- ((a + b + x) / (a - 1))
  P2 <- 1 - ((alpha+T.cal)/(alpha+T.cal+T.star))^(r+x) * hyperg_2F1(r+x, b+x+1, b+x+a, T.star/(alpha+T.cal+T.star))
  P3 <- cbgcnbd.PAlive(params, x, t.x, T.cal)
  return (P1 * P2 * P3)
}


#' CBD/CNBD-k Unconditional Expected Transactions
#' 
#' Uses CBG/CNBD-k model parameters to return the number of transactions a 
#' random customer is expected to make in a given time period.
#' 
#' @param params CBG/CNBD-k parameters - a vector with 'k', 'r', 'alpha', 'a'
#'   and 'b' in that order.
#' @param t length of time for which we are calculating the expected number of 
#'   transactions.
#' @return Number of transactions a customer is expected to make in a time 
#'   period of length t. If any of the input parameters has a length greater 
#'   than 1, this will be a vector of expected number of transactions.
#' @export
#' @seealso cbgcnbd.EstimateParameters
cbgcnbd.UnconditionalExpectedTransactions <- function(params, t) {
  dc.check.model.params(c("k", "r", "alpha", "a", "b"), params, 
    "cbgcnbd.UnconditionalExpectedTransactions")
  if (params[1] != floor(params[1]) | params[1] < 1)
    stop("k must be integer being greater or equal to 1.")
  if (params[1]>1)
    message("Results for k>1 are approximative")
  k <- params[1]
  r <- params[2]
  alpha <- params[3]
  a <- params[4]
  b <- params[5]
  # for k>1 we simply scale rate parameter alpha accordingly to approximate results
  alpha <- alpha * k
  P1 <- ((a + b + x) / (a - 1))
  P2 <- 1 - ((alpha+T.cal)/(alpha+T.cal+T.star))^(r+x) * hyperg_2F1(r+x, b+x+1, b+x+a, T.star/(alpha+T.cal+T.star))
  P3 <- cbgcnbd.PAlive(params, x, t.x, T.cal)
  return (P1 * P2 * P3)
}


#' CBD/CNBD-k Unconditional Probability Distribution of Transactions
#' 
#' Uses CBG/CNBD-k model parameters to return the probability distribution of
#' purchase frequencies for a random customer in a given time period, i.e.
#' P(X(t)=x|r,alpha,a,b)
#' 
#' @param params CBG/CNBD-k parameters - a vector with 'k', 'r', 'alpha', 'a' and 'b'
#'   in that order.
#' @param t length of time for which we are calculating the expected number of 
#'   transactions.
#' @param x number of transactions for which probability is calculated.
#' @return P(X(t)=x|r,alpha,a,b). If any of the input parameters has a length
#'   greater than 1, this will be a vector of expected number of transactions.
#' @export
#' @seealso cbgcnbd.EstimateParameters
cbgcnbd.Px <- function(params, t, x) {
  dc.check.model.params(c("k", "r", "alpha", "a", "b"), params, 
                        "cbgcnbd.Px")
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


#' Simulate data according to CBG/CNBD-k model assumptions
#' 
#' @param n number of customers
#' @param T.cal length of calibration period
#' @param T.star length of holdout period
#' @param params CBG/CNBD-k parameters - a vector with 'k', 'r', 'alpha', 'a' and 'b'
#'   in that order.
#' @param return.elog boolean - if TRUE then the event log is returned in 
#'   addition to the CBS summary
#' @return list with elements 'cbs' and 'elog' containing data.frames
#' @export
cbgcnbd.GenerateData <- function(n, T.cal, T.star=T.cal, params, return.elog=F) {
  # check model parameters
  dc.check.model.params(c("k", "r", "alpha", "a", "b"), params,
                        "cbgcnbd.GenerateData")
  if (params[1] != floor(params[1]) | params[1] < 1)
    stop("k must be integer being greater or equal to 1.")  
  
  k <- params[1]
  r <- params[2]
  alpha <- params[3]
  a <- params[4]
  b <- params[5]  
  
  # sample intertransaction timings parameter lambda for each customer
  lambdas <- rgamma(n, shape=r, rate=alpha)
  
  # sample churn-probability p for each customer
  ps <- rbeta(n, a, b)
  
  # sample intertransaction timings & churn
  cbs <- data.frame()
  elog <- data.frame(cust=numeric(0), t=numeric(0))
  for (i in 1:n) {
    p <- ps[i]
    lambda <- lambdas[i]
    # sample no. of transactions until churn
    churn <- which.max(rbinom(min(10000, 10/p), 1, p)) - 1
    # sample transaction times
    times <- cumsum(c(0, rgamma(churn, shape=k, rate=lambda)))
    if (return.elog) elog <- rbind(elog, data.frame(cust=i, t=times[times<(T.cal+T.star)]))
    # determine frequency, recency, etc.
    ts.cal <- times[times<T.cal]
    ts.star <- times[times>=T.cal & times<(T.cal+T.star)]
    cbs[i, "x"] <- length(ts.cal)-1
    cbs[i, "t.x"] <- max(ts.cal)
    cbs[i, "T.cal"] <- T.cal
    cbs[i, "litt"] <- sum(log(diff(ts.cal)))
    cbs[i, "alive"] <- churn > cbs[i, "x"]
    cbs[i, "T.star"] <- T.star
    cbs[i, "x.star"] <- length(ts.star)
    cbs[i, "p"] <- p
    cbs[i, "lambda"] <- lambda
    cbs[i, "churn"] <- churn
    cbs[i, "k"] <- k
    cbs[i, "x.max"] <- max(times)
  }
  out <- list(cbs=cbs)
  if (return.elog) out$elog <- elog
  return(out)
}
