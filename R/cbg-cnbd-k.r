
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
#' @param k specified degree of regularity for Erlang-k distributed 
#'   interpurchase times
#' @param par.start initial CBG/CNBD-k parameters - a vector with 'r', 'alpha', 
#'   'a' and 'b' in that order.
#' @param max.param.value the upper bound on parameters
#' @return list of estimated parameters
#' @import BTYD
#' @export
#' @references Platzer, Michael. "Stochastic models of noncontractual consumer
#'   relationships." Master of Science in Business Administration thesis, Vienna
#'   University of Economics and Business Administration, Austria (2008).
#'   \url{https://sites.google.com/site/michaelplatzer/stochastic-models-of-noncontractual-consumer-relationships}
#' @examples
#' # generate artificial data 
#' k <- 3
#' params <- list(r=0.53, alpha=1.25, a=1.29, b=1.85)
#' data <- cbgcnbd.GenerateData(n=1000, k=k, T.cal=32, T.star=32, params=params, return.elog=T)
#' cbs <- data$cbs
#' elog <- data$elog
#' 
#' # estimate parameters, and compare to true parameters
#' k.est <- cbgcnbd.EstimateRegularity(elog)
#' # 3.175759; -> CBG/CNBD-3 should be modelled
#' est <- cbgcnbd.EstimateParameters(cbs, round(k.est))
#' rbind(params, est=round(est, 3))
#' #        r     alpha a     b    
#' # params 0.53  1.25  1.29  1.85 
#' # est    0.443 1.217 1.324 2.695
#' 
#' # estimate future transactions, and compare to actuals from holdout period
#' cbs$x.est <- cbgcnbd.ConditionalExpectedTransactions(est, k=round(k.est), cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)
#' sqrt(mean((cbs$x.star-cbs$x.est)^2))
#' # 1.168548
cbgcnbd.EstimateParameters <- function(cal.cbs, k, par.start = c(1, k, 1, 1), max.param.value = 10000) {
  dc.check.model.params(c("r", "alpha", "a", "b"), par.start, 
    "cbgcnbd.EstimateParameters")
  if (!"litt" %in% colnames(cal.cbs)) cal.cbs[, "litt"] <- 0 # for finding optimal parameters the litt-variable is not mandatory
  cbgcnbd.eLL <- function(params, k, cal.cbs, max.param.value) {
    params <- exp(params)
    params[params > max.param.value] <- max.param.value
    return(-1 * cbgcnbd.cbs.LL(params, k, cal.cbs))
  }
  logparams <- log(par.start)
  results <- optim(logparams, cbgcnbd.eLL, cal.cbs = cal.cbs, k = k,
    max.param.value = max.param.value, method = "L-BFGS-B")
  estimated.params <- exp(results$par)
  estimated.params[estimated.params > max.param.value] <- max.param.value
  return(estimated.params)
}


#' Calculate the log-likelihood of the CBG/CNBD-k model
#' 
#' @param params CBG/CNBD-k parameters - a vector with 'r', 'alpha', 'a' and 'b'
#'   in that order.
#' @param k specified degree of regularity for Erlang-k distributed
#'   interpurchase times
#' @param cal.cbs calibration period CBS. It must contain columns for frequency 
#'   'x', for recency 't.x.' and total time observed 'T.cal'. Optionally a 
#'   column 'custs' can be provided, which represents number of customers with a
#'   specific combination of frequency 'x' and 'T.cal'.
#' @return the total log-likelihood for the provided data.
#' @export
cbgcnbd.cbs.LL <- function(params, k, cal.cbs) {
  dc.check.model.params(c("r", "alpha", "a", "b"), params, 
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
  return(sum(custs * cbgcnbd.LL(params, k, x, t.x, T.cal, litt)))
}


#' Calculate the log-likelihood of the CBG/CNBD-k model
#' 
#' @param params CBG/CNBD-k parameters - a vector with 'r', 'alpha', 'a' and 'b'
#'   in that order.
#' @param k specified degree of regularity for Erlang-k distributed 
#'   interpurchase times
#' @param cal.cbs calibration period CBS. It must contain columns for frequency 
#'   'x', for recency 't.x.' and total time observed 'T.cal'. Optionally a 
#'   column 'custs' can be provided, which represents number of customers with a
#'   specific combination of frequency 'x' and 'T.cal'.
#' @return a vector of log-likelihoods
#' @export
cbgcnbd.LL <- function(params, k, x, t.x, T.cal, litt) {
  max.length <- max(length(x), length(t.x), length(T.cal))
  if (max.length%%length(x)) 
    warning("Maximum vector length not a multiple of the length of x")
  if (max.length%%length(t.x)) 
    warning("Maximum vector length not a multiple of the length of t.x")
  if (max.length%%length(T.cal)) 
    warning("Maximum vector length not a multiple of the length of T.cal")
  if (max.length%%length(litt)) 
    warning("Maximum vector length not a multiple of the length of litt")
  dc.check.model.params(c("r", "alpha", "a", "b"), params, 
   "cbgcnbd.LL")
  if (length(k) != 1 | !is.numeric(k) | k < 1 | k != floor(k))
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
  r <- params[1]
  alpha <- params[2]
  a <- params[3]
  b <- params[4]
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
#' @param params CBG/CNBD-k parameters - a vector with 'r', 'alpha', 'a' and 'b'
#'   in that order.
#' @param k specified degree of regularity for Erlang-k distributed
#'   interpurchase times
#' @param x number of repeat transactions in the calibration period T.cal, or a 
#'   vector of calibration period frequencies.
#' @param t.x recency, i.e. length between first and last transaction during 
#'   calibration period.
#' @param T.cal length of calibration period, or a vector of calibration period 
#'   lengths.
#' @return Probability that the customer is still alive at the end of the
#'   calibration period.
#' @export
cbgcnbd.PAlive <- function(params, k, x, t.x, T.cal) {
  max.length <- max(length(x), length(t.x), length(T.cal))
  if (max.length%%length(x)) 
    warning("Maximum vector length not a multiple of the length of x")
  if (max.length%%length(t.x)) 
    warning("Maximum vector length not a multiple of the length of t.x")
  if (max.length%%length(T.cal)) 
    warning("Maximum vector length not a multiple of the length of T.cal")
  dc.check.model.params(c("r", "alpha", "a", "b"), params, 
   "cbgcnbd.PAlive")  
  if (length(k) != 1 | !is.numeric(k) | k < 1 | k != floor(k))
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
  r <- params[1]
  alpha <- params[2]
  a <- params[3]
  b <- params[4]  
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
#' Uses CBG/CNBD-k model parameters and a customer's past transaction behavior to 
#' return the number of transactions they are expected to make in a given time 
#' period.
#' 
#' @param params CBG/CNBD-k parameters - a vector with 'r', 'alpha', 'a' and 'b' in
#'   that order.
#' @param k specified degree of regularity for Erlang-k distributed
#'   interpurchase times
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
cbgcnbd.ConditionalExpectedTransactions <- function(params, k=1, T.star, x, t.x, T.cal) {
  if (k>1)
    message("Results for k>1 are approximative")
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
    "cbgcnbd.ConditionalExpectedTransactions")
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
  # for k>1 we simply scale rate parameter alpha accordingly to approximate results
  alpha <- alpha * k
  P1 <- ((a + b + x) / (a - 1))
  P2 <- 1 - ((alpha+T.cal)/(alpha+T.cal+T.star))^(r+x) * hyperg_2F1(r+x, b+x+1, b+x+a, T.star/(alpha+T.cal+T.star))
  P3 <- cbgcnbd.PAlive(params, k, x, t.x, T.cal)
  return (P1 * P2 * P3)
}

#' CBD/CNBD-k Unconditional Expected Transactions
#' 
#' Uses CBG/CNBD-k model parameters to return the number of transactions a 
#' random customer is expected to make in a given time period.
#' 
#' @param params CBG/CNBD-k parameters - a vector with 'r', 'alpha', 'a' and 'b'
#'   in that order.
#' @param k specified degree of regularity for Erlang-k distributed 
#'   interpurchase times
#' @param t length of time for which we are calculating the expected number of
#'   transactions.
#' @return Number of transactions a customer is expected to make in a time 
#'   period of length t. If any of the input parameters has a length greater
#'   than 1, this will be a vector of expected number of transactions.
#' @export
cbgcnbd.UnconditionalExpectedTransactions <- function(params, k=1, t) {
  if (k>1)
    message("Results for k>1 are approximative")
  dc.check.model.params(c("r", "alpha", "a", "b"), params, 
    "cbgcnbd.UnconditionalExpectedTransactions")
  r <- params[1]
  alpha <- params[2]
  a <- params[3]
  b <- params[4]
  # for k>1 we simply scale rate parameter alpha accordingly to approximate results
  alpha <- alpha * k
  P1 <- ((a + b + x) / (a - 1))
  P2 <- 1 - ((alpha+T.cal)/(alpha+T.cal+T.star))^(r+x) * hyperg_2F1(r+x, b+x+1, b+x+a, T.star/(alpha+T.cal+T.star))
  P3 <- cbgcnbd.PAlive(params, k, x, t.x, T.cal)
  return (P1 * P2 * P3)
}


#' CBD/CNBD-k Unconditional Probability Distribution of Transactions
#' 
#' Uses CBG/CNBD-k model parameters to return the probability distribution of
#' purchase frequencies for a random customer in a given time period, i.e.
#' P(X(t)=x|r,alpha,a,b)
#' 
#' @param params CBG/CNBD-k parameters - a vector with 'r', 'alpha', 'a' and 'b'
#'   in that order.
#' @param k specified degree of regularity for Erlang-k distributed 
#'   interpurchase times
#' @param t length of time for which we are calculating the expected number of 
#'   transactions.
#' @param x number of transactions for which probability is calculated.
#' @return P(X(t)=x|r,alpha,a,b). If any of the input parameters has a length
#'   greater than 1, this will be a vector of expected number of transactions.
#' @export
cbgcnbd.Px <- function(params, k=1, t, x) {
  nbd.Px <- function(params, t, x) {
    r <- params[1]
    alpha <- params[2]
    a <- params[3]
    b <- params[4]
    return ((gamma(r+x)*alpha^r*t^x)/(factorial(x)*gamma(r)*(alpha+t)^(r+x)))
  }
  r <- params[1]
  alpha <- params[2]
  a <- params[3]
  b <- params[4]
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
#' @param params CBG/CNBD-k parameters - a list with 'r', 'alpha', 'a' and 'b'
#' @param return.elog boolean - if TRUE then the event log is returned in 
#'   addition to the CBS summary
#' @return list with elements 'cbs' and 'elog' containing data.frames
#' @export
cbgcnbd.GenerateData <- function(n, k, T.cal, T.star=T.cal, params, return.elog=F) {
  # sample intertransaction timings parameter lambda for each customer
  lambdas <- rgamma(n, shape=params$r, rate=params$alpha)
  
  # sample churn-probability p for each customer
  ps <- rbeta(n, params$a, params$b)
  
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
    cbs[i, "alive"] <- churn>length(ts.cal)
    cbs[i, "T.star"] <- T.star
    cbs[i, "x.star"] <- length(ts.star)
    cbs[i, "p"] <- p
    cbs[i, "lambda"] <- lambda
    cbs[i, "churn"] <- churn
    cbs[i, "x.max"] <- max(times)
  }
  elog$date <- as.Date("2001/01/01") + elog$t
  out <- list(cbs=cbs)
  if (return.elog) {
    out$elog <- transform(elog, date=as.Date("2001/01/01") + t)
  }
  return(out)
}
