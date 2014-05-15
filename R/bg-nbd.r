
##########################################
### Implementation of the BG/NBD model ###
##########################################

bgnbd.EstimateParameters <- function (cal.cbs, par.start = c(1, 1, 1, 1), max.param.value = 10000) 
{
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
  return(estimated.params)
}


bgnbd.cbs.LL <- function (params, cal.cbs) 
{
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


bgnbd.LL <- function (params, x, t.x, T.cal) 
{
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



bgnbd.PAlive <- function (params, x, t.x, T.cal) 
{
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


bgnbd.ConditionalExpectedTransactions <- function (params, T.star, x, t.x, T.cal) 
{
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


bgnbd.GenerateData <- function (n, T.cal, T.star, params, return.elog=F) 
{
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
  }
  elog$date <- as.Date("2001/01/01") + elog$t
  out <- list(cbs=cbs)
  if (return.elog) {
    out$elog <- transform(elog, date=as.Date("2001/01/01") + t)
  }
  return(out)
}
