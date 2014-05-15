
#######################################
### Implementation of the NBD model ###
#######################################

nbd.EstimateParameters <- function (cal.cbs, par.start = c(1, 1), max.param.value = 10000) 
{
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


nbd.cbs.LL <- function (params, cal.cbs) 
{
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


nbd.LL <- function (params, x, T.cal) 
{
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


nbd.ConditionalExpectedTransactions <- function (params, T.star, x, T.cal) 
{
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


nbd.GenerateData <- function (n, T.cal, T.star, params, return.elog=F) 
{
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
  elog$date <- as.Date("2001/01/01") + elog$t
  out <- list(cbs=cbs)
  if (return.elog) {
    out$elog <- transform(elog, date=as.Date("2001/01/01") + t)
  }
  return(out)
}
