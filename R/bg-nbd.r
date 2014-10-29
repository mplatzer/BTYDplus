
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
bgnbd.GenerateData <- function(n, T.cal, T.star, params, return.elog=FALSE) {
  # check model parameters
  dc.check.model.params(c("r", "alpha", "a", "b"), params,
                        "bgnbd.GenerateData")
  
  r <- params[1]
  alpha <- params[2]
  a <- params[3]
  b <- params[4]
  
  if (length(T.cal)==1) T.cal <- rep(T.cal, n)  
  if (length(T.star)==1) T.star <- rep(T.star, n)  
  
  # sample intertransaction timings parameter lambda for each customer
  lambdas <- rgamma(n, shape=r, rate=alpha)

  # sample dropout probability p for each customer
  ps <- rbeta(n, a, b)
  
  # sample intertransaction timings & churn
  cbs_list <- list()
  elog_list <- list()
  for (i in 1:n) {
    p <- ps[i]
    lambda <- lambdas[i]
    # sample no. of transactions until churn
    churn <- which.max(rbinom(10/p, 1, p))
    # sample transaction times
    times <- cumsum(c(0, rexp(churn, rate=lambda)))
    if (return.elog)
      elog_list[[i]] <- data.frame(cust=i, t=times[times<(T.cal[i]+T.star[i])])    
    # determine frequency, recency, etc.
    ts.cal <- times[times<T.cal[i]]
    ts.star <- times[times>=T.cal[i] & times<(T.cal[i]+T.star[i])]
    cbs_list[[i]] <- list(cust   = i,
                          x      = length(ts.cal)-1,
                          t.x    = max(ts.cal),
                          churn  = churn,
                          alive  = churn > (length(ts.cal)-1),
                          x.star = length(ts.star))
  }
  cbs <- do.call(rbind.data.frame, cbs_list)
  cbs$lambda <- lambdas
  cbs$p      <- ps
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
