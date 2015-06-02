
#' Parameter Estimation for the MBG/NBD model
#'
#' @seealso \code{\link{mbgcnbd.EstimateParameters}}
mbgnbd.EstimateParameters <- function(cal.cbs, par.start = c(1, 1, 1, 1), max.param.value = 10000) {
  params <- mbgcnbd.EstimateParameters(cal.cbs, par.start, max.param.value, k=1)
  return(params[2:5])
}


#' Calculates the log-likelihood of the MBG/NBD model
#' 
#' @seealso \code{\link{mbgcnbd.cbs.LL}}
mbgnbd.cbs.LL <- function(params, cal.cbs) {
  if (length(params)==4) params <- c(k=1, params)
  return(mbgcnbd.cbs.LL(params, cal.cbs))
}


#' Calculates the log-likelihood of the MBG/NBD model
#' 
#' @seealso \code{\link{mbgcnbd.LL}}
mbgnbd.LL <- function(params, x, t.x, T.cal) {
  if (length(params)==4) params <- c(k=1, params)
  return(mbgcnbd.LL(params, x, t.x, T.cal))
}


#' MBG/NBD P(alive)
#' 
#' @seealso \code{\link{mbgcnbd.PAlive}}
mbgnbd.PAlive <- function(params, x, t.x, T.cal) {
  if (length(params)==4) params <- c(k=1, params)
  return(mbgcnbd.PAlive(params, x, t.x, T.cal))
}


#' MBG/NBD Conditional Expected Transactions
#' 
#' @seealso \code{\link{mbgcnbd.ConditionalExpectedTransactions}}
mbgnbd.ConditionalExpectedTransactions <- function(params, T.star, x, t.x, T.cal) {
  if (length(params)==4) params <- c(k=1, params)
  return(mbgcnbd.ConditionalExpectedTransactions(params, T.star, x, t.x, T.cal))
}


#' Simulate data according to MBG/NBD model assumptions
#' 
#' @seealso \code{\link{mbgcnbd.GenerateData}}
mbgnbd.GenerateData <- function(n, T.cal, T.star, params, return.elog=FALSE) {
  if (length(params)==4) params <- c(k=1, params)
  return(mbgcnbd.GenerateData(n, T.cal, T.star, params, return.elog))
}
