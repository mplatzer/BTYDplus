
#' Parameter Estimation for the CBG/NBD model
#'
#' @seealso cbgcnbd.EstimateParameters
cbgnbd.EstimateParameters <- function(cal.cbs, par.start = c(1, 1, 1, 1), max.param.value = 10000) {
  params <- cbgcnbd.EstimateParameters(cal.cbs, par.start, max.param.value, k=1)
  return(params[2:5])
}


#' Calculates the log-likelihood of the CBG/NBD model
#' 
#' @seealso cbgcnbd.cbs.LL
cbgnbd.cbs.LL <- function(params, cal.cbs) {
  if (length(params)==4) params <- c(k=1, params)
  return(cbgcnbd.cbs.LL(params, cal.cbs))
}


#' Calculates the log-likelihood of the CBG/NBD model
#' 
#' @seealso cbgcnbd.LL
cbgnbd.LL <- function(params, x, t.x, T.cal) {
  if (length(params)==4) params <- c(k=1, params)
  return(cbgcnbd.LL(params, x, t.x, T.cal))
}


#' CBG/NBD P(alive)
#' 
#' @seealso cbgcnbd.PAlive
cbgnbd.PAlive <- function(params, x, t.x, T.cal) {
  if (length(params)==4) params <- c(k=1, params)
  return(cbgcnbd.PAlive(params, x, t.x, T.cal))
}


#' CBG/NBD Conditional Expected Transactions
#' 
#' @seealso cbgcnbd.ConditionalExpectedTransactions
cbgnbd.ConditionalExpectedTransactions <- function(params, T.star, x, t.x, T.cal) {
  if (length(params)==4) params <- c(k=1, params)
  return(cbgcnbd.ConditionalExpectedTransactions(params, T.star, x, t.x, T.cal))
}


#' Simulate data according to CBG/NBD model assumptions
#' 
#' @seealso cbgcnbd.GenerateData
cbgnbd.GenerateData <- function(n, T.cal, T.star, params, return.elog=F) {
  if (length(params)==4) params <- c(k=1, params)
  return(cbgcnbd.GenerateData(n, T.cal, T.star, params, return.elog))
}
