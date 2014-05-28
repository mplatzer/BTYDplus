
#' Parameter Estimation for the CBG/NBD model
#' 
#' @seealso cbgcnbd.EstimateParameters
cbgnbd.EstimateParameters <- function(cal.cbs, par.start = c(1, 1, 1, 1), max.param.value = 10000) {
  return(cbgcnbd.EstimateParameters(cal.cbs, k=1, par.start, max.param.value))
}


#' Calculates the log-likelihood of the CBG/NBD model
#' 
#' @seealso cbgcnbd.cbs.LL
cbgnbd.cbs.LL <- function(params, cal.cbs) {
  return(cbgcnbd.cbs.LL(params, k=1, cal.cbs))
}


#' Calculates the log-likelihood of the CBG/NBD model
#' 
#' @seealso cbgcnbd.LL
cbgnbd.LL <- function(params, x, t.x, T.cal) {
  return(cbgcnbd.LL(params, k=1, x, t.x, T.cal))
}


#' CBG/NBD P(alive)
#' 
#' @seealso cbgcnbd.PAlive
cbgnbd.PAlive <- function(params, x, t.x, T.cal) {
  return(cbgcnbd.PAlive(params, k=1, x, t.x, T.cal))
}


#' CBG/NBD Conditional Expected Transactions
#' 
#' @seealso cbgcnbd.ConditionalExpectedTransactions
cbgnbd.ConditionalExpectedTransactions <- function(params, T.star, x, t.x, T.cal) {
  return(cbgcnbd.ConditionalExpectedTransactions(params, k=1, T.star, x, t.x, T.cal))
}


#' Simulate data according to CBG/NBD model assumptions
#' 
#' @seealso cbgcnbd.GenerateData
cbgnbd.GenerateData <- function(n, T.cal, T.star, params, return.elog=F) {
  return(cbgcnbd.GenerateData(n, k=1, T.cal, T.star, params, return.elog))
}
