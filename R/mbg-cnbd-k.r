
#' Parameter Estimation for the MBG/CNBD-k model
#' 
#' @param ... Arguments to be passed to bgcnbd.EstimateParameters
#' @seealso \code{\link{bgcnbd.EstimateParameters}}
#' @export
mbgcnbd.EstimateParameters <- function(...) {
  bgcnbd.EstimateParameters(..., dropout_at_zero = TRUE)
}


#' Calculate the log-likelihood of the MBG/CNBD-k model
#' 
#' @param ... Arguments to be passed to bgcnbd.cbs.LL
#' @seealso \code{\link{bgcnbd.cbs.LL}}
#' @export
mbgcnbd.cbs.LL <- function(...) {
  bgcnbd.cbs.LL(..., dropout_at_zero = TRUE)
}


#' Calculate the log-likelihood of the MBG/CNBD-k model
#' 
#' @param ... Arguments to be passed to bgcnbd.LL
#' @seealso \code{\link{bgcnbd.LL}}
#' @export
mbgcnbd.LL <- function(...) {
  bgcnbd.LL(..., dropout_at_zero = TRUE)
}


#' MBG/CNBD-k P(alive)
#' 
#' @param ... Arguments to be passed to bgcnbd.PAlive
#' @seealso \code{\link{bgcnbd.PAlive}}
#' @export
mbgcnbd.PAlive <- function(...) {
  bgcnbd.PAlive(..., dropout_at_zero = TRUE)
}


#' MBG/CNBD-k Conditional Expected Transactions
#' 
#' @param ... Arguments to be passed to bgcnbd.ConditionalExpectedTransactions
#' @seealso \code{\link{bgcnbd.ConditionalExpectedTransactions}}
#' @export
mbgcnbd.ConditionalExpectedTransactions <- function(...) {
  bgcnbd.ConditionalExpectedTransactions(..., dropout_at_zero = TRUE)
}


#' MBG/CNBD-k Unconditional Probability Distribution of Transactions
#' 
#' @param ... Arguments to be passed to bgcnbd.pmf
#' @seealso \code{\link{bgcnbd.pmf}}
#' @export
mbgcnbd.pmf <- function(...) {
  bgcnbd.pmf(..., dropout_at_zero = TRUE)
}


#' Simulate data according to MBG/CNBD-k model assumptions
#' 
#' @param ... Arguments to be passed to bgcnbd.GenerateData
#' @seealso \code{\link{bgcnbd.GenerateData}}
#' @export
mbgcnbd.GenerateData <- function(...) {
  bgcnbd.GenerateData(..., dropout_at_zero = TRUE)
}


#' MBG/CNBD-k Plot Frequency in Calibration Period
#' 
#' @param ... Arguments to be passed to bgcnbd.PlotFrequencyInCalibration
#' @seealso \code{\link{bgcnbd.PlotFrequencyInCalibration}}
#' @export
mbgcnbd.PlotFrequencyInCalibration <- function(...) {
  bgcnbd.PlotFrequencyInCalibration(..., dropout_at_zero = TRUE)
}


#' MBG/CNBD-k Expectation
#' 
#' @param ... Arguments to be passed to bgcnbd.Expectation
#' @seealso \code{\link{bgcnbd.Expectation}}
#' @export
mbgcnbd.Expectation <- function(...) {
  bgcnbd.Expectation(..., dropout_at_zero = TRUE)  
}


#' MBG/CNBD-k Expected Cumulative Transactions
#' 
#' @param ... Arguments to be passed to bgnbd.ExpectedCumulativeTransactions
#' @seealso \code{\link{bgcnbd.ExpectedCumulativeTransactions}}
#' @export
mbgcnbd.ExpectedCumulativeTransactions <- function(...) {
  bgcnbd.ExpectedCumulativeTransactions(..., dropout_at_zero = TRUE)
}


#' MBG/CNBD-k Tracking Cumulative Transactions Plot
#' 
#' @param ... Arguments to be passed to bgcnbd.PlotTrackingCum
#' @seealso \code{\link{bgcnbd.PlotTrackingCum}}
#' @export
mbgcnbd.PlotTrackingCum <- function(...) {
  bgcnbd.PlotTrackingCum(..., dropout_at_zero = TRUE)
}


#' MBG/CNBD-k Tracking Incremental Transactions Comparison
#' 
#' @param ... Arguments to be passed to bgcnbd.PlotTrackingInc
#' @seealso \code{\link{bgcnbd.PlotTrackingInc}}
#' @export
mbgcnbd.PlotTrackingInc <- function(...) {
  bgcnbd.PlotTrackingInc(..., dropout_at_zero = TRUE)
}
