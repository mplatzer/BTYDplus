
#' Parameter Estimation for the MBG/CNBD-k model
#' 
#' @param ... Arguments to be passed to \code{\link{bgcnbd.EstimateParameters}}
#' @references Platzer Michael, and Thomas Reutterer (forthcoming)
#' @export
mbgcnbd.EstimateParameters <- function(...) {
  bgcnbd.EstimateParameters(..., dropout_at_zero = TRUE)
}


#' Calculate the log-likelihood of the MBG/CNBD-k model
#' 
#' @param ... Arguments to be passed to \code{\link{bgcnbd.cbs.LL}}
#' @references Platzer Michael, and Thomas Reutterer (forthcoming)
#' @export
mbgcnbd.cbs.LL <- function(...) {
  bgcnbd.cbs.LL(..., dropout_at_zero = TRUE)
}


#' Calculate the log-likelihood of the MBG/CNBD-k model
#' 
#' @param ... Arguments to be passed to \code{\link{bgcnbd.LL}}
#' @references Platzer Michael, and Thomas Reutterer (forthcoming)
#' @export
mbgcnbd.LL <- function(...) {
  bgcnbd.LL(..., dropout_at_zero = TRUE)
}


#' MBG/CNBD-k P(alive)
#' 
#' @param ... Arguments to be passed to \code{\link{bgcnbd.PAlive}}
#' @references Platzer Michael, and Thomas Reutterer (forthcoming)
#' @export
mbgcnbd.PAlive <- function(...) {
  bgcnbd.PAlive(..., dropout_at_zero = TRUE)
}


#' MBG/CNBD-k Conditional Expected Transactions
#' 
#' @param ... Arguments to be passed to \code{\link{bgcnbd.ConditionalExpectedTransactions}}
#' @references Platzer Michael, and Thomas Reutterer (forthcoming)
#' @export
mbgcnbd.ConditionalExpectedTransactions <- function(...) {
  bgcnbd.ConditionalExpectedTransactions(..., dropout_at_zero = TRUE)
}


#' MBG/CNBD-k Unconditional Probability Distribution of Transactions
#' 
#' @param ... Arguments to be passed to \code{\link{bgcnbd.pmf}}
#' @references Platzer Michael, and Thomas Reutterer (forthcoming)
#' @export
mbgcnbd.pmf <- function(...) {
  bgcnbd.pmf(..., dropout_at_zero = TRUE)
}


#' Simulate data according to MBG/CNBD-k model assumptions
#' 
#' @param ... Arguments to be passed to \code{\link{bgcnbd.GenerateData}}
#' @references Platzer Michael, and Thomas Reutterer (forthcoming)
#' @export
mbgcnbd.GenerateData <- function(...) {
  bgcnbd.GenerateData(..., dropout_at_zero = TRUE)
}


#' MBG/CNBD-k Plot Frequency in Calibration Period
#' 
#' @param ... Arguments to be passed to \code{\link{bgcnbd.PlotFrequencyInCalibration}}
#' @references Platzer Michael, and Thomas Reutterer (forthcoming)
#' @export
mbgcnbd.PlotFrequencyInCalibration <- function(...) {
  bgcnbd.PlotFrequencyInCalibration(..., dropout_at_zero = TRUE)
}


#' MBG/CNBD-k Expectation
#' 
#' @param ... Arguments to be passed to \code{\link{bgcnbd.Expectation}}
#' @references Platzer Michael, and Thomas Reutterer (forthcoming)
#' @export
mbgcnbd.Expectation <- function(...) {
  bgcnbd.Expectation(..., dropout_at_zero = TRUE)
}


#' MBG/CNBD-k Expected Cumulative Transactions
#' 
#' @param ... Arguments to be passed to \code{\link{bgnbd.ExpectedCumulativeTransactions}}
#' @references Platzer Michael, and Thomas Reutterer (forthcoming)
#' @export
mbgcnbd.ExpectedCumulativeTransactions <- function(...) {
  bgcnbd.ExpectedCumulativeTransactions(..., dropout_at_zero = TRUE)
}


#' MBG/CNBD-k Tracking Cumulative Transactions Plot
#' 
#' @param ... Arguments to be passed to \code{\link{bgcnbd.PlotTrackingCum}}
#' @references Platzer Michael, and Thomas Reutterer (forthcoming)
#' @export
mbgcnbd.PlotTrackingCum <- function(...) {
  bgcnbd.PlotTrackingCum(..., dropout_at_zero = TRUE)
}


#' MBG/CNBD-k Tracking Incremental Transactions Comparison
#' 
#' @param ... Arguments to be passed to \code{\link{bgcnbd.PlotTrackingInc}}
#' @references Platzer Michael, and Thomas Reutterer (forthcoming)
#' @export
mbgcnbd.PlotTrackingInc <- function(...) {
  bgcnbd.PlotTrackingInc(..., dropout_at_zero = TRUE)
}
