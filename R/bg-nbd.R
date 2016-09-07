
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
bgnbd.GenerateData <- function(n, T.cal, T.star, params, return.elog = FALSE) {
  bgcnbd.GenerateData(n = n, T.cal = T.cal, T.star = T.star, params = c(k = 1, params), return.elog = return.elog)
}
