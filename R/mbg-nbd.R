
#' Parameter Estimation for the MBG/NBD model
#' 
#' Estimates parameters for the MBG/NBD model via Maximum Likelihood 
#' Estimation.
#' 
#' @param cal.cbs calibration period CBS. It must contain columns for frequency 
#'   \code{x}, for recency \code{t.x} and total time observed \code{T.cal}. Optionally a 
#'   column \code{custs} can be provided, which represents number of customers with a
#'   specific combination of frequency \code{x}, recency \code{t.x} and \code{T.cal}.
#' @param par.start initial MBG/NBD parameters - a vector with \code{r}, \code{alpha}, 
#'   \code{a} and \code{b} in that order.
#' @param max.param.value the upper bound on parameters
#' @return list of estimated parameters
#' @export
#' @seealso \code{\link{bgcnbd.EstimateParameters}}
#' @references Batislam, E.P., M. Denizel, A. Filiztekin. 2007. Empirical
#'   validation and comparison of models for customer base analysis.
#'   International Journal of Research in Marketing 24(3) 201-209. - Hoppe,
#'   Daniel, and Udo Wagner. 'Customer base analysis: The case for a central
#'   variant of the Betageometric/NBD Model.' Marketing Journal of Research and
#'   Management 3.2 (2007): 75-90.
mbgnbd.EstimateParameters <- function(cal.cbs, par.start = c(1, 1, 1, 1), max.param.value = 10000) {
  params <- mbgcnbd.EstimateParameters(cal.cbs, par.start, max.param.value, k = 1)
  return(params[2:5])
}


#' Calculate the log-likelihood of the MBG/NBD model
#' 
#' @param params MBG/NBD parameters - a vector with \code{r}, \code{alpha}, \code{a} 
#'   and \code{b} in that order.
#' @param cal.cbs calibration period CBS. It must contain columns for frequency 
#'   \code{x}, for recency \code{t.x}, for sum of logarithmic interpurchase times \code{litt}
#'   and total time observed \code{T.cal}. Optionally a column \code{custs} can be
#'   provided, which represents number of customers with a specific combination
#'   of frequency \code{x} and \code{T.cal}.
#' @return the total log-likelihood for the provided data.
#' @export
#' @seealso \code{\link{bgcnbd.cbs.LL}}
#' @seealso \code{\link{mbgnbd.EstimateParameters}}
#' @references Batislam, E.P., M. Denizel, A. Filiztekin. 2007. Empirical
#'   validation and comparison of models for customer base analysis.
#'   International Journal of Research in Marketing 24(3) 201-209. - Hoppe,
#'   Daniel, and Udo Wagner. 'Customer base analysis: The case for a central
#'   variant of the Betageometric/NBD Model.' Marketing Journal of Research and
#'   Management 3.2 (2007): 75-90.
mbgnbd.cbs.LL <- function(params, cal.cbs) {
  if (length(params) == 4) 
    params <- c(k = 1, params)
  return(mbgcnbd.cbs.LL(params, cal.cbs))
}


#' Calculate the log-likelihood of the MBG/NBD model
#' 
#' @param params MBG/NBD parameters - a vector with \code{r}, \code{alpha}, \code{a}
#'   and \code{b} in that order.
#' @param x frequency, i.e. number of re-purchases
#' @param t.x recency, i.e. time elapsed from first purchase to last purchase
#' @param T.cal total time of observation period
#' @return a vector of log-likelihoods
#' @export
#' @seealso \code{\link{bgcnbd.LL}}
#' @seealso \code{\link{mbgnbd.EstimateParameters}}
#' @references Batislam, E.P., M. Denizel, A. Filiztekin. 2007. Empirical
#'   validation and comparison of models for customer base analysis.
#'   International Journal of Research in Marketing 24(3) 201-209. - Hoppe,
#'   Daniel, and Udo Wagner. 'Customer base analysis: The case for a central
#'   variant of the Betageometric/NBD Model.' Marketing Journal of Research and
#'   Management 3.2 (2007): 75-90.
mbgnbd.LL <- function(params, x, t.x, T.cal) {
  if (length(params) == 4) 
    params <- c(k = 1, params)
  return(mbgcnbd.LL(params, x, t.x, T.cal))
}


#' MBG/NBD P(alive)
#' 
#' Uses MBG/NBD model parameters and a customer's past transaction behavior 
#' to return the probability that they are still alive at the end of the 
#' calibration period.
#' 
#' @param params MBG/NBD parameters - a vector with \code{r}, \code{alpha}, \code{a}
#'   and \code{b} in that order.
#' @param x number of repeat transactions in the calibration period T.cal, or a 
#'   vector of calibration period frequencies.
#' @param t.x recency, i.e. length between first and last transaction during 
#'   calibration period.
#' @param T.cal length of calibration period, or a vector of calibration period 
#'   lengths.
#' @return Probability that the customer is still alive at the end of the 
#'   calibration period.
#' @export
#' @seealso \code{\link{bgcnbd.PAlive}}
#' @seealso \code{\link{mbgnbd.EstimateParameters}}
#' @references Batislam, E.P., M. Denizel, A. Filiztekin. 2007. Empirical
#'   validation and comparison of models for customer base analysis.
#'   International Journal of Research in Marketing 24(3) 201-209. - Hoppe,
#'   Daniel, and Udo Wagner. 'Customer base analysis: The case for a central
#'   variant of the Betageometric/NBD Model.' Marketing Journal of Research and
#'   Management 3.2 (2007): 75-90.
mbgnbd.PAlive <- function(params, x, t.x, T.cal) {
  if (length(params) == 4) 
    params <- c(k = 1, params)
  return(mbgcnbd.PAlive(params, x, t.x, T.cal))
}


#' MBG/NBD Conditional Expected Transactions
#' 
#' Uses MBG/NBD model parameters and a customer's past transaction behavior 
#' to return the number of transactions they are expected to make in a given 
#' time period.
#' 
#' @param params MBG/NBD parameters - a vector with \code{r}, \code{alpha}, \code{a}
#'   and \code{b} in that order.
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
#' @seealso \code{\link{bgcnbd.ConditionalExpectedTransactions}}
#' @seealso \code{\link{mbgnbd.EstimateParameters}}
#' @references Batislam, E.P., M. Denizel, A. Filiztekin. 2007. Empirical
#'   validation and comparison of models for customer base analysis.
#'   International Journal of Research in Marketing 24(3) 201-209. - Hoppe,
#'   Daniel, and Udo Wagner. 'Customer base analysis: The case for a central
#'   variant of the Betageometric/NBD Model.' Marketing Journal of Research and
#'   Management 3.2 (2007): 75-90.
mbgnbd.ConditionalExpectedTransactions <- function(params, T.star, x, t.x, T.cal) {
  if (length(params) == 4) 
    params <- c(k = 1, params)
  return(mbgcnbd.ConditionalExpectedTransactions(params, T.star, x, t.x, T.cal))
}


#' MBG/NBD Unconditional Probability Distribution of Transactions
#' 
#' Uses MBG/NBD model parameters to return the probability distribution of
#' purchase frequencies for a random customer in a given time period, i.e.
#' P(X(t)=x|r,alpha,a,b)
#' 
#' @param params MBG/NBD parameters - a vector with \code{r}, \code{alpha}, \code{a} and \code{b}
#'   in that order.
#' @param t length of time for which we are calculating the expected number of 
#'   transactions.
#' @param x number of transactions for which probability is calculated.
#' @return P(X(t)=x|r,alpha,a,b). If any of the input parameters has a length
#'   greater than 1, this will be a vector of expected number of transactions.
#' @export
#' @seealso \code{\link{bgcnbd.pmf}}
#' @seealso \code{\link{mbgnbd.EstimateParameters}}
#' @references Batislam, E.P., M. Denizel, A. Filiztekin. 2007. Empirical
#'   validation and comparison of models for customer base analysis.
#'   International Journal of Research in Marketing 24(3) 201-209. - Hoppe,
#'   Daniel, and Udo Wagner. 'Customer base analysis: The case for a central
#'   variant of the Betageometric/NBD Model.' Marketing Journal of Research and
#'   Management 3.2 (2007): 75-90.
mbgnbd.pmf <- function(params, t, x) {
  if (length(params) == 4) 
    params <- c(k = 1, params)
  mbgcnbd.pmf(params = params, t = t, x = x)
}


#' Simulate data according to MBG/NBD model assumptions
#' 
#' @param n number of customers
#' @param T.cal length of calibration period
#' @param T.star length of holdout period
#' @param params BG/CNBD parameters - a vector with \code{r}, \code{alpha}, \code{a} and \code{b}
#'   in that order.
#' @param return.elog boolean - if \code{TRUE} then the event log is returned in 
#'   addition to the CBS summary
#' @return list with elements \code{cbs} and \code{elog} containing data.frames
#' @export
#' @seealso \code{\link{bgcnbd.GenerateData}}
#' @references Batislam, E.P., M. Denizel, A. Filiztekin. 2007. Empirical
#'   validation and comparison of models for customer base analysis.
#'   International Journal of Research in Marketing 24(3) 201-209. - Hoppe,
#'   Daniel, and Udo Wagner. 'Customer base analysis: The case for a central
#'   variant of the Betageometric/NBD Model.' Marketing Journal of Research and
#'   Management 3.2 (2007): 75-90.
mbgnbd.GenerateData <- function(n, T.cal, T.star, params, return.elog = FALSE) {
  if (length(params) == 4) 
    params <- c(k = 1, params)
  return(mbgcnbd.GenerateData(n = n, T.cal = T.cal, T.star = T.star, params = params, return.elog = return.elog))
} 
