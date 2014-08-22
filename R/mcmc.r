
#' Calculates P(alive) based on MCMC draws
#'
#' @param cal.cbs data.frame with column \code{T.cal}
#' @param draws MCMC draws returned by \code{\link{pnbd.mcmc.DrawParameters}}, \code{\link{pcnbd.mcmc.DrawParameters}}, \code{\link{abe.mcmc.DrawParameters}}
#' @return numeric vector of probabilities
#' @export
mcmc.PAlive <- function(cal.cbs, draws) {
  
  nr_of_cust <- length(draws$level_1)
  if (nr_of_cust != nrow(cal.cbs))
    stop("mismatch between number of customers in parameters 'cal.cbs' and 'draws'")
  
  p.alives <- sapply(1:nr_of_cust, function(i) mean(as.matrix(draws$level_1[[i]][, "z"])))
  return(p.alives)
}


#' Samples number of future transactions based on drawn parameters
#'
#' @param cal.cbs data.frame with column \code{t.x} and \code{T.cal}
#' @param draws MCMC draws returned by \code{\link{pnbd.mcmc.DrawParameters}}, \code{\link{pcnbd.mcmc.DrawParameters}}, \code{\link{abe.mcmc.DrawParameters}}
#' @param T.star length of period for which future transactions are counted
#' @return 2-dim array [draw x cust] with sampled future transactions
#' @export
mcmc.DrawFutureTransactions <- function(cal.cbs, draws, T.star=cal.cbs$T.star) {
  
  nr_of_draws <- niter(draws$level_2) * nchain(draws$level_2)
  nr_of_cust <- length(draws$level_1)
  parameters <- varnames(draws$level_1[[1]])
  
  if (nr_of_cust != nrow(cal.cbs))
    stop("mismatch between number of customers in parameters 'cal.cbs' and 'draws'")
  if (is.null(T.star))
    stop("T.star is missing")
  
  x.stars <- array(NA_real_, dim=c(nr_of_draws, nr_of_cust))
  if (length(T.star)==1) T.star <- rep(T.star, nr_of_cust)
  
  draw_left_truncated_gamma <- function(lower, k, lambda) {
    rand <- runif(1, pgamma(lower, k, k*lambda), 1)
    qgamma(rand, k, k*lambda)
  }
  
  for (cust in 1:nrow(cal.cbs)) {
    Tcal    <- cal.cbs[cust, "T.cal"]
    Tstar   <- T.star[cust]
    tx      <- cal.cbs[cust, "t.x"]
    taus    <- as.matrix(draws$level_1[[cust]][, "tau"])
    ks      <- if ("k" %in% parameters) as.matrix(draws$level_1[[cust]][, "k"]) else rep(1, nr_of_draws)
    lambdas <- as.matrix(draws$level_1[[cust]][, "lambda"])
    alive   <- (taus>Tcal)
    
    # Case: customer alive
    for (draw in which(alive)) {
      # sample itt which is larger than (Tcal-tx)
      itts <- draw_left_truncated_gamma(Tcal-tx, ks[draw], lambdas[draw])
      # sample 'sufficiently' large amount of inter-transaction times
      minT <- pmin(Tcal + Tstar - tx, taus[draw] - tx)
      nr_of_itt_draws <- pmax(10, round(minT * lambdas[draw]))
      itts <- c(itts, rgamma(nr_of_itt_draws * 2, shape=ks[draw], rate=ks[draw]*lambdas[draw]))
      if (sum(itts)<minT) itts <- c(itts, rgamma(nr_of_itt_draws * 4, shape=ks[draw], rate=ks[draw]*lambdas[draw]))
      if (sum(itts)<minT) itts <- c(itts, rgamma(nr_of_itt_draws * 800, shape=ks[draw], rate=ks[draw]*lambdas[draw]))
      if (sum(itts)<minT)
        stop("not enough inter-transaction times sampled! cust:", cust, " draw:", draw, " ", sum(itts), " < ", minT)
      x.stars[draw, cust] <- sum(cumsum(itts)<minT)
    }
    
    # Case: customer churned
    if (any(!alive)) {
      x.stars[!alive, cust] <- 0
    }
  }
  return(x.stars)
}


#' (Re-)set burnin of MCMC chains.
#'
#' @param draws MCMC draws returned by \code{\link{pnbd.mcmc.DrawParameters}}, \code{\link{pcnbd.mcmc.DrawParameters}}, \code{\link{abe.mcmc.DrawParameters}}
#' @param burnin new start index
#' @return 2-element list with MCMC draws
#' @export
mcmc.setBurnin <- function(draws, burnin) {
  if (burnin < start(draws$level_2) | burnin > end(draws$level_2))
    stop("specified burnin is out of bound: ", start(draws$level_2), " - ", end(draws$level_2))
  draws$level_2 <- window(draws$level_2, start=burnin)
  draws$level_1 <- lapply(draws$level_1, function(draw) window(draw, start=burnin))
  return(draws)
}
