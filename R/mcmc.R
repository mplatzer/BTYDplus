
#' Calculates P(alive) based on MCMC parameter draws
#'
#' @param draws MCMC draws as returned by \code{*.mcmc.DrawParameters}
#' @return Numeric vector with the customers' probabilities of being still alive
#'   at end of calibration period
#' @export
#' @examples
#' data("groceryElog")
#' cbs <- elog2cbs(groceryElog)
#' param.draws <- pnbd.mcmc.DrawParameters(cbs,
#'   mcmc = 100, burnin = 50, thin = 10, chains = 1) # short MCMC to run demo fast
#' palive <- mcmc.PAlive(param.draws)
#' head(palive)
#' mean(palive)
mcmc.PAlive <- function(draws) {
  nr_of_cust <- length(draws$level_1)
  p.alives <- sapply(1:nr_of_cust, function(i) mean(as.matrix(draws$level_1[[i]][, "z"])))
  return(p.alives)
}


#' Draws number of future transactions based on MCMC parameter draws
#'
#' For each customer and each provided MCMC parameter draw this method will
#' sample the number of transactions during the holdout period \code{T.star}. If
#' argument \code{size} is provided then it returns a flexible number of draws,
#' whereas for each customer and each draw it will first make a draw from the
#' parameter draws.
#'
#' @param cal.cbs Calibration period customer-by-sufficient-statistic (CBS)
#'   data.frame.
#' @param draws MCMC draws as returned by \code{*.mcmc.DrawParameters}
#' @param T.star Length of period for which future transactions are counted.
#' @param sample_size Number of samples to draw. Defaults to the same number of
#'   parameter draws that are passed to \code{draws}.
#' @return 2-dim matrix [draw x customer] with sampled future transactions.
#' @export
#' @examples
#' data("groceryElog")
#' cbs <- elog2cbs(groceryElog, T.cal = "2006-12-31")
#' param.draws <- pnbd.mcmc.DrawParameters(cbs,
#'   mcmc = 100, burnin = 50, thin = 10, chains = 1) # short MCMC to run demo fast
#' xstar.draws <- mcmc.DrawFutureTransactions(cbs, param.draws)
#' cbs$xstar.est <- apply(xstar.draws, 2, mean)
#' cbs$pactive <- mcmc.PActive(xstar.draws)
#' head(cbs)
mcmc.DrawFutureTransactions <- function(cal.cbs, draws, T.star = cal.cbs$T.star, sample_size = NULL) {

  if (is.null(sample_size)) {
    nr_of_draws <- niter(draws$level_2) * nchain(draws$level_2)
  } else {
    stopifnot(is.numeric(sample_size))
    nr_of_draws <- as.integer(sample_size)
  }
  stopifnot(nr_of_draws >= 1)
  nr_of_cust <- length(draws$level_1)
  parameters <- varnames(draws$level_1[[1]])

  if (nr_of_cust != nrow(cal.cbs))
    stop("mismatch between number of customers in parameters 'cal.cbs' and 'draws'")
  if (is.null(T.star))
    stop("T.star is missing")

  x.stars <- array(NA_real_, dim = c(nr_of_draws, nr_of_cust))
  if (length(T.star) == 1)
    T.star <- rep(T.star, nr_of_cust)

  draw_left_truncated_gamma <- function(lower, k, lambda) {
    rand <- runif(1, pgamma(lower, k, k * lambda), 1)
    qgamma(rand, k, k * lambda)
  }

  for (cust in 1:nrow(cal.cbs)) {
    Tcal <- cal.cbs$T.cal[cust]
    Tstar <- T.star[cust]
    tx <- cal.cbs$t.x[cust]
    taus <- drop(as.matrix(draws$level_1[[cust]][, "tau"]))
    if ("k" %in% parameters) {
      ks <- drop(as.matrix(draws$level_1[[cust]][, "k"]))
    } else {
      ks <- rep(1, length(taus))
    }
    lambdas <- drop(as.matrix(draws$level_1[[cust]][, "lambda"]))
    stopifnot(length(taus) == length(ks) && length(taus) == length(lambdas))
    if (!is.null(sample_size)) {
      idx <- sample(length(taus), size = sample_size, replace = TRUE)
      taus <- taus[idx]
      ks <- ks[idx]
      lambdas <- lambdas[idx]
    }
    alive <- (taus > Tcal)

    # Case: customer alive
    for (draw in which(alive)) {
      # sample itt which is larger than (Tcal-tx)
      itts <- draw_left_truncated_gamma(Tcal - tx, ks[draw], lambdas[draw])
      # sample 'sufficiently' large amount of inter-transaction times
      minT <- pmin(Tcal + Tstar - tx, taus[draw] - tx)
      nr_of_itt_draws <- pmax(10, round(minT * lambdas[draw]))
      itts <- c(itts, rgamma(nr_of_itt_draws * 2, shape = ks[draw], rate = ks[draw] * lambdas[draw]))
      if (sum(itts) < minT)
        itts <- c(itts, rgamma(nr_of_itt_draws * 4, shape = ks[draw], rate = ks[draw] * lambdas[draw]))
      if (sum(itts) < minT)
        itts <- c(itts, rgamma(nr_of_itt_draws * 800, shape = ks[draw], rate = ks[draw] * lambdas[draw]))
      if (sum(itts) < minT)
        stop("not enough inter-transaction times sampled! cust:", cust, " draw:", draw, " ", sum(itts),
          " < ", minT)
      x.stars[draw, cust] <- sum(cumsum(itts) < minT)
    }

    # Case: customer churned
    if (any(!alive)) {
      x.stars[!alive, cust] <- 0
    }
  }
  return(x.stars)
}


#' Calculates P(active) based on drawn future transactions.
#'
#' @param xstar Future transaction draws as returned by
#'   \code{\link{mcmc.DrawFutureTransactions}}.
#' @return numeric A vector with the customers' probabilities of being active
#'   during the holdout period.
#' @export
#' @examples
#' data("groceryElog")
#' cbs <- elog2cbs(groceryElog, T.cal = "2006-12-31")
#' param.draws <- pnbd.mcmc.DrawParameters(cbs,
#'   mcmc = 100, burnin = 50, thin = 10, chains = 1) # short MCMC to run demo fast
#' xstar.draws <- mcmc.DrawFutureTransactions(cbs, param.draws)
#' cbs$pactive <- mcmc.PActive(xstar.draws)
#' head(cbs)
mcmc.PActive <- function(xstar) {
  return(apply(xstar, 2, function(x) mean(x > 0)))
}


#' (Re-)set burnin of MCMC chains.
#'
#' @param draws MCMC draws as returned by \code{*.mcmc.DrawParameters}
#' @param burnin New start index.
#' @return 2-element list with MCMC draws
#' @export
#' @examples
#' data("groceryElog")
#' cbs <- elog2cbs(groceryElog)
#' param.draws <- pnbd.mcmc.DrawParameters(cbs,
#'   mcmc = 100, burnin = 50, thin = 10, chains = 1) # short MCMC to run demo fast
#' param.draws.stable <- mcmc.setBurnin(param.draws, burnin = 80)
mcmc.setBurnin <- function(draws, burnin) {
  if (burnin < start(draws$level_2) | burnin > end(draws$level_2))
    stop("specified burnin is out of bound: ", start(draws$level_2), " - ", end(draws$level_2))
  draws$level_2 <- window(draws$level_2, start = burnin)
  draws$level_1 <- lapply(draws$level_1, function(draw) window(draw, start = burnin))
  return(draws)
}


#' Draw diagnostic plot to inspect error in P(active).
#'
#' @param cbs A data.frame with column \code{x} and \code{x.star}.
#' @param xstar Future transaction draws as returned by
#'   \code{\link{mcmc.DrawFutureTransactions}}.
#' @param title Plot title.
#' @export
#' @examples
#' data("groceryElog")
#' cbs <- elog2cbs(groceryElog, T.cal = "2006-12-31")
#' param.draws <- pnbd.mcmc.DrawParameters(cbs,
#'   mcmc = 100, burnin = 50, thin = 10, chains = 1) # short MCMC to run demo fast
#' xstar.draws <- mcmc.DrawFutureTransactions(cbs, param.draws)
#' mcmc.plotPActiveDiagnostic(cbs, xstar.draws)
mcmc.plotPActiveDiagnostic <- function(cbs, xstar, title = "Diagnostic Plot for P(active)") {
  pactive <- mcmc.PActive(xstar)
  x.star <- cbs$x.star
  op <- par(mar = c(4, 4, 2, 2), mgp = c(2.5, 1, 0))
  cuts <- unique(quantile(c(0, pactive, 1), seq(0, 1, 0.1)))
  spls.y <- sapply(split(x.star > 0, cut(pactive, breaks = cuts, include.lowest = T)), mean)
  spls.x <- sapply(split(pactive, cut(pactive, breaks = cuts, include.lowest = T)), mean)
  plot(spls.x, spls.y, typ = "b", xlim = c(0, 1), ylim = c(0, 1), frame = 0, axes = F,
       xlab = "Estimated P(active)", ylab = "Actual Share of Actives", main = title)
  axis(side = 1, at = seq(0, 1, 0.1), pos = 0, labels = paste(100 * seq(0, 1, 0.1), "%"))
  axis(side = 2, at = seq(0, 1, 0.1), pos = 0, labels = paste(100 * seq(0, 1, 0.1), "%"), las = 2)
  abline(0, 1)
  abline(h = seq(0, 1, 0.1), col = "lightgray", lty = "dotted")
  abline(v = seq(0, 1, 0.1), col = "lightgray", lty = "dotted")
  points(mean(pactive[cbs$x == 0]), mean(x.star[cbs$x == 0] > 0), col = "red", pch = "0")
  par(op)
  invisible()
}


#' Probability Mass Function for Pareto/GGG, Pareto/NBD (HB) and Pareto/NBD (Abe)
#'
#' Return the probability distribution of purchase frequencies for a random
#' customer in a given time period, i.e. \eqn{P(X(t)=x)}. This is estimated by
#' generating \code{sample_size} number of random customers that follow the
#' provided parameter draws. Due to this sampling, the return result varies from
#' one call to another.
#'
#' @param draws MCMC draws as returned by \code{*.mcmc.DrawParameters}
#' @param t Length of time for which we are calculating the expected number of
#'   transactions. May also be a vector.
#' @param x Number of transactions for which probability is calculated. May also
#'   be a vector.
#' @param sample_size Sample size for estimating the probability distribution.
#' @param covariates (optional) Matrix of covariates, for Pareto/NBD (Abe)
#'   model, passed to \code{\link{abe.GenerateData}} for simulating data.
#' @return \eqn{P(X(t)=x)}. If either \code{t} or \code{x} is a vector, then the
#'   output will be a vector as well. If both are vectors, the output will be a
#'   matrix.
#' @export
#' @examples
#' data("groceryElog")
#' cbs <- elog2cbs(groceryElog)
#' param.draws <- pnbd.mcmc.DrawParameters(cbs,
#'   mcmc = 100, burnin = 50, thin = 10, chains = 1) # short MCMC to run demo fast
#' mcmc.pmf(param.draws, t = c(26, 52), x = 0:6)
mcmc.pmf <- function(draws, t, x, sample_size = 10000, covariates = NULL) {
  cohort_draws <- as.matrix(draws$level_2)
  nr_of_draws <- nrow(cohort_draws)
  # use posterior mean
  draw_idx_cnt <- table(sample(nr_of_draws, size = sample_size, replace = TRUE))
  model <- ifelse(all(c("r", "alpha") %in% colnames(cohort_draws)), "pggg", "abe")
  cbs <- rbindlist(lapply(names(draw_idx_cnt), function(idx) {
    n <- unname(draw_idx_cnt[idx])
    if (model == "pggg") {
      params <- as.list(cohort_draws[as.integer(idx), ])
      pggg.GenerateData(n = n, T.cal = 0, T.star = unique(t), params = params)$cbs
    } else if (model == "abe") {
      p <- cohort_draws[as.integer(idx), ]
      params <- list()
      params$beta  <- matrix(p[grepl("^log\\_", names(p))], byrow = TRUE, ncol = 2)
      params$gamma <- matrix(c(p["var_log_lambda"],
                               p["cov_log_lambda_log_mu"],
                               p["cov_log_lambda_log_mu"],
                               p["var_log_mu"]),
                             ncol = 2)
      abe.GenerateData(n = n, T.cal = 0, T.star = unique(t), params = params,
                       covariates = covariates)$cbs
    }
  }))
  pmf <- sapply(1:length(t), function(idx) {
    col <- ifelse(uniqueN(t) == 1, "x.star", paste0("x.star", t[idx]))
    stopifnot(col %in% names(cbs))
    sapply(x, function(x) sum(cbs[[col]] == x) / sample_size)
  })
  if (length(x) == 1) pmf <- t(pmf)
  rownames(pmf) <- x
  colnames(pmf) <- t
  drop(pmf)
}


#' Unconditional Expectation for Pareto/GGG, Pareto/NBD (HB) and Pareto/NBD (Abe)
#'
#' Uses model parameter draws to return the expected number of repeat
#' transactions that a randomly chosen customer (for whom we have no prior
#' information) is expected to make in a given time period. \deqn{E(X(t))}.
#'
#' The expected transactions need to be sampled. Due to this sampling, the
#' return result varies from one call to another. Larger values of
#' \code{sample_size} will generate more stable results.
#'
#' @param draws MCMC draws as returned by \code{*.mcmc.DrawParameters}
#' @param t Length of time for which we are calculating the expected number of
#'   transactions. May also be a vector.
#' @param sample_size Sample size for estimating the probability distribution.
#' @return Number of repeat transactions a customer is expected to make in a
#'   time period of length t.
#' @export
#' @examples
#' data("groceryElog")
#' cbs <- elog2cbs(groceryElog)
#' param.draws <- pnbd.mcmc.DrawParameters(cbs,
#'   mcmc = 100, burnin = 50, thin = 10, chains = 1) # short MCMC to run demo fast
#' mcmc.Expectation(param.draws, t = c(26, 52))
mcmc.Expectation <- function(draws, t, sample_size = 10000) {
  pmf <- mcmc.pmf(draws, t, 0:100, sample_size = sample_size)
  apply(pmf * matrix(rep(0:100, length(t)), ncol = length(t)), 2, sum)
}


#' Expected Cumulative Transactions for Pareto/GGG, Pareto/NBD (HB) and
#' Pareto/NBD (Abe)
#'
#' Uses model parameter draws to return the expected number of repeat
#' transactions that a randomly chosen customer (for whom we have no prior
#' information) is expected to make in a given time period.
#'
#' The expected transactions need to be sampled. Due to this sampling, the
#' return result varies from one call to another. Larger values of
#' \code{sample_size} will generate more stable results.
#'
#' @param draws MCMC draws as returned by \code{*.mcmc.DrawParameters}
#' @param T.cal A vector to represent customers' calibration period lengths (in
#'   other words, the \code{T.cal} column from a
#'   customer-by-sufficient-statistic matrix). Considering rounding in order to
#'   speed up calculations.
#' @param T.tot End of holdout period. Must be a single value, not a vector.
#' @param n.periods.final Number of time periods in the calibration and holdout
#'   periods.
#' @param sample_size Sample size for estimating the probability distribution.
#' @param covariates (optional) Matrix of covariates, for Pareto/NBD (Abe)
#'   model, passed to \code{\link{abe.GenerateData}} for simulating data.
#' @return Numeric vector of expected cumulative total repeat transactions by
#'   all customers.
#' @export
#' @examples
#' data("groceryElog")
#' cbs <- elog2cbs(groceryElog)
#' param.draws <- pnbd.mcmc.DrawParameters(cbs,
#'   mcmc = 100, burnin = 50, thin = 10, chains = 1) # short MCMC to run demo fast
#' # Returns a vector containing expected cumulative repeat transactions for 104
#' # weeks, with every eigth week being reported.
#' mcmc.ExpectedCumulativeTransactions(param.draws,
#'   T.cal = cbs$T.cal, T.tot = 104, n.periods.final = 104/8, sample_size = 1000)
mcmc.ExpectedCumulativeTransactions <- function(draws, T.cal, T.tot, n.periods.final,
                                                sample_size = 10000, covariates = NULL) {
  if (any(T.cal < 0) || !is.numeric(T.cal))
    stop("T.cal must be numeric and may not contain negative numbers.")
  if (length(T.tot) > 1 || T.tot < 0 || !is.numeric(T.tot))
    stop("T.cal must be a single numeric value and may not be negative.")
  if (length(n.periods.final) > 1 || n.periods.final < 0 || !is.numeric(n.periods.final))
    stop("n.periods.final must be a single numeric value and may not be negative.")

  cohort_draws <- as.matrix(draws$level_2)
  nr_of_draws <- nrow(cohort_draws)
  model <- ifelse(all(c("r", "alpha") %in% colnames(cohort_draws)), "pggg", "abe")
  elog <- rbindlist(lapply(1:nr_of_draws, function(i) {
    n <- ceiling(sample_size / nr_of_draws)
    if (model == "pggg") {
      params <- as.list(cohort_draws[i, ])
      elog <- pggg.GenerateData(n = n, T.cal = T.tot, T.star = 0, params = params)$elog
    } else if (model == "abe") {
      p <- as.list(cohort_draws[i, ])
      params <- list()
      params$beta  <- matrix(as.numeric(p[grepl("^log\\_", names(p))]), byrow = TRUE, ncol = 2)
      params$gamma <- matrix(as.numeric(c(p["var_log_lambda"],
                               p["cov_log_lambda_log_mu"],
                               p["cov_log_lambda_log_mu"],
                               p["var_log_mu"])),
                             ncol = 2)
      elog <- abe.GenerateData(n = n, T.cal = T.tot, T.star = 0, params = params,
                               covariates = covariates)$elog
    }
    setDT(elog)
    elog$cust <- paste0(elog$cust, "_", i)
    elog <- elog[t > 0] # drop initial transaction
    elog
  }))
  setkey(elog, t)

  intervals <- seq(T.tot / n.periods.final, T.tot, length.out = n.periods.final)
  cust.birth.periods <- max(T.cal) - T.cal
  expected.transactions <- sapply(intervals, function(interval) {
    if (interval <= min(cust.birth.periods))
      return(0)
    t <- interval - cust.birth.periods[cust.birth.periods < interval]
    uEs <- sapply(unique(t), function(ut) elog[t < ut, .N] / sample_size)
    names(uEs) <- unique(t)
    sum(uEs[as.character(t)])
  })
  return(expected.transactions)
}


#' Tracking Cumulative Transactions Plot for Pareto/GGG, Pareto/NBD (HB) and
#' Pareto/NBD (Abe)
#'
#' Plots the actual and expected cumulative total repeat transactions by all
#' customers for the calibration and holdout periods, and returns this
#' comparison in a matrix.
#'
#' The expected transactions need to be sampled. Due to this sampling, the
#' return result varies from one call to another. Larger values of
#' \code{sample_size} will generate more stable results.
#'
#' @param draws MCMC draws as returned by \code{*.mcmc.DrawParameters}
#' @param T.cal A vector to represent customers' calibration period lengths (in
#'   other words, the \code{T.cal} column from a
#'   customer-by-sufficient-statistic matrix). Considering rounding in order to
#'   speed up calculations.
#' @param T.tot End of holdout period. Must be a single value, not a vector.
#' @param actual.cu.tracking.data A vector containing the cumulative number of
#'   repeat transactions made by customers for each period in the total time
#'   period (both calibration and holdout periods).
#' @param xlab Descriptive label for the x axis.
#' @param ylab Descriptive label for the y axis.
#' @param xticklab A vector containing a label for each tick mark on the x axis.
#' @param title Title placed on the top-center of the plot.
#' @param ymax Upper boundary for y axis.
#' @param sample_size Sample size for estimating the probability distribution.
#'   See \code{\link{mcmc.ExpectedCumulativeTransactions}}.
#' @param covariates (optional) Matrix of covariates, for Pareto/NBD (Abe)
#'   model, passed to \code{\link{abe.GenerateData}} for simulating data.
#' @param legend plot legend, defaults to `Actual` and `Model`.
#' @return Matrix containing actual and expected cumulative repeat transactions.
#' @export
#' @seealso \code{\link{mcmc.PlotTrackingInc}}
#'   \code{\link{mcmc.ExpectedCumulativeTransactions}} \code{\link{elog2cum}}
#' @examples
#' \dontrun{
#' data("groceryElog")
#' cbs <- elog2cbs(groceryElog, T.cal = "2006-12-31")
#' cum <- elog2cum(groceryElog)
#' param.draws <- pnbd.mcmc.DrawParameters(cbs)
#' mat <- mcmc.PlotTrackingCum(param.draws,
#'   T.cal = cbs$T.cal,
#'   T.tot = max(cbs$T.cal + cbs$T.star),
#'   actual.cu.tracking.data = cum)
#' }
mcmc.PlotTrackingCum <- function(draws, T.cal, T.tot, actual.cu.tracking.data,
                                 xlab = "Week", ylab = "Cumulative Transactions",
                                 xticklab = NULL, title = "Tracking Cumulative Transactions",
                                 ymax = NULL, sample_size = 10000, covariates = NULL,
                                 legend = c("Actual", "Model")) {

  actual <- actual.cu.tracking.data
  expected <- mcmc.ExpectedCumulativeTransactions(draws, T.cal, T.tot, length(actual),
                                                  sample_size = sample_size, covariates = covariates)

  dc.PlotTracking(actual = actual, expected = expected, T.cal = T.cal,
                  xlab = xlab, ylab = ylab, title = title,
                  xticklab = xticklab, ymax = ymax,
                  legend = legend)
}


#' Tracking Incremental Transactions Plot for Pareto/GGG, Pareto/NBD (HB) and
#' Pareto/NBD (Abe)
#'
#' Plots the actual and expected incremental total repeat transactions by all
#' customers for the calibration and holdout periods, and returns this
#' comparison in a matrix.
#'
#' The expected transactions need to be sampled. Due to this sampling, the
#' return result varies from one call to another. Larger values of
#' \code{sample_size} will generate more stable results.
#'
#' @param draws MCMC draws as returned by \code{*.mcmc.DrawParameters}
#' @param T.cal A vector to represent customers' calibration period lengths (in
#'   other words, the \code{T.cal} column from a
#'   customer-by-sufficient-statistic matrix). Considering rounding in order to
#'   speed up calculations.
#' @param T.tot End of holdout period. Must be a single value, not a vector.
#' @param actual.inc.tracking.data A vector containing the incremental number of
#'   repeat transactions made by customers for each period in the total time
#'   period (both calibration and holdout periods).
#' @param xlab Descriptive label for the x axis.
#' @param ylab Descriptive label for the y axis.
#' @param xticklab A vector containing a label for each tick mark on the x axis.
#' @param title Title placed on the top-center of the plot.
#' @param ymax Upper boundary for y axis.
#' @param sample_size Sample size for estimating the probability distribution.
#'   See \code{\link{mcmc.ExpectedCumulativeTransactions}}.
#' @param covariates (optional) Matrix of covariates, for Pareto/NBD (Abe)
#'   model, passed to \code{\link{abe.GenerateData}} for simulating data.
#' @param legend plot legend, defaults to `Actual` and `Model`.
#' @return Matrix containing actual and expected incremental repeat
#'   transactions.
#' @export
#' @seealso \code{\link{mcmc.PlotTrackingCum}}
#'   \code{\link{mcmc.ExpectedCumulativeTransactions}} \code{\link{elog2inc}}
#' @examples
#' \dontrun{
#' data("groceryElog")
#' cbs <- elog2cbs(groceryElog, T.cal = "2006-12-31")
#' inc <- elog2inc(groceryElog)
#' param.draws <- pnbd.mcmc.DrawParameters(cbs)
#' mat <- mcmc.PlotTrackingInc(param.draws,
#'   T.cal = cbs$T.cal,
#'   T.tot = max(cbs$T.cal + cbs$T.star),
#'   actual.inc.tracking.data = inc)
#' }
mcmc.PlotTrackingInc <- function(draws, T.cal, T.tot, actual.inc.tracking.data,
                                 xlab = "Week", ylab = "Transactions",
                                 xticklab = NULL, title = "Tracking Weekly Transactions",
                                 ymax = NULL, sample_size = 10000, covariates = NULL,
                                 legend = c("Actual", "Model")) {

  actual <- actual.inc.tracking.data
  expected_cum <- mcmc.ExpectedCumulativeTransactions(draws, T.cal, T.tot, length(actual),
                                                      sample_size = sample_size, covariates = covariates)
  expected <- BTYD::dc.CumulativeToIncremental(expected_cum)

  dc.PlotTracking(actual = actual, expected = expected, T.cal = T.cal,
                  xlab = xlab, ylab = ylab, title = title,
                  xticklab = xticklab, ymax = ymax,
                  legend = legend)
}


#' Frequency in Calibration Period for Pareto/GGG, Pareto/NBD (HB) and Pareto/NBD (Abe)
#'
#' Plots a histogram and returns a matrix comparing the actual and expected
#' number of customers who made a certain number of repeat transactions in the
#' calibration period, binned according to calibration period frequencies.
#'
#' The method \code{\link{mcmc.pmf}} is called to calculate the expected numbers
#' based on the corresponding model.
#'
#' @param draws MCMC draws as returned by \code{*.mcmc.DrawParameters}
#' @param cal.cbs Calibration period customer-by-sufficient-statistic (CBS)
#'   data.frame. It must contain columns for frequency ('x') and total time
#'   observed ('T.cal').
#' @param censor Cutoff point for number of transactions in plot.
#' @param xlab Descriptive label for the x axis.
#' @param ylab Descriptive label for the y axis.
#' @param title Title placed on the top-center of the plot.
#' @param sample_size Sample size for estimating the probability distribution.
#'   See \code{\link{mcmc.pmf}}.
#' @return Calibration period repeat transaction frequency comparison matrix
#'   (actual vs. expected).
#' @export
#' @seealso \code{\link{mcmc.pmf}}
#' @examples
#' \dontrun{
#' data("groceryElog")
#' cbs <- elog2cbs(groceryElog, T.cal = "2006-12-31")
#' param.draws <- pnbd.mcmc.DrawParameters(cbs,
#'   mcmc = 100, burnin = 50, thin = 10, chains = 1) # short MCMC to run demo fast
#' mcmc.PlotFrequencyInCalibration(param.draws, cbs, sample_size = 100)
#' }
mcmc.PlotFrequencyInCalibration <- function(draws, cal.cbs, censor = 7,
                                            xlab = "Calibration period transactions",
                                            ylab = "Customers",
                                            title = "Frequency of Repeat Transactions",
                                            sample_size = 1000) {

  # actual
  x_act <- cal.cbs$x
  x_act[x_act > censor] <- censor
  x_act <- table(x_act)

  # expected
  x_est <- sapply(unique(cal.cbs$T.cal), function(tcal) {
    n <- sum(cal.cbs$T.cal == tcal)
    prop <- mcmc.pmf(draws, t = tcal, x = 0:(censor - 1), sample_size = sample_size)
    prop <- c(prop, 1 - sum(prop))
    prop * (n / nrow(cal.cbs))
  })
  x_est <- apply(x_est, 1, sum) * nrow(cal.cbs)

  mat <- matrix(c(x_act, x_est), nrow = 2, ncol = censor + 1, byrow = TRUE)
  rownames(mat) <- c("n.x.actual", "n.x.expected")
  colnames(mat) <- c(0:(censor - 1), paste0(censor, "+"))

  barplot(mat, beside = TRUE, col = 1:2, main = title, xlab = xlab, ylab = ylab, ylim = c(0, max(x_act) * 1.1))
  legend("topright", legend = c("Actual", "Model"), col = 1:2, lty = 1:2, lwd = 1, xjust = 1)

  colnames(mat) <- paste0("freq.", colnames(mat))
  mat
}
