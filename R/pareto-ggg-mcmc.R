
#' Pareto/GGG Parameter Draws
#'
#' Returns draws from the posterior distributions of the Pareto/GGG
#' parameters, on cohort as well as on customer level.
#'
#' See \code{demo('pareto-ggg')} for how to apply this model.
#'
#' @param cal.cbs Calibration period customer-by-sufficient-statistic (CBS)
#'   data.frame. It must contain a row for each customer, and columns \code{x}
#'   for frequency, \code{t.x} for recency , \code{T.cal} for the total time
#'   observed, as well as the sum over logarithmic intertransaction times
#'   \code{litt}. A correct format can be easily generated based on the complete
#'   event log of a customer cohort with \code{\link{elog2cbs}}.
#' @param mcmc Number of MCMC steps.
#' @param burnin Number of initial MCMC steps which are discarded.
#' @param thin Only every \code{thin}-th MCMC step will be returned.
#' @param chains Number of MCMC chains to be run.
#' @param mc.cores Number of cores to use in parallel (Unix only). Defaults to \code{min(chains, detectCores())}.
#' @param param_init List of start values for cohort-level parameters.
#' @param trace Print logging statement every \code{trace}-th iteration. Not available for \code{mc.cores > 1}.
#' @return List of length 2:
#' \item{\code{level_1}}{list of \code{\link{mcmc.list}}s, one for each customer, with draws for customer-level parameters \code{k}, \code{lambda}, \code{tau}, \code{z}, \code{mu}}
#' \item{\code{level_2}}{\code{\link{mcmc.list}}, with draws for cohort-level parameters \code{r}, \code{alpha}, \code{s}, \code{beta}, \code{t}, \code{gamma}}
#' @export
#' @references Platzer, M., & Reutterer, T. (2016). Ticking away the moments:
#'   Timing regularity helps to better predict customer activity. Marketing
#'   Science, 35(5), 779-799. \doi{10.1287/mksc.2015.0963}
#' @seealso \code{\link{pggg.GenerateData} } \code{\link{mcmc.PAlive} } \code{\link{mcmc.DrawFutureTransactions} }
#' @examples
#' data("groceryElog")
#' cbs <- elog2cbs(groceryElog, T.cal = "2006-12-31")
#' param.draws <- pggg.mcmc.DrawParameters(cbs,
#'   mcmc = 20, burnin = 10, thin = 2, chains = 1) # short MCMC to run demo fast
#'
#' # cohort-level parameter draws
#' as.matrix(param.draws$level_2)
#' # customer-level parameter draws for customer with ID '4'
#' as.matrix(param.draws$level_1[["4"]])
#'
#' # estimate future transactions
#' xstar.draws <- mcmc.DrawFutureTransactions(cbs, param.draws, cbs$T.star)
#' xstar.est <- apply(xstar.draws, 2, mean)
#' head(xstar.est)
pggg.mcmc.DrawParameters <- function(cal.cbs, mcmc = 2500, burnin = 500, thin = 50, chains = 2, mc.cores = NULL,
  param_init = NULL, trace = 100) {

  # ** methods to sample heterogeneity parameters {r, alpha, s, beta, t, gamma} **

  draw_gamma_params <- function(type, level_1, level_2, hyper_prior) {
    if (type == "lambda") {
      x <- level_1["lambda", ]
      cur_params <- c(level_2["r"], level_2["alpha"])
      hyper <- unlist(hyper_prior[c("r_1", "r_2", "alpha_1", "alpha_2")])
    } else if (type == "mu") {
      x <- level_1["mu", ]
      cur_params <- c(level_2["s"], level_2["beta"])
      hyper <- unlist(hyper_prior[c("s_1", "s_2", "beta_1", "beta_2")])
    } else if (type == "k") {
      x <- level_1["k", ]
      cur_params <- c(level_2["t"], level_2["gamma"])
      hyper <- unlist(hyper_prior[c("t_1", "t_2", "gamma_1", "gamma_2")])
    }
    slice_sample_gamma_parameters(x, cur_params, hyper, steps = 200, w = 0.1)
  }

  # ** methods to sample individual-level parameters **

  draw_k <- function(data, level_1, level_2) {
    pggg_slice_sample("k",
                      x = data$x, tx = data$t.x, Tcal = data$T.cal, litt = data$litt,
                      k = level_1["k", ], lambda = level_1["lambda", ],
                      mu = level_1["mu", ], tau = level_1["tau", ],
                      t = level_2["t"], gamma = level_2["gamma"],
                      r = level_2["r"], alpha = level_2["alpha"],
                      s = level_2["s"], beta = level_2["beta"])
  }

  draw_lambda <- function(data, level_1, level_2) {
    pggg_slice_sample("lambda",
                      x = data$x, tx = data$t.x, Tcal = data$T.cal, litt = data$litt,
                      k = level_1["k", ], lambda = level_1["lambda", ],
                      mu = level_1["mu", ], tau = level_1["tau", ],
                      t = level_2["t"], gamma = level_2["gamma"],
                      r = level_2["r"], alpha = level_2["alpha"],
                      s = level_2["s"], beta = level_2["beta"])
  }

  draw_mu <- function(data, level_1, level_2) {
    N <- nrow(data)
    tau <- level_1["tau", ]
    s <- level_2["s"]
    beta <- level_2["beta"]

    mu <- rgamma(n = N, shape = s + 1, rate = beta + tau)
    mu[mu == 0 | log(mu) < -30] <- exp(-30)  # avoid numeric overflow
    return(mu)
  }

  draw_tau <- function(data, level_1, level_2) {
    N <- nrow(data)
    x <- data$x
    tx <- data$t.x
    Tcal <- data$T.cal
    lambda <- level_1["lambda", ]
    k <- level_1["k", ]
    mu <- level_1["mu", ]

    # sample z
    p_alive <- pggg_palive(x, tx, Tcal, k, lambda, mu)
    alive <- p_alive > runif(n = N)

    # sample tau
    tau <- numeric(N)

    # Case: still alive - left truncated exponential distribution -> [Tcal, Inf]
    if (any(alive)) {
      tau[alive] <- Tcal[alive] + rexp(sum(alive), mu[alive])
    }

    # Case: churned - distribution of tau truncated to [tx, pmin(tx+1, Tcal)]
    if (any(!alive)) {
      tau[!alive] <- pggg_slice_sample("tau", x = data$x[!alive], tx = data$t.x[!alive], Tcal = data$T.cal[!alive],
        litt = data$litt[!alive], k = level_1["k", !alive], lambda = level_1["lambda", !alive], mu = level_1["mu",
          !alive], tau = level_1["tau", !alive], t = level_2["t"], gamma = level_2["gamma"], r = level_2["r"],
        alpha = level_2["alpha"], s = level_2["s"], beta = level_2["beta"])
    }

    return(tau)
  }

  run_single_chain <- function(chain_id, data, hyper_prior) {

    ## initialize arrays for storing draws ##

    nr_of_cust <- nrow(data)
    nr_of_draws <- (mcmc - 1) %/% thin + 1
    level_2_draws <- array(NA_real_, dim = c(nr_of_draws, 6))
    dimnames(level_2_draws)[[2]] <- c("t", "gamma", "r", "alpha", "s", "beta")
    level_1_draws <- array(NA_real_, dim = c(nr_of_draws, 5, nr_of_cust))
    dimnames(level_1_draws)[[2]] <- c("k", "lambda", "mu", "tau", "z")

    ## initialize parameters ##

    level_2 <- level_2_draws[1, ]
    level_2["t"] <- param_init$t
    level_2["gamma"] <- param_init$gamma
    level_2["r"] <- param_init$r
    level_2["alpha"] <- param_init$alpha
    level_2["s"] <- param_init$s
    level_2["beta"] <- param_init$beta

    level_1 <- level_1_draws[1, , ] # nolint
    level_1["k", ] <- 1
    level_1["lambda", ] <- mean(data$x) / mean(ifelse(data$t.x == 0, data$T.cal, data$t.x))
    level_1["tau", ] <- data$t.x + 0.5 / level_1["lambda", ]
    level_1["z", ] <- as.numeric(level_1["tau", ] > data$T.cal)
    level_1["mu", ] <- 1 / level_1["tau", ]

    ## run MCMC chain ##

    for (step in 1:(burnin + mcmc)) {
      if (step %% trace == 0)
        cat("chain:", chain_id, "step:", step, "of", (burnin + mcmc), "\n")

      # store
      if ( (step - burnin) > 0 & (step - 1 - burnin) %% thin == 0) {
        idx <- (step - 1 - burnin) %/% thin + 1
        level_1_draws[idx, , ] <- level_1 # nolint
        level_2_draws[idx, ] <- level_2
      }

      # draw individual-level parameters
      level_1["k", ] <- draw_k(data, level_1, level_2)
      level_1["lambda", ] <- draw_lambda(data, level_1, level_2)
      level_1["mu", ] <- draw_mu(data, level_1, level_2)
      level_1["tau", ] <- draw_tau(data, level_1, level_2)
      level_1["z", ] <- as.numeric(level_1["tau", ] > data$T.cal)

      # draw heterogeneity parameters
      level_2[c("t", "gamma")] <- draw_gamma_params("k", level_1, level_2, hyper_prior)
      level_2[c("r", "alpha")] <- draw_gamma_params("lambda", level_1, level_2, hyper_prior)
      level_2[c("s", "beta")] <- draw_gamma_params("mu", level_1, level_2, hyper_prior)
    }

    # convert MCMC draws into coda::mcmc objects
    return(list(
      "level_1" = lapply(1:nr_of_cust,
                         function(i) mcmc(level_1_draws[, , i], start = burnin, thin = thin)), # nolint
      "level_2" = mcmc(level_2_draws, start = burnin, thin = thin)))
  }

  # set hyper priors
  hyper_prior <- list(r_1 = 0.001, r_2 = 0.001,
                      alpha_1 = 0.001, alpha_2 = 0.001,
                      s_1 = 0.001, s_2 = 0.001,
                      beta_1 = 0.001, beta_2 = 0.001,
                      t_1 = 0.001, t_2 = 0.001,
                      gamma_1 = 0.001, gamma_2 = 0.001)

  # set param_init (if not passed as argument)
  if (is.null(param_init)) {
    try({
        df <- cal.cbs[sample(nrow(cal.cbs), min(nrow(cal.cbs), 1000)), ]
        param_init <- c(1, 1, BTYD::pnbd.EstimateParameters(df))
        names(param_init) <- c("t", "gamma", "r", "alpha", "s", "beta")
        param_init <- as.list(param_init)
      },
      silent = TRUE)
    if (is.null(param_init))
      param_init <- list(t = 1, gamma = 1, r = 1, alpha = 1, s = 1, beta = 1)
    cat("set param_init:", paste(round(unlist(param_init), 4), collapse = ", "), "\n")
  }

  # check whether input data meets requirements
  stopifnot(is.data.frame(cal.cbs))
  stopifnot(all(c("x", "t.x", "T.cal", "litt") %in% names(cal.cbs)))
  stopifnot(all(is.finite(cal.cbs$litt)))

  # run multiple chains - executed in parallel on Unix
  ncores <- ifelse(!is.null(mc.cores), min(chains, mc.cores), ifelse(.Platform$OS.type == "windows", 1, min(chains,
    detectCores())))
  if (ncores > 1)
    cat("running in parallel on", ncores, "cores\n")
  draws <- mclapply(1:chains, function(i) run_single_chain(i, cal.cbs, hyper_prior), mc.cores = ncores)

  # merge chains into code::mcmc.list objects
  out <- list(level_1 = lapply(1:nrow(cal.cbs), function(i) mcmc.list(lapply(draws, function(draw) draw$level_1[[i]]))),
    level_2 = mcmc.list(lapply(draws, function(draw) draw$level_2)))
  if ("cust" %in% names(cal.cbs))
    names(out$level_1) <- cal.cbs$cust
  return(out)
}


#' Pareto/GGG Plot Regularity Rate Heterogeneity
#'
#' Plots and returns the estimated gamma distribution of k (customers'
#' regularity in interpurchase times).
#'
#' @param draws MCMC draws as returned by \code{\link{pggg.mcmc.DrawParameters}}.
#' @param xmax Upper bound for x-scale.
#' @param fn Optional function to summarize individual-level draws for k, e.g. 'mean'.
#' @param title Plot title.
#'
#' @references Platzer, M., & Reutterer, T. (2016). Ticking away the moments:
#'   Timing regularity helps to better predict customer activity. Marketing
#'   Science, 35(5), 779-799. \doi{10.1287/mksc.2015.0963}
#' @export
#' @examples
#' data("groceryElog")
#' cbs <- elog2cbs(groceryElog, T.cal = "2006-12-31")
#' param.draws <- pggg.mcmc.DrawParameters(cbs,
#'   mcmc = 20, burnin = 10, thin = 2, chains = 1) # short MCMC to run demo fast
#' pggg.plotRegularityRateHeterogeneity(param.draws)
pggg.plotRegularityRateHeterogeneity <- function(draws, xmax = NULL, fn = NULL,
                                                 title = "Distribution of Regularity Rate k") {
  stopifnot("k" %in% colnames(as.matrix(draws$level_1[[1]])))
  ks <- sapply(draws$level_1, function(draw) as.matrix(draw[, "k"]))
  if (!is.null(fn))
    ks <- apply(ks, 2, fn)
  if (is.null(xmax))
    xmax <- min(10, quantile(ks, 0.95) * 1.5)
  mar_top <- ifelse(title != "", 2.5, 1)
  op <- par(mar = c(2.5, 2.5, mar_top, 2.5))
  plot(density(ks, from = 0), xlim = c(0, xmax), main = title, xlab = "", ylab = "", frame = FALSE)
  abline(v = 1, lty = 3)
  abline(v = median(ks), col = "red")
  par(op)
  invisible()
}


#' Simulate data according to Pareto/GGG model assumptions
#'
#' @param n Number of customers.
#' @param T.cal Length of calibration period. If a vector is provided, then it
#'   is assumed that customers have different 'birth' dates, i.e.
#'   \eqn{max(T.cal)-T.cal}.
#' @param T.star Length of holdout period. This may be a vector.
#' @param params A list of model parameters \code{r},
#'   \code{alpha}, \code{s}, \code{beta}, \code{t} and \code{gamma}.
#' @param date.zero Initial date for cohort start. Can be of class character, Date or POSIXt.
#' @return List of length 2:
#' \item{\code{cbs}}{A data.frame with a row for each customer and the summary statistic as columns.}
#' \item{\code{elog}}{A data.frame with a row for each transaction, and columns \code{cust}, \code{date} and \code{t}.}
#' @export
#' @references Platzer, M., & Reutterer, T. (2016). Ticking away the moments:
#'   Timing regularity helps to better predict customer activity. Marketing
#'   Science, 35(5), 779-799. \doi{10.1287/mksc.2015.0963}
#' @examples
#' params <- list(t = 4.5, gamma = 1.5, r = 5, alpha = 10, s = 0.8, beta = 12)
#' data <- pggg.GenerateData(n = 200, T.cal = 32, T.star = 32, params)
#' cbs <- data$cbs  # customer by sufficient summary statistic - one row per customer
#' elog <- data$elog  # Event log - one row per event/purchase
pggg.GenerateData <- function(n, T.cal, T.star, params, date.zero = "2000-01-01") {

  # set start date for each customer, so that they share same T.cal date
  T.cal.fix <- max(T.cal)
  T.cal <- rep(T.cal, length.out = n)
  T.zero <- T.cal.fix - T.cal
  date.zero <- as.POSIXct(date.zero)

  # sample regularity parameter k for each customer
  if (all(c("t", "gamma") %in% names(params))) {
    # Case A: regularity parameter k is gamma-distributed across customers
    ks <- rgamma(n, shape = params$t, rate = params$gamma)
    ks <- pmax(0.1, ks)  # ensure that k is not too small, otherwise itt can be 0

  } else if ("k" %in% names(params)) {
    # Case B: regularity parameter k is fixed across customers
    ks <- rep(params$k, n)

  } else {
    # Case C: k=1 is assumed, i.e. Pareto/NBD
    ks <- rep(1, n)
  }

  # sample intertransaction timings parameter lambda for each customer
  lambdas <- rgamma(n, shape = params$r, rate = params$alpha)

  # sample lifetime for each customer
  mus <- rgamma(n, shape = params$s, rate = params$beta)
  taus <- rexp(n, rate = mus)

  # sample intertransaction timings
  elog_list <- lapply(1:n, function(i) {
    # sample 'sufficiently' large amount of inter-transaction times
    minT <- min(T.cal[i] + max(T.star), taus[i])
    itt_draws <- max(10, round(minT * lambdas[i] * 1.5))
    itt_fn <- function(n) rgamma(n, shape = ks[i], rate = ks[i] * lambdas[i])
    itts <- itt_fn(itt_draws)
    if (sum(itts) < minT) itts <- c(itts, itt_fn(itt_draws * 4))
    if (sum(itts) < minT) itts <- c(itts, itt_fn(itt_draws * 800))
    if (sum(itts) < minT) stop("not enough inter-transaction times sampled: ", sum(itts), " < ", minT)
    ts <- cumsum(c(0, itts))
    ts <- ts[ts <= taus[i]] # trim to lifetime
    ts <- T.zero[i] + ts # shift by T_0
    ts <- ts[ts <= (T.cal.fix + max(T.star))] # trim to observation length
    return(ts)
  })

  # build elog
  elog <- data.table("cust" = rep(1:n, sapply(elog_list, length)), "t" = unlist(elog_list))
  elog[["date"]] <- date.zero + elog[["t"]] * 3600 * 24 * 7

  # build cbs
  date.cal <- date.zero + T.cal.fix * 3600 * 24 * 7
  date.tot <- date.cal + T.star * 3600 * 24 * 7
  cbs <- elog2cbs(elog, T.cal = date.cal)
  if (length(T.star) == 1) set(cbs, j = "T.star", value = T.star[1])
  xstar.cols <- if (length(T.star) == 1) "x.star" else paste0("x.star", T.star)
  for (j in 1:length(date.tot)) {
    set(cbs, j = xstar.cols[j],
        value = sapply(elog_list, function(t) sum(t > T.cal.fix & t <= T.cal.fix + T.star[j])))
  }
  set(cbs, j = "k", value = ks)
  set(cbs, j = "lambda", value = lambdas)
  set(cbs, j = "mu", value = mus)
  set(cbs, j = "tau", value = taus)
  set(cbs, j = "alive", value = (T.zero + taus) > T.cal.fix)

  return(list("cbs" = setDF(cbs), "elog" = setDF(elog)))
}
