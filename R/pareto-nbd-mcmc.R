
#' Pareto/NBD (HB) Parameter Draws
#'
#' Returns draws from the posterior distributions of the Pareto/NBD (HB)
#' parameters, on cohort as well as on customer level.
#'
#' See \code{demo('pareto-ggg')} for how to apply this model.
#'
#' method 1) If \code{use_data_augmentation==TRUE} MCMC scheme takes advantage of
#' conjugate priors for drawing lambda and mu, by augmenting the parameter space
#' with unobserved lifetime 'tau' and activity status 'z'. See technical appendix
#' to (Abe 2009).
#'
#' method 2) If \code{use_data_augmentation==FALSE} then implementation follows
#' Shao-Hui Ma & Jin-Lan Liu paper
#' \url{http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=4344404}, i.e. no data
#' augmentation and draws on individual level need to be done via slice
#' sampling. As such it is 10x slower than method 1)
#'
#' @param cal.cbs Calibration period customer-by-sufficient-statistic (CBS)
#'   data.frame. It must contain a row for each customer, and columns \code{x}
#'   for frequency, \code{t.x} for recency and \code{T.cal} for the total time
#'   observed. A correct format can be easily generated based on the complete
#'   event log of a customer cohort with \code{\link{elog2cbs}}.
#' @param mcmc Number of MCMC steps.
#' @param burnin Number of initial MCMC steps which are discarded.
#' @param thin Only every \code{thin}-th MCMC step will be returned.
#' @param chains Number of MCMC chains to be run.
#' @param mc.cores Number of cores to use in parallel (Unix only). Defaults to \code{min(chains, detectCores())}.
#' @param use_data_augmentation determines MCMC method to be used
#' @param param_init List of start values for cohort-level parameters.
#' @param trace Print logging statement every \code{trace}-th iteration. Not available for \code{mc.cores > 1}.
#' @return 2-element list:
#' \itemize{
#'  \item{\code{level_1 }}{list of \code{\link{mcmc.list}}s, one for each customer, with draws for customer-level parameters \code{lambda}, \code{tau}, \code{z}, \code{mu}}
#'  \item{\code{level_2 }}{\code{\link{mcmc.list}}, with draws for cohort-level parameters \code{r}, \code{alpha}, \code{s}, \code{beta}}
#' }
#' @export
#' @seealso \code{\link{pnbd.GenerateData} } \code{\link{mcmc.DrawFutureTransactions} } \code{\link{mcmc.PAlive} }
#' @references Ma, Shao-Hui, and Jin-Lan Liu. 'The MCMC approach for solving the Pareto/NBD model and possible extensions.' Natural Computation, 2007. ICNC 2007. Third International Conference on. Vol. 2. IEEE, 2007. \url{http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=4344404}
#' @references Abe, Makoto. 'Counting your customers one by one: A hierarchical Bayes extension to the Pareto/NBD model.' Marketing Science 28.3 (2009): 541-553.
#' @examples
#' data("groceryElog")
#' cbs <- elog2cbs(groceryElog, T.cal = "2006-12-31")
#' param.draws <- pnbd.mcmc.DrawParameters(cbs,
#'   mcmc = 200, burnin = 100, thin = 20, chains = 1) # short MCMC to run demo fast
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
pnbd.mcmc.DrawParameters <- function(cal.cbs, mcmc = 2500, burnin = 500, thin = 50, chains = 2, mc.cores = NULL,
  use_data_augmentation = TRUE, param_init = NULL, trace = 100) {

  # ** methods to sample heterogeneity parameters {r, alpha, s, beta} **

  draw_gamma_params <- function(type, level_1, level_2, hyper_prior) {
    if (type == "lambda") {
      x <- level_1["lambda", ]
      cur_params <- c(level_2["r"], level_2["alpha"])
      hyper <- unlist(hyper_prior[c("r_1", "r_2", "alpha_1", "alpha_2")])
    } else if (type == "mu") {
      x <- level_1["mu", ]
      cur_params <- c(level_2["s"], level_2["beta"])
      hyper <- unlist(hyper_prior[c("s_1", "s_2", "beta_1", "beta_2")])
    }
    slice_sample_gamma_parameters(x, cur_params, hyper, steps = 50, w = 0.1)
  }

  # ** methods to sample individual-level parameters (with data augmentation) **

  draw_lambda <- function(data, level_1, level_2) {
    N <- nrow(data)
    x <- data$x
    T.cal <- data$T.cal
    tau <- level_1["tau", ]
    r <- level_2["r"]
    alpha <- level_2["alpha"]

    lambda <- rgamma(n = N, shape = r + x, rate = alpha + pmin(tau, T.cal))
    lambda[lambda == 0 | log(lambda) < -30] <- exp(-30)  # avoid numeric overflow
    return(lambda)
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

  draw_tau <- function(data, level_1) {
    N <- nrow(data)
    tx <- data$t.x
    Tcal <- data$T.cal
    lambda <- level_1["lambda", ]
    mu <- level_1["mu", ]

    mu_lam <- mu + lambda
    t_diff <- Tcal - tx

    # sample z
    p_alive <- 1 / (1 + (mu / mu_lam) * (exp(mu_lam * t_diff) - 1))
    alive <- p_alive > runif(n = N)

    # sample tau
    tau <- numeric(N)

    # Case: still alive - left truncated exponential distribution -> [Tcal, Inf]
    if (any(alive)) {
      tau[alive] <- Tcal[alive] + rexp(sum(alive), mu[alive])
    }

    # Case: churned - double truncated exponential distribution -> [tx, Tcal]
    if (any(!alive)) {
      mu_lam_tx <- pmin(700, mu_lam[!alive] * tx[!alive])
      mu_lam_Tcal <- pmin(700, mu_lam[!alive] * Tcal[!alive])
      # sample with http://en.wikipedia.org/wiki/Inverse_transform_sampling
      rand <- runif(n = sum(!alive))
      tau[!alive] <- -log((1 - rand) * exp(-mu_lam_tx) + rand * exp(-mu_lam_Tcal)) / mu_lam[!alive] # nolint
    }

    return(tau)
  }

  # ** methods to sample individual-level parameters (without data augmentation) **

  draw_lambda_ma_liu <- function(data, level_1, level_2) {
    slice_sample_ma_liu("lambda",
                        x = data$x, tx = data$t.x, Tcal = data$T.cal,
                        lambda = level_1["lambda", ], mu = level_1["mu", ],
                        r = level_2["r"], alpha = level_2["alpha"],
                        s = level_2["s"], beta = level_2["beta"])
  }

  draw_mu_ma_liu <- function(data, level_1, level_2) {
    slice_sample_ma_liu("mu",
                        x = data$x, tx = data$t.x, Tcal = data$T.cal,
                        lambda = level_1["lambda", ], mu = level_1["mu", ],
                        r = level_2["r"], alpha = level_2["alpha"],
                        s = level_2["s"], beta = level_2["beta"])
  }

  run_single_chain <- function(chain_id = 1, data, hyper_prior) {

    ## initialize arrays for storing draws ##

    nr_of_cust <- nrow(data)
    nr_of_draws <- (mcmc - 1) %/% thin + 1
    level_2_draws <- array(NA_real_, dim = c(nr_of_draws, 4))
    dimnames(level_2_draws)[[2]] <- c("r", "alpha", "s", "beta")
    level_1_draws <- array(NA_real_, dim = c(nr_of_draws, 4, nr_of_cust))
    dimnames(level_1_draws)[[2]] <- c("lambda", "mu", "tau", "z")

    ## initialize parameters ##

    level_2 <- level_2_draws[1, ]
    level_2["r"] <- param_init$r
    level_2["alpha"] <- param_init$alpha
    level_2["s"] <- param_init$s
    level_2["beta"] <- param_init$beta

    level_1 <- level_1_draws[1, , ] # nolint
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
      draw_lambda <- if (use_data_augmentation)
        draw_lambda else draw_lambda_ma_liu
      draw_mu <- if (use_data_augmentation)
        draw_mu else draw_mu_ma_liu
      level_1["lambda", ] <- draw_lambda(data, level_1, level_2)
      level_1["mu", ] <- draw_mu(data, level_1, level_2)
      level_1["tau", ] <- draw_tau(data, level_1)
      level_1["z", ] <- as.numeric(level_1["tau", ] > data$T.cal)

      # draw heterogeneity parameters
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
                      beta_1 = 0.001, beta_2 = 0.001)

  # set param_init (if not passed as argument)
  if (is.null(param_init)) {
    try({
        df <- cal.cbs[sample(nrow(cal.cbs), min(nrow(cal.cbs), 1000)), ]
        param_init <- BTYD::pnbd.EstimateParameters(df)
        names(param_init) <- c("r", "alpha", "s", "beta")
        param_init <- as.list(param_init)
      },
      silent = TRUE)
    if (is.null(param_init))
      param_init <- list(r = 1, alpha = 1, s = 1, beta = 1)
    cat("set param_init:", paste(round(unlist(param_init), 4), collapse = ", "), "\n")
  }

  # check whether input data meets requirements
  stopifnot(is.data.frame(cal.cbs))
  stopifnot(all(c("x", "t.x", "T.cal") %in% names(cal.cbs)))
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


#' Simulate data according to Pareto/NBD model assumptions
#'
#' @param n Number of customers.
#' @param T.cal Length of calibration period. If a vector is provided, then it
#'   is assumed that customers have different 'birth' dates, i.e.
#'   \eqn{max(T.cal)-T.cal}.
#' @param T.star Length of holdout period. This may be a vector.
#' @param params A list of model parameters \code{r},
#'   \code{alpha}, \code{s}, \code{beta}.
#' @param date.zero Initial date for cohort start. Can be of class character, Date or POSIXt.
#' @return List of length 2:
#' \item{\code{cbs}}{A data.frame with a row for each customer and the summary statistic as columns.}
#' \item{\code{elog}}{A data.frame with a row for each transaction, and columns \code{cust}, \code{date} and \code{t}.}
#' @export
#' @examples
#' params <- list(r = 5, alpha = 10, s = 0.8, beta = 12)
#' data <- pnbd.GenerateData(n = 1000, T.cal = 32, T.star = 32, params)
#' cbs <- data$cbs  # customer by sufficient summary statistic - one row per customer
#' elog <- data$elog  # Event log - one row per event/purchase
pnbd.GenerateData <- function(n, T.cal, T.star, params, date.zero = "2000-01-01") {
  params$k <- 1
  pggg.GenerateData(n, T.cal, T.star, params, date.zero)
}
