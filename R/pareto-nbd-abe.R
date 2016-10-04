
#' Pareto/NBD (Abe) Parameter Draws
#'
#' Returns draws from the posterior distributions of the Pareto/NBD (Abe)
#' parameters, on cohort as well as on customer level.
#'
#' See \code{demo('pareto-abe')} for how to apply this model.
#'
#' @param cal.cbs Calibration period customer-by-sufficient-statistic (CBS)
#'   data.frame. It must contain a row for each customer, and columns \code{x}
#'   for frequency, \code{t.x} for recency and \code{T.cal} for the total time
#'   observed. A correct format can be easily generated based on the complete
#'   event log of a customer cohort with \code{\link{elog2cbs}}.
#' @param covariates A vector of columns of \code{cal.cbs} which contain customer-level covariates.
#' @param mcmc Number of MCMC steps.
#' @param burnin Number of initial MCMC steps which are discarded.
#' @param thin Only every \code{thin}-th MCMC step will be returned.
#' @param chains Number of MCMC chains to be run.
#' @param mc.cores Number of cores to use in parallel (Unix only). Defaults to \code{min(chains, detectCores())}.
#' @param trace Print logging statement every \code{trace}-th iteration. Not available for \code{mc.cores > 1}.
#' @return List of length 2:
#' \item{\code{level_1}}{list of \code{\link{mcmc.list}}s, one for each customer, with draws for customer-level parameters \code{k}, \code{lambda}, \code{tau}, \code{z}, \code{mu}}
#' \item{\code{level_2}}{\code{\link{mcmc.list}}, with draws for cohort-level parameters}
#' @export
#' @seealso \code{\link{abe.GenerateData} } \code{\link{mcmc.PAlive} } \code{\link{mcmc.DrawFutureTransactions} }
#' @references Abe, Makoto. 'Counting your customers one by one: A hierarchical Bayes extension to the Pareto/NBD model.' Marketing Science 28.3 (2009): 541-553.
#' @examples
#' data("groceryElog")
#' cbs <- elog2cbs(groceryElog, T.cal = "2006-12-31")
#' cbs$cov1 <- as.integer(cbs$cust) %% 2 # create dummy covariate
#' param.draws <- abe.mcmc.DrawParameters(cbs, c("cov1"),
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
abe.mcmc.DrawParameters <- function(cal.cbs, covariates = c(), mcmc = 2500, burnin = 500, thin = 50, chains = 2,
  mc.cores = NULL, trace = 100) {

  # ** methods to sample heterogeneity parameters {beta, gamma} **

  draw_level_2 <- function(covars, level_1, hyper_prior) {
    # standard multi-variate normal regression update
    draw <- bayesm::rmultireg(Y = log(t(level_1[c("lambda", "mu"), ])),
                              X = covars,
                              Bbar = hyper_prior$beta_0,
                              A = hyper_prior$A_0,
                              nu = hyper_prior$nu_00,
                              V = hyper_prior$gamma_00)
    return(list(beta = t(draw$B), gamma = draw$Sigma))
  }

  # ** methods to sample individual-level parameters **

  draw_z <- function(data, level_1) {
    tx <- data$t.x
    Tcal <- data$T.cal
    lambda <- level_1["lambda", ]
    mu <- level_1["mu", ]

    mu_lam <- mu + lambda
    t_diff <- Tcal - tx

    prob <- 1 / (1 + (mu / mu_lam) * (exp(mu_lam * t_diff) - 1))
    z <- as.numeric(runif(length(prob)) < prob)
    return(z)
  }

  draw_tau <- function(data, level_1) {
    N <- nrow(data)
    tx <- data$t.x
    Tcal <- data$T.cal
    lambda <- level_1["lambda", ]
    mu <- level_1["mu", ]
    mu_lam <- mu + lambda
    z <- level_1["z", ]

    alive <- z == 1
    tau <- numeric(N)

    # Case: still alive - left truncated exponential distribution -> [T.cal, Inf]
    if (any(alive)) {
      tau[alive] <- Tcal[alive] + rexp(sum(alive), mu[alive])
    }

    # Case: churned - double truncated exponential distribution -> [tx, T.cal]
    if (any(!alive)) {
      mu_lam_tx <- pmin(700, mu_lam[!alive] * tx[!alive])
      mu_lam_Tcal <- pmin(700, mu_lam[!alive] * Tcal[!alive])
      rand <- runif(n = sum(!alive))
      tau[!alive] <- -log( (1 - rand) * exp(-mu_lam_tx) + rand * exp(-mu_lam_Tcal)) / mu_lam[!alive]
    }
    return(tau)
  }

  draw_level_1 <- function(data, covars, level_1, level_2) {
    # sample (lambda, mu) given (z, tau, beta, gamma)
    N <- nrow(data)
    x <- data$x
    Tcal <- data$T.cal
    z <- level_1["z", ]
    tau <- level_1["tau", ]
    mvmean <- covars[, ] %*% t(level_2$beta)
    gamma <- level_2$gamma
    inv_gamma <- solve(gamma)

    cur_lambda <- level_1["lambda", ]
    cur_mu <- level_1["mu", ]

    log_post <- function(log_theta) {
      log_lambda <- log_theta[1, ]
      log_mu <- log_theta[2, ]
      diff_lambda <- log_lambda - mvmean[, 1]
      diff_mu <- log_mu - mvmean[, 2]
      likel <- x * log_lambda + (1 - z) * log_mu - (exp(log_lambda) + exp(log_mu)) * (z * Tcal + (1 - z) *
        tau)
      prior <- -0.5 * (diff_lambda ^ 2 * inv_gamma[1, 1] +
                         2 * diff_lambda * diff_mu * inv_gamma[1, 2] +
                         diff_mu ^ 2 * inv_gamma[2, 2])
      post <- likel + prior
      post[log_mu > 5] <- -Inf  # cap !!
      return(post)
    }

    # current state
    cur_log_theta <- rbind(log(cur_lambda), log(cur_mu))
    cur_post <- log_post(cur_log_theta)

    step <- function(cur_log_theta, cur_post) {
      # new proposal
      new_log_theta <- cur_log_theta + rbind(gamma[1, 1] * rt(N, df = 3), gamma[2, 2] * rt(n = N, df = 3))
      new_log_theta[1, ] <- pmax(pmin(new_log_theta[1, ], 70), -70)
      new_log_theta[2, ] <- pmax(pmin(new_log_theta[2, ], 70), -70)
      new_post <- log_post(new_log_theta)

      # accept/reject new proposal
      mhratio <- exp(new_post - cur_post)
      accepted <- mhratio > runif(n = N)

      cur_log_theta[, accepted] <- new_log_theta[, accepted]
      cur_post[accepted] <- new_post[accepted]

      list(cur_log_theta = cur_log_theta, cur_post = cur_post)
    }

    iter <- 1  # how high do we need to set this? 1/5/10/100?
    for (i in 1:iter) {
      draw <- step(cur_log_theta, cur_post)
      cur_log_theta <- draw$cur_log_theta
      cur_post <- draw$cur_post
    }
    cur_theta <- exp(cur_log_theta)

    return(list(lambda = cur_theta[1, ], mu = cur_theta[2, ]))
  }


  run_single_chain <- function(chain_id, data, hyper_prior) {

    ## initialize arrays for storing draws ##
    nr_of_cust <- nrow(data)
    nr_of_draws <- (mcmc - 1) %/% thin + 1

    level_1_draws <- array(NA_real_, dim = c(nr_of_draws, 4, nr_of_cust))
    dimnames(level_1_draws)[[2]] <- c("lambda", "mu", "tau", "z")

    level_2_draws <- array(NA_real_, dim = c(nr_of_draws, 2 * K + 3))
    nm <- c("log_lambda", "log_mu")
    if (K > 1)
      nm <- paste(rep(nm, times = K), rep(colnames(covars), each = 2), sep = "_")
    dimnames(level_2_draws)[[2]] <- c(nm, "var_log_lambda", "cov_log_lambda_log_mu", "var_log_mu")

    ## initialize parameters ##

    level_1 <- level_1_draws[1, , ] # nolint
    level_1["lambda", ] <- mean(data$x) / mean(ifelse(data$t.x == 0, data$T.cal, data$t.x))
    level_1["mu", ] <- 1 / (data$t.x + 0.5 / level_1["lambda", ])

    ## run MCMC chain ##

    hyper_prior$beta_0[1, "log_lambda"] <- log(mean(level_1["lambda", ]))
    hyper_prior$beta_0[1, "log_mu"] <- log(mean(level_1["mu", ]))

    for (step in 1:(burnin + mcmc)) {
      if (step %% trace == 0)
        cat("chain:", chain_id, "step:", step, "of", (burnin + mcmc), "\n")

      # draw individual-level parameters
      level_1["z", ] <- draw_z(data, level_1)
      level_1["tau", ] <- draw_tau(data, level_1)

      level_2 <- draw_level_2(covars, level_1, hyper_prior)

      draw <- draw_level_1(data, covars, level_1, level_2)
      level_1["lambda", ] <- draw$lambda
      level_1["mu", ] <- draw$mu

      # store
      if ( (step - burnin) > 0 & (step - 1 - burnin) %% thin == 0) {
        idx <- (step - 1 - burnin) %/% thin + 1
        level_1_draws[idx, , ] <- level_1 # nolint
        level_2_draws[idx, ] <- c(level_2$beta, level_2$gamma[1, 1], level_2$gamma[1, 2], level_2$gamma[2,
          2])
      }
    }

    # convert MCMC draws into coda::mcmc objects
    return(list(
      "level_1" = lapply(1:nr_of_cust,
                         function(i) mcmc(level_1_draws[, , i], start = burnin, thin = thin)), # nolint
      "level_2" = mcmc(level_2_draws, start = burnin, thin = thin)))
  }

  # check whether input data meets requirements
  stopifnot(is.data.frame(cal.cbs))
  stopifnot(all(c("x", "t.x", "T.cal") %in% names(cal.cbs)))
  stopifnot(all(covariates %in% names(cal.cbs)))

  # Setup Regressors (Covariates) for location of 1st-stage prior, i.e. beta = [log(lambda), log(mu)]
  cal.cbs[, "intercept"] <- 1
  covariates <- c("intercept", covariates)
  K <- length(covariates)  # number of covars
  covars <- as.matrix(subset(cal.cbs, select = covariates))

  # set hyper priors
  beta_0 <- matrix(0, nrow = K, ncol = 2, dimnames = list(NULL, c("log_lambda", "log_mu")))
  A_0 <- diag(rep(0.01, K), ncol = K, nrow = K)  # diffuse precision matrix
  # set diffuse hyper-parameters for 2nd-stage prior of gamma_0; follows defaults from rmultireg example
  nu_00 <- 3 + K  # 30
  gamma_00 <- nu_00 * diag(2)
  hyper_prior <- list(beta_0 = beta_0, A_0 = A_0, nu_00 = nu_00, gamma_00 = gamma_00)

  # run multiple chains - executed in parallel on Unix
  ncores <- ifelse(!is.null(mc.cores), min(chains, mc.cores), ifelse(.Platform$OS.type == "windows", 1, min(chains,
    detectCores())))
  if (ncores > 1)
    cat("running in parallel on", ncores, "cores\n")
  draws <- mclapply(1:chains, function(i) run_single_chain(i, cal.cbs, hyper_prior = hyper_prior), mc.cores = ncores)

  # merge chains into code::mcmc.list objects
  out <- list(level_1 = lapply(1:nrow(cal.cbs), function(i) mcmc.list(lapply(draws, function(draw) draw$level_1[[i]]))),
    level_2 = mcmc.list(lapply(draws, function(draw) draw$level_2)))
  if ("cust" %in% names(cal.cbs))
    names(out$level_1) <- cal.cbs$cust
  return(out)
}


#' Simulate data according to Pareto/NBD (Abe) model assumptions
#'
#' @param n Number of customers.
#' @param T.cal Length of calibration period. If a vector is provided, then it
#'   is assumed that customers have different 'birth' dates, i.e.
#'   \eqn{max(T.cal)-T.cal}.
#' @param T.star Length of holdout period. This may be a vector.
#' @param params A list of model parameters: \code{beta} and \code{gamma}.
#' @param date.zero Initial date for cohort start. Can be of class character, Date or POSIXt.
#' @return List of length 2:
#' \item{\code{cbs}}{A data.frame with a row for each customer and the summary statistic as columns.}
#' \item{\code{elog}}{A data.frame with a row for each transaction, and columns \code{cust}, \code{date} and \code{t}.}
#' @export
#' @examples
#' # generate artificial Pareto/NBD (Abe) with 2 covariates
#' params <- list()
#' params$beta  <- matrix(c(0.18, -2.5, 0.5, -0.3, -0.2, 0.8), byrow = TRUE, ncol = 2)
#' params$gamma <- matrix(c(0.05, 0.1, 0.1, 0.2), ncol = 2)
#' data <- abe.GenerateData(n = 2000, T.cal = 32, T.star = 32, params)
#' cbs <- data$cbs  # customer by sufficient summary statistic - one row per customer
#' elog <- data$elog  # Event log - one row per event/purchase
abe.GenerateData <- function(n, T.cal, T.star, params, date.zero = "2000-01-01") {

  # set start date for each customer, so that they share same T.cal date
  T.cal.fix <- max(T.cal)
  T.cal <- rep(T.cal, length.out = n)
  T.zero <- T.cal.fix - T.cal
  date.zero <- as.POSIXct(date.zero)

  if (!is.matrix(params$beta))
    params$beta <- matrix(params$beta, nrow = 1, ncol = 2)

  nr_covars <- nrow(params$beta)
  covars <- matrix(c(rep(1, n), runif( (nr_covars - 1) * n, -1, 1)), nrow = n, ncol = nr_covars)
  colnames(covars) <- paste("covariate", 0:(nr_covars - 1), sep = "_")
  colnames(covars)[1] <- "intercept"

  # sample log-normal distributed parameters lambda/mu for each customer
  thetas <- exp( (covars %*% params$beta) + mvtnorm::rmvnorm(n, mean = c(0, 0), sigma = params$gamma))
  lambdas <- thetas[, 1]
  mus <- thetas[, 2]

  # sample lifetime for each customer
  taus <- rexp(n, rate = mus)

  # sample intertransaction timings
  elog_list <- lapply(1:n, function(i) {
    # sample 'sufficiently' large amount of inter-transaction times
    minT <- min(T.cal[i] + max(T.star), taus[i])
    itt_draws <- max(10, round(minT * lambdas[i] * 1.5))
    itt_fn <- function(n) rexp(n, rate = lambdas[i])
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
  set(cbs, j = "lambda", value = lambdas)
  set(cbs, j = "mu", value = mus)
  set(cbs, j = "tau", value = taus)
  set(cbs, j = "alive", value = (T.zero + taus) > T.cal.fix)
  cbs <- cbind(cbs, covars)

  return(list("cbs" = setDF(cbs), "elog" = setDF(elog)))
}
