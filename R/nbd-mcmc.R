
#' NBD (HB) Parameter Draws
#'
#' Returns draws from the posterior distributions of the NBD (HB)
#' parameters, on cohort as well as on customer level.
#'
#' @param cal.cbs Calibration period customer-by-sufficient-statistic (CBS)
#'   data.frame. It must contain a row for each customer, and columns \code{x}
#'   for frequency,  and \code{T.cal} for the total time observed. A correct
#'   format can be easily generated based on the complete event log of a
#'   customer cohort with \code{\link{elog2cbs}}.
#' @param mcmc Number of MCMC steps.
#' @param burnin Number of initial MCMC steps which are discarded.
#' @param thin Only every \code{thin}-th MCMC step will be returned.
#' @param chains Number of MCMC chains to be run.
#' @param mc.cores Number of cores to use in parallel (Unix only). Defaults to \code{min(chains, detectCores())}.
#' @param param_init List of start values for cohort-level parameters.
#' @param trace Print logging statement every \code{trace}-th iteration. Not available for \code{mc.cores > 1}.
#' @return 2-element list:
#' \itemize{
#'  \item{\code{level_1 }}{list of \code{\link{mcmc.list}}s, one for each customer, with draws for customer-level parameters \code{lambda}}
#'  \item{\code{level_2 }}{\code{\link{mcmc.list}}, with draws for cohort-level parameters \code{r}, \code{alpha}}
#' }
#' @export
#' @seealso \code{\link{nbd.GenerateData} }
nbd.mcmc.DrawParameters <- function(cal.cbs, mcmc = 2500, burnin = 500, thin = 50, chains = 2, mc.cores = NULL,
  param_init = NULL, trace = 100) {

  # ** methods to sample heterogeneity parameters {r, alpha} **

  draw_gamma_params <- function(type, level_1, level_2, hyper_prior) {
    x <- level_1["lambda", ]
    cur_params <- c(level_2["r"], level_2["alpha"])
    hyper <- unlist(hyper_prior[c("r_1", "r_2", "alpha_1", "alpha_2")])
    slice_sample_gamma_parameters(x, cur_params, hyper, steps = 50, w = 0.1)
  }

  draw_lambda <- function(data, level_1, level_2) {
    N <- nrow(data)
    x <- data$x
    T.cal <- data$T.cal
    r <- level_2["r"]
    alpha <- level_2["alpha"]

    lambda <- rgamma(n = N, shape = r + x, rate = alpha + T.cal)
    lambda[lambda == 0 | log(lambda) < -30] <- exp(-30)  # avoid numeric overflow
    return(lambda)
  }

  run_single_chain <- function(chain_id = 1, data, hyper_prior) {

    ## initialize arrays for storing draws ##

    nr_of_cust <- nrow(data)
    nr_of_draws <- (mcmc - 1) %/% thin + 1
    level_2_draws <- array(NA_real_, dim = c(nr_of_draws, 2))
    dimnames(level_2_draws) <- list(NULL, c("r", "alpha"))
    level_1_draws <- array(NA_real_, dim = c(nr_of_draws, 1, nr_of_cust))
    dimnames(level_1_draws) <- list(NULL, c('lambda'), NULL)

    ## initialize parameters ##

    level_2 <- level_2_draws[1, ]
    level_2["r"] <- param_init$r
    level_2["alpha"] <- param_init$alpha

    level_1 <- array(level_1_draws[1, , ], dim = dim(level_1_draws)[-1], dimnames = list("lambda", NULL))
    level_1["lambda", ] <- mean(data$x) / mean(data$T.cal)

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
      level_1["lambda", ] <- draw_lambda(data, level_1, level_2)

      # draw heterogeneity parameters
      level_2[c("r", "alpha")] <- draw_gamma_params("lambda", level_1, level_2, hyper_prior)
    }

    # convert MCMC draws into coda::mcmc objects
    return(list(
      "level_1" = lapply(1:nr_of_cust, function(i) {
                           level_1_draws_i <- array(level_1_draws[, , i], dim = dim(level_1_draws)[-3], dimnames = list(NULL, "lambda"))
                           mcmc(level_1_draws_i, start = burnin, thin = thin)
                         }), # nolint
      "level_2" = mcmc(level_2_draws, start = burnin, thin = thin)))
  }

  # set hyper priors
  hyper_prior <- list(r_1 = 0.001, r_2 = 0.001,
                      alpha_1 = 0.001, alpha_2 = 0.001)

  # set param_init (if not passed as argument)
  if (is.null(param_init)) {
    try({
        df <- cal.cbs[sample(nrow(cal.cbs), min(nrow(cal.cbs), 1000)), ]
        param_init <- nbd.EstimateParameters(df)
        names(param_init) <- c("r", "alpha")
        param_init <- as.list(param_init)
      },
      silent = TRUE)
    if (is.null(param_init))
      param_init <- list(r = 1, alpha = 1)
    cat("set param_init:", paste(round(unlist(param_init), 4), collapse = ", "), "\n")
  }

  # check whether input data meets requirements
  stopifnot(is.data.frame(cal.cbs))
  stopifnot(all(c("x", "T.cal") %in% names(cal.cbs)))

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
