
#' Hierarchical Bayes variant of Pareto/NBD
#'
#' \code{pnbd.mcmc.DrawParameters} samples parameters via MCMC for a given CBS
#' matrix
#' 
#' method 1) If \code{use_data_augmentation==TRUE} then implementation follows chapter
#' 3.2 of Sandeep Conoor's dissertation
#' \url{http://gradworks.umi.com/34/02/3402149.html}, i.e. parameter space is expanded
#' with unobserved lifetime tau_i. Note, however, that  we follow the notation
#' of original SMC paper with alpha and beta being the 'rate' parameters of the
#' Gamma distributions, instead of 'scale' parameter.
#' 
#' method 2) If \code{use_data_augmentation==FALSE} then implementation follows
#' Shao-Hui Ma & Jin-Lan Liu paper 
#' \url{http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=4344404}, i.e. no data
#' augmentation and draws on individual level need to be done via slice
#' sampling. As such it is 10x slower than method 1)
#' 
#' Estimating parameters via Pareto/NBD MCMC can be 10x slower than Pareto/NBD
#' MLE, which itself can be 10x slower than BG/NBD. Both methods exhibit highly
#' autocorrelated draws of {r, alpha, s, beta} and hence need to be run long, to
#' generate 'enough' draws
#'
#' @param cal.cbs data.frame with columns \code{x}, \code{t.x}, \code{T.cal}
#' @param mcmc number of MCMC steps
#' @param burnin number of initial MCMC steps which are discarded
#' @param thin only every thin-th MCMC step will be returned
#' @param chains number of MCMC chains to be run
#' @param use_data_augmentation determines MCMC method to be used
#' @param param_init list of 2nd-level parameter start values
#' @param trace print logging step every \code{trace} iteration
#' @return 2-element list:
##' \itemize{
##'  \item{\code{level_1}}{list of \code{\link{mcmc.list}} objects; one for each customer, containing individual-level draws}
##'  \item{\code{level_2}}{\code{\link{mcmc.list}} object containing draws of heterogeneity parameters}
##' }
#' @import coda parallel
#' @export
#' @example demo/pareto-cnbd.r
#' @seealso \code{\link{pcnbd.GenerateData}} \code{\link{pcnbd.DrawFutureTransactions}}
#' @references Ma, Shao-Hui, and Jin-Lan Liu. "The MCMC approach for solving the Pareto/NBD model and possible extensions." Natural Computation, 2007. ICNC 2007. Third International Conference on. Vol. 2. IEEE, 2007. \url{http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=4344404}
#' @references Conoor, Sandeep S. Customer-base analysis in noncontractual settings. Diss. NORTHWESTERN UNIVERSITY, 2010. \url{http://gradworks.umi.com/34/02/3402149.html}
pnbd.mcmc.DrawParameters <-
  function(cal.cbs,
           mcmc = 1500, burnin = 500, thin = 50, chains = 2,
           use_data_augmentation = TRUE, param_init = NULL, trace = 100) {
  
  ## methods to sample heterogeneity parameters {r, alpha, s, beta} ##
  
  draw_gamma_params <- function(type, level_1, level_2, hyper_prior) {
    if (type=="lambda") {
      x <- level_1["lambda",]
      cur_params <- c(level_2["r"], level_2["alpha"])
      hyper <- unlist(hyper_prior[c("r_1", "r_2", "alpha_1", "alpha_2")])
    } else if (type=="mu") {
      x <- level_1["mu",]
      cur_params <- c(level_2["s"], level_2["beta"])
      hyper <- unlist(hyper_prior[c("s_1", "s_2", "beta_1", "beta_2")])
    }
    slice_sample_gamma_parameters(x, cur_params, hyper, steps=50, w=0.1)
  }
  
  ## methods to sample individual-level parameters (with data augmentation) ##  
  
  draw_lambda_conoor <- function(data, level_1, level_2) {
    N      <- nrow(data)
    x      <- data[, "x"]
    T.cal  <- data[, "T.cal"]
    tau    <- level_1["tau", ]
    r      <- level_2["r"]
    alpha  <- level_2["alpha"]
    
    lambda <- rgamma(n     = N,
                     shape = r + x,
                     rate  = alpha + pmin(tau, T.cal))
    lambda[lambda==0 | log(lambda) < -70] <- exp(-70) # avoid numeric overflow
    return(lambda)
  }
  
  draw_mu_conoor <- function(data, level_1, level_2) {
    N      <- nrow(data)
    tau    <- level_1["tau", ]  
    s      <- level_2["s"]
    beta   <- level_2["beta"]
    
    mu <- rgamma(n     = N, 
                 shape = s + 1, 
                 rate  = beta + tau)
    mu[mu==0 | log(mu) < -70] <- exp(-70) # avoid numeric overflow
    return(mu)
  }
  
  draw_tau <- function(data, level_1) {
    N      <- nrow(data)
    tx     <- data[, "t.x"]
    T.cal  <- data[, "T.cal"]
    lambda <- level_1["lambda", ]
    mu     <- level_1["mu", ]
    
    mu_lam <- mu + lambda
    t_diff <- T.cal - tx
    
    # sample z
    p_alive <- 1 / (1+(mu/mu_lam)*(exp(mu_lam*t_diff)-1))
    alive   <- p_alive > runif(n=N)
    
    tau <- numeric(N)
    
    # Case: still alive - left truncated exponential distribution -> [T.cal, Inf]
    if (any(alive)) {
      tau[alive]  <- T.cal[alive] + rexp(sum(alive), mu[alive])
    }
    
    # Case: churned     - double truncated exponential distribution -> [tx, T.cal]
    if (any(!alive)) {
      mu_lam_tx   <- pmin(700, mu_lam[!alive] * tx[!alive])
      mu_lam_Tcal <- pmin(700, mu_lam[!alive] * T.cal[!alive])
      # sample with http://en.wikipedia.org/wiki/Inverse_transform_sampling
      rand        <- runif(n=sum(!alive))
      tau[!alive] <- -log( (1-rand)*exp(-mu_lam_tx) + rand*exp(-mu_lam_Tcal)) / mu_lam[!alive]
    }
    
    return(tau)
  }
  
  ## methods to sample individual-level parameters (without data augmentation) ##  
  
  draw_lambda_ma_liu <- function(data, level_1, level_2) {
    slice_sample_ma_liu("lambda", 
                        x = data[,"x"], tx = data[,"t.x"], Tcal = data[,"T.cal"], 
                        lambda = level_1["lambda",], mu = level_1["mu",], 
                        r = level_2["r"], alpha = level_2["alpha"], 
                        s = level_2["s"], beta = level_2["beta"])
  }
  
  draw_mu_ma_liu <- function(data, level_1, level_2) {
    slice_sample_ma_liu("mu", 
                        x = data[,"x"], tx = data[,"t.x"], Tcal = data[,"T.cal"], 
                        lambda = level_1["lambda",], mu = level_1["mu",], 
                        r = level_2["r"], alpha = level_2["alpha"], 
                        s = level_2["s"], beta = level_2["beta"])
  }
  
  run_single_chain <- function(chain_id=1, data) {
    
    ## initialize arrays for storing draws ##
    
    nr_of_cust <- nrow(data)
    nr_of_draws <- (mcmc-1) %/% thin + 1
    level_2_draws <- array(NA_real_, dim=c(nr_of_draws, 4))
    dimnames(level_2_draws)[[2]] <- c("r", "alpha", "s", "beta")
    level_1_draws <- array(NA_real_, dim=c(nr_of_draws, 3, nr_of_cust))
    dimnames(level_1_draws)[[2]] <- c("lambda", "mu", "tau")
    
    ## initialize parameters ##
    
    level_2          <- level_2_draws[1,]
    level_2["r"]     <- param_init$r
    level_2["alpha"] <- param_init$alpha
    level_2["s"]     <- param_init$s
    level_2["beta"]  <- param_init$beta
    
    level_1            <- level_1_draws[1,,]
    level_1["lambda",] <- mean(data$x) / mean(ifelse(data$t.x==0, data$T.cal, data$t.x))
    level_1["tau",]    <- data$t.x + 0.5/level_1["lambda",]
    level_1["mu",]     <- 1/level_1["tau",]  
    
    ## run MCMC chain ##
    
    for (step in 1:(burnin+mcmc)) {
      if (step%%trace==0) cat("chain:", chain_id, "step:", step, "of", (burnin+mcmc), "\n")
      
      # store
      if ((step-burnin)>0 & (step-1-burnin)%%thin==0) {
        idx <- (step-1-burnin)%/%thin + 1
        level_1_draws[idx,,] <- level_1
        level_2_draws[idx,]  <- level_2
      }
      
      # draw individual-level parameters
      draw_lambda <- if (use_data_augmentation) draw_lambda_conoor else draw_lambda_ma_liu
      draw_mu     <- if (use_data_augmentation) draw_mu_conoor else draw_mu_ma_liu
      level_1["lambda", ] <- draw_lambda(data, level_1, level_2)
      level_1["mu", ]     <- draw_mu(data, level_1, level_2)
      level_1["tau", ]    <- draw_tau(data, level_1)
      
      # draw heterogeneity parameters
      level_2[c("r", "alpha")] <- draw_gamma_params("lambda", level_1, level_2, hyper_prior)
      level_2[c("s", "beta")]  <- draw_gamma_params("mu", level_1, level_2, hyper_prior)
    }
    
    # convert MCMC draws into coda::mcmc objects
    return(list(
      level_1 = lapply(1:nr_of_cust, function(i) mcmc(level_1_draws[,,i], start=burnin, thin=thin)),
      level_2 = mcmc(level_2_draws, start=burnin, thin=thin)))
  }
  
  # set hyper priors
  hyper_prior <- list(r_1 = 1e-3, r_2 = 1e-3,
                      alpha_1 = 1e-3, alpha_2 = 1e-3,
                      s_1 = 1e-3, s_2 = 1e-3,
                      beta_1 = 1e-3, beta_2 = 1e-3)
  
  # set param_init (if not passed as argument)
  if (is.null(param_init)) {
    tryCatch({
      df <- cal.cbs[sample(nrow(cal.cbs), min(nrow(cal.cbs), 1000)),]
      param_init <- pnbd.EstimateParameters(df)
      names(param_init) <- c("r", "alpha", "s", "beta")
      param_init <- as.list(param_init)
    }, error = function(e) param_init <- list(r=1, alpha=1, s=1, beta=1))
    cat("set param_init:", paste(round(unlist(param_init), 4), collapse=", "), "\n")
  }
  
  # check whether input data meets requirements
  stopifnot(is.data.frame(cal.cbs))
  stopifnot(all(c("x", "t.x", "T.cal") %in% names(cal.cbs)))
  stopifnot(all(is.finite(cal.cbs$litt)))
  
  # run multiple chains - executed in parallel on Unix
  ncores <- ifelse(.Platform$OS.type=="windows", 1, min(chains, detectCores()))
  if (ncores>1) cat("running in parallel on", ncores, "cores\n")
  draws <- mclapply(1:chains, function(i) run_single_chain(i, cal.cbs), mc.cores=ncores)

  # merge chains into code::mcmc.list objects
  return(list(
    level_1 = lapply(1:nrow(cal.cbs), function(i) mcmc.list(lapply(draws, function(draw) draw$level_1[[i]]))),
    level_2 = mcmc.list(lapply(draws, function(draw) draw$level_2))))
}
