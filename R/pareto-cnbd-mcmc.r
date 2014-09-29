
#' Hierarchical Bayes implementation of Pareto/CNBD
#'
#' Returns 2-element list
#'   level_1:  3-dim array [draw x parameter x cust] wrapped as coda::mcmc.list object
#'   level_2:  2-dim array [draw x parameter] wrapped as coda::mcmc.list object
#'
#' @param cal.cbs data.frame with columns \code{x}, \code{t.x}, \code{T.cal}, \code{litt}; e.g. output of \code{\link{elog2cbs}}
#' @param mcmc number of MCMC steps
#' @param burnin number of initial MCMC steps which are discarded
#' @param thin only every thin-th MCMC step will be returned
#' @param chains number of MCMC chains to be run
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
#' @seealso \code{link{pcnbd.GenerateData}} \code{\link{mcmc.PAlive}} \code{\link{mcmc.DrawFutureTransactions}} \code{\link{elog2cbs}}
pcnbd.mcmc.DrawParameters <-
  function(cal.cbs,
           mcmc = 1500, burnin = 500, thin = 50, chains = 2,
           param_init = NULL, trace = 100) {
  
  ## methods to sample heterogeneity parameters {r, alpha, s, beta, t, gamma} ##
    
  draw_gamma_params <- function(type, level_1, level_2, hyper_prior) {
    if (type=="lambda") {
      x <- level_1["lambda",]
      cur_params <- c(level_2["r"], level_2["alpha"])
      hyper <- unlist(hyper_prior[c("r_1", "r_2", "alpha_1", "alpha_2")])
    } else if (type=="mu") {
      x <- level_1["mu",]
      cur_params <- c(level_2["s"], level_2["beta"])
      hyper <- unlist(hyper_prior[c("s_1", "s_2", "beta_1", "beta_2")])
    } else if (type=="k") {
      x <- level_1["k",]
      cur_params <- c(level_2["t"], level_2["gamma"])
      hyper <- unlist(hyper_prior[c("t_1", "t_2", "gamma_1", "gamma_2")])      
    }
    slice_sample_gamma_parameters(x, cur_params, hyper, steps=200, w=0.1)
  }
  
  ## methods to sample individual-level parameters ##
  
  draw_k <- function(data, level_1, level_2) {
    pcnbd_slice_sample("k",
                       x = data$x, tx = data$t.x, Tcal = data$T.cal, litt = data$litt, 
                       k = level_1["k",], lambda = level_1["lambda",], mu = level_1["mu",], tau = level_1["tau",],
                       t = level_2["t"], gamma = level_2["gamma"],
                       r = level_2["r"], alpha = level_2["alpha"], 
                       s = level_2["s"], beta  = level_2["beta"])
  }
  
  draw_lambda <- function(data, level_1, level_2) {
    pcnbd_slice_sample("lambda",
                       x = data$x, tx = data$t.x, Tcal = data$T.cal, litt = data$litt, 
                       k = level_1["k",], lambda = level_1["lambda",], mu = level_1["mu",], tau = level_1["tau",],
                       t = level_2["t"], gamma = level_2["gamma"],
                       r = level_2["r"], alpha = level_2["alpha"], 
                       s = level_2["s"], beta  = level_2["beta"])
  }
  
  draw_mu <- function(data, level_1, level_2) {
    N      <- nrow(data)
    tau    <- level_1["tau", ]  
    s      <- level_2["s"]
    beta   <- level_2["beta"]
    
    mu <- rgamma(n     = N, 
                 shape = s + 1, 
                 rate  = beta + tau)
    mu[mu==0 | log(mu) < -10] <- exp(-10) # avoid numeric overflow
    return(mu)
  }
  
  draw_tau <- function(data, level_1, level_2) {
    N      <- nrow(data)
    x      <- data$x
    tx     <- data$t.x
    Tcal   <- data$T.cal
    lambda <- level_1["lambda", ]
    k      <- level_1["k", ]
    mu     <- level_1["mu", ]
    
    # sample z
    p_alive <- pcnbd_palive(x, tx, Tcal, k, lambda, mu)
    alive   <- p_alive > runif(n=N)
    
    # sample tau
    tau <- numeric(N)
    
    # Case: still alive - left truncated exponential distribution -> [Tcal, Inf]
    if (any(alive)) {
      tau[alive]  <- Tcal[alive] + rexp(sum(alive), mu[alive])
    }
    
    # Case: churned     - distribution of min(t_(x+1), tau), truncated to [tx, Tcal]
    if (any(!alive)) {
      tau[!alive] <- pcnbd_slice_sample("tau",
                                        x = data$x[!alive], tx = data$t.x[!alive], Tcal = data$T.cal[!alive], litt = data$litt[!alive],
                                        k = level_1["k",!alive], lambda = level_1["lambda",!alive], mu = level_1["mu",!alive], tau = level_1["tau",!alive],
                                        t = level_2["t"], gamma = level_2["gamma"],
                                        r = level_2["r"], alpha = level_2["alpha"], 
                                        s = level_2["s"], beta  = level_2["beta"])
    }
    
    return(tau)
  }
  
  run_single_chain <- function(chain_id, data) {
    
    ## initialize arrays for storing draws ##
    
    nr_of_cust <- nrow(data)
    nr_of_draws <- (mcmc-1) %/% thin + 1
    level_2_draws <- array(NA_real_, dim=c(nr_of_draws, 6))
    dimnames(level_2_draws)[[2]] <- c("t", "gamma", "r", "alpha", "s", "beta")
    level_1_draws <- array(NA_real_, dim=c(nr_of_draws, 5, nr_of_cust))
    dimnames(level_1_draws)[[2]] <- c("k", "lambda", "mu", "tau", "z")
    
    ## initialize parameters ##
    
    level_2          <- level_2_draws[1,]
    level_2["t"]     <- param_init$t
    level_2["gamma"] <- param_init$gamma
    level_2["r"]     <- param_init$r
    level_2["alpha"] <- param_init$alpha
    level_2["s"]     <- param_init$s
    level_2["beta"]  <- param_init$beta
    
    level_1            <- level_1_draws[1,,]
    level_1["k",]      <- 1
    level_1["lambda",] <- mean(data$x) / mean(ifelse(data$t.x==0, data$T.cal, data$t.x))
    level_1["tau",]    <- data$t.x + 0.5/level_1["lambda",]
    level_1["z",]      <- as.numeric(level_1["tau",] > data$T.cal)
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
      level_1["k", ]      <- draw_k(data, level_1, level_2)
      level_1["lambda", ] <- draw_lambda(data, level_1, level_2)
      level_1["mu", ]     <- draw_mu(data, level_1, level_2)
      level_1["tau", ]    <- draw_tau(data, level_1, level_2)
      level_1["z", ]      <- as.numeric(level_1["tau",] > data$T.cal)
      
      # draw heterogeneity parameters
      level_2[c("t", "gamma")] <- draw_gamma_params("k", level_1, level_2, hyper_prior)
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
                      beta_1 = 1e-3, beta_2 = 1e-3,
                      t_1 = 1e-3, t_2 = 1e-3,
                      gamma_1 = 1e-3, gamma_2 = 1e-3)  
  
  # set param_init (if not passed as argument)
  if (is.null(param_init)) {
    try({
      df <- cal.cbs[sample(nrow(cal.cbs), min(nrow(cal.cbs), 1000)),]
      param_init <- c(1, 1, BTYD::pnbd.EstimateParameters(df))
      names(param_init) <- c("t", "gamma", "r", "alpha", "s", "beta")
      param_init <- as.list(param_init)
    }, silent=TRUE)
    if (is.null(param_init)) param_init <- list(t=1, gamma=1, r=1, alpha=1, s=1, beta=1)
    cat("set param_init:", paste(round(unlist(param_init), 4), collapse=", "), "\n")
  }
  
  # check whether input data meets requirements
  stopifnot(is.data.frame(cal.cbs))
  stopifnot(all(c("x", "t.x", "T.cal", "litt") %in% names(cal.cbs)))
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


#' Generate artificial data which follows Pareto/CNBD model assumptions
#'
#' Returns 2-element list
#' * cbs: data.frame with 'cust', \code{x}, \code{t.x}, \code{T.cal}, 'T.star', 'x.star' 
#'        this is the summary statistics data.frame which contains all 
#'        needed information for parameter estimation
#' * elog: data.frame with 'cust', \code{t}
#'
#' @param n number of customers
#' @param T.cal length of calibration period
#' @param T.star length of holdout period
#' @param params list of parameters: {r, alpha, s, beta, t, gamma}
#' @param return.elog if \code{TRUE} then the event-log is returned as well; decreases performance
#' 
#' @return 2-elemnt list
#' @export
pcnbd.GenerateData <- function(n, T.cal, T.star, params, return.elog=FALSE) {
  
  if (length(T.cal)==1) T.cal <- rep(T.cal, n)
  if (length(T.star)==1) T.star <- rep(T.star, n)
  
  # sample regularity parameter k for each customer
  if (all(c("t", "gamma") %in% names(params))) {
    # Case A: regularity parameter k is gamma-distributed across customers
    ks <- rgamma(n, shape=params$t, rate=params$gamma)
    ks <- pmax(0.01, ks) # ensure that k is not too small, otherwise itt can be 0
    
  } else if ("k" %in% params) {
    # Case B: regularity parameter k is fixed across customers
    ks <- rep(param$k, n)
    
  } else {
    # Case C: k=1 is assumed, i.e. Pareto/NBD
    ks <- rep(1, n)
  }
  
  # sample intertransaction timings parameter lambda for each customer
  lambdas <- rgamma(n, shape=params$r, rate=params$alpha)
  
  # sample lifetime for each customer
  mus <- rgamma(n, shape=params$s, rate=params$beta)
  taus <- rexp(n, rate=mus)
  
  # sample intertransaction timings & churn
  cbs_list <- list()
  elog_list <- list()
  for (i in 1:n) {
    k      <- ks[i]
    lambda <- lambdas[i]
    mu     <- mus[i]
    tau    <- taus[i]
    # sample 'sufficiently' large amount of inter-transaction times
    minT <- min(T.cal[i]+T.star[i], tau)
    nr_of_itt_draws <- max(10, round(minT * lambda))
    itts <- rgamma(nr_of_itt_draws * 2, shape=k, rate=k*lambda)
    if (sum(itts)<minT) itts <- c(itts, rgamma(nr_of_itt_draws * 4, shape=k, rate=k*lambda))
    if (sum(itts)<minT) itts <- c(itts, rgamma(nr_of_itt_draws * 800, shape=k, rate=k*lambda))
    if (sum(itts)<minT)
      stop("not enough inter-transaction times sampled: ", sum(itts), " < ", minT)
    times <- cumsum(c(0, itts))
    times <- times[times<tau]
    if (return.elog)
      elog_list[[i]] <- data.frame(cust=i, t=times[times<(T.cal[i]+T.star[i])])
    # determine frequency, recency, etc.
    ts.cal   <- times[times<T.cal[i]]
    ts.star  <- times[times>=T.cal[i] & times<(T.cal[i]+T.star[i])]
    cbs_list[[i]] <- list(cust   = i,
                          x      = length(ts.cal)-1,
                          t.x    = max(ts.cal),
                          litt   = ifelse(length(ts.cal)-1==0, 0, sum(log(itts[1:(length(ts.cal)-1)]))),
                          alive  = tau>T.cal[i],
                          x.star = length(ts.star))
  }
  cbs <- do.call(rbind.data.frame, cbs_list)
  cbs$k      <- ks
  cbs$lambda <- lambdas
  cbs$mu     <- mus
  cbs$tau    <- taus
  cbs$T.cal  <- T.cal
  cbs$T.star <- T.star
  rownames(cbs) <- NULL
  out <- list(cbs=cbs)
  if (return.elog) {
    elog <- do.call(rbind.data.frame, elog_list)
    out$elog <- elog
  }
  return(out)
}


#' Pareto/CNBD Plot Regularity Rate Heterogeneity
#' 
#' Plots and returns the estimated gamma distribution of k (customers' regularity in interpurchase times).
#' 
#' @param draws MCMC draws returned by \code{\link{pcnbd.mcmc.DrawParameters}}
#' @param xmax upper bound for x-scale
#'
#' @export
pcnbd.mcmc.plotRegularityRateHeterogeneity <- function(draws, xmax=NULL) {
  ks <- sapply(draws$level_1, function(draw) as.matrix(draw[, "k"]))
  if (is.null(xmax)) xmax <- min(10, quantile(ks, 0.95)*1.5)
  plot(density(ks), xlim=c(0, xmax), main="Distribution of Regularity Rate k", xlab="", ylab="", frame=FALSE)
  #plot(density(apply(ks, 2, mean)), xlim=c(0, xmax), main="Distribution of Regularity Rate k", xlab="", ylab="", frame=FALSE)
  abline(v=1, lty=3)
  abline(v=median(ks), col="red")
}
