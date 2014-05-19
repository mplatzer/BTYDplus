
#' Hierarchical Bayes implementation of Pareto/CNBD
#'
#' Returns 2-element list
#'   level_1:  3-dim array [draw x parameter x cust] wrapped as coda::mcmc.list object
#'   level_2:  2-dim array [draw x parameter] wrapped as coda::mcmc.list object
#'
#' @param data data.frame with columns 'x', 't.x', 'T.cal'
#' @param mcmc number of MCMC steps
#' @param burnin number of initial MCMC steps which are discarded
#' @param thin only every thin-th MCMC step will be returned
#' @param chains number of MCMC chains to be run
#' @param use_data_augmentation determines MCMC method to be used
#' @param param_init list of 2nd-level parameter start values
#' @param hyper_prior list of hyper parameters for 2nd-level parameters
#' @return list
#' @import coda Rcpp parallel
#' @export
#' @examples
#' #params <- list(r=1.4, alpha=1.3, s=0.7, beta=7, t=1.5, gamma=1)
#' #cbs <- pcnbd.GenerateData(1000, 10, 5, params)$cbs
#' #draws <- pcnbd.mcmc.DrawParameters(cbs, mcmc=1000, chains=2)
#' #plot(draws$level_2, ask=F)
#' #rbind("actual"=unlist(params), "estimated"=summary(draws$level_2, quantiles=0.5)$quantiles)
#' @seealso pcnbd.GenerateData
pcnbd.mcmc.DrawParameters <-
  function(data,
           mcmc = 10000, burnin = 0, thin = 1, chains = 1,           
           param_init = list(t=1, gamma=1, r=1, alpha=1, s=1, beta=1),
           mc.cores = 1,
           hyper_prior = list(t_1=1/1000, t_2=1/1000,
                              gamma_1=1/1000, gamma_2=1/1000,
                              r_1=1/1000, r_2=1/1000,
                              alpha_1=1/1000, alpha_2=1/1000,
                              s_1=1/1000, s_2=1/1000,
                              beta_1=1/1000, beta_2=1/1000)) {
    
  ## methods to sample heterogeneity parameters {r, alpha, s, beta, t, gamma} ##
  
  draw_shape <- function(x, shape, rate, prior1, prior2, steps=20) {
    
    calc_prior <- function(r) (prior1 - 1) * log(r) - (r * prior2)
    calc_likel <- function(r) (r - 1) * sum(log(x)) + length(x) * (r * log(rate) - lgamma(r))
    
    cur_shape <- shape
    cur_prior <- calc_prior(cur_shape)
    cur_likel <- calc_likel(cur_shape)
    
    scale <- mean(x)
    
    for (i in 1:steps) {
      
      # generate new proposal
      new_shape <- cur_shape * exp(max(-100, min(100, scale * rt(1, df=3))))
      new_prior <- calc_prior(new_shape)
      new_likel <- calc_likel(new_shape)
      
      # accept/reject new proposal
      mhratio <- exp(new_prior+new_likel-cur_prior-cur_likel)
      if (mhratio > runif(n=1))
        cur_shape <- new_shape
      
    }
    return(cur_shape)
  }
  
  draw_rate <- function(x, shape, rate, prior1, prior2) {
    rgamma(n     = 1, 
           shape = prior1 + length(x) * shape,
           rate  = prior2 + sum(x))
  }
  
  draw_t <- function(level_1, level_2, hyper_prior) {
    draw_shape(x      = level_1["k",], 
               shape  = level_2["t"],
               rate   = level_2["gamma"],
               prior1 = hyper_prior$t_1,
               prior2 = hyper_prior$t_2)
  }
  
  draw_gamma <- function(level_1, level_2, hyper_prior) {
    draw_rate(x      = level_1["k",],
              shape  = level_2["t"],
              rate   = level_2["gamma"],
              prior1 = hyper_prior$gamma_1,
              prior2 = hyper_prior$gamma_2)
  }
  
  draw_r <- function(level_1, level_2, hyper_prior) {
    draw_shape(x      = level_1["lambda",], 
               shape  = level_2["r"],
               rate   = level_2["alpha"],
               prior1 = hyper_prior$r_1,
               prior2 = hyper_prior$r_2)
  }
  
  draw_alpha <- function(level_1, level_2, hyper_prior) {
    draw_rate(x      = level_1["lambda",],
              shape  = level_2["r"],
              rate   = level_2["alpha"],
              prior1 = hyper_prior$alpha_1,
              prior2 = hyper_prior$alpha_2)
  }
  
  draw_s <- function(level_1, level_2, hyper_prior) {
    draw_shape(x      = level_1["mu",], 
               shape  = level_2["s"],
               rate   = level_2["beta"],
               prior1 = hyper_prior$s_1,
               prior2 = hyper_prior$s_2)
  }
  
  draw_beta <- function(level_1, level_2, hyper_prior) {
    draw_rate(x      = level_1["mu",],
              shape  = level_2["s"],
              rate   = level_2["beta"],
              prior1 = hyper_prior$beta_1,
              prior2 = hyper_prior$beta_2)
  }
    
  ## methods to sample individual-level parameters ##
  
  draw_k <- function(data, level_1, level_2) {
    pcnbd_slice_sample("k",
                       x = data[,"x"], tx = data[,"t.x"], Tcal = data[,"T.cal"], litt = data[,"litt"], 
                       k = level_1["k",], lambda = level_1["lambda",], mu = level_1["mu",], tau = level_1["tau",],
                       t = level_2["t"], gamma = level_2["gamma"],
                       r = level_2["r"], alpha = level_2["alpha"], 
                       s = level_2["s"], beta  = level_2["beta"])
  }
  
  draw_lambda <- function(data, level_1, level_2) {
    pcnbd_slice_sample("lambda",
                       x = data[,"x"], tx = data[,"t.x"], Tcal = data[,"T.cal"], litt = data[,"litt"], 
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
    mu[mu==0 | log(mu) < -70] <- exp(-70) # avoid numeric overflow
    return(mu)
  }
  
  draw_tau <- function(data, level_1, level_2) {
    N      <- nrow(data)
    x      <- data[, "x"]
    tx     <- data[, "t.x"]
    Tcal   <- data[, "T.cal"]
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
                         x = data[!alive,"x"], tx = data[!alive,"t.x"], Tcal = data[!alive,"T.cal"], litt = data[!alive,"litt"],
                         k = level_1["k",!alive], lambda = level_1["lambda",!alive], mu = level_1["mu",!alive], tau = level_1["tau",!alive],
                         t = level_2["t"], gamma = level_2["gamma"],
                         r = level_2["r"], alpha = level_2["alpha"], 
                         s = level_2["s"], beta  = level_2["beta"])
    }
    
    return(tau)
  }
  
  run_single_chain <- function(chain_id, df) {
    
    ## initialize arrays for storing draws ##
    
    nr_of_cust <- nrow(df)
    nr_of_draws <- (mcmc-1) %/% thin + 1
    level_2_draws <- array(NA_real_, dim=c(nr_of_draws, 6))
    dimnames(level_2_draws)[[2]] <- c("t", "gamma", "r", "alpha", "s", "beta")
    level_1_draws <- array(NA_real_, dim=c(nr_of_draws, 4, nr_of_cust))
    dimnames(level_1_draws)[[2]] <- c("k", "lambda", "mu", "tau")
    
    # scale time-dimension so that T~1
    
#     scale <- 1 / mean(df$T.cal)
#     cat("scale1", scale, "\n")
#     df <- transform(df,
#                     t.x = t.x * scale,
#                     T.cal = T.cal * scale,
#                     T.star = T.star * scale,
#                     litt = litt + x * log(scale))
    
    ## initialize parameters ##
    
    level_2          <- level_2_draws[1,]
    level_2["t"]     <- param_init$t
    level_2["gamma"] <- param_init$gamma
    level_2["r"]     <- param_init$r
    level_2["alpha"] <- param_init$alpha # FIXME
    level_2["s"]     <- param_init$s
    level_2["beta"]  <- param_init$beta
    
    level_1            <- level_1_draws[1,,]
    level_1["k",]      <- 1
    level_1["lambda",] <- mean(df$x) / mean(ifelse(df$t.x==0, df$T.cal, df$t.x))
    level_1["tau",]    <- df$t.x + 0.5/level_1["lambda",]
    level_1["mu",]     <- 1/level_1["tau",]
    
    ## run MCMC chain ##
    
    for (step in 1:(burnin+mcmc)) {
      if (step%%100==0) cat("chain:", chain_id, "step:", step, "of", (burnin+mcmc), "\n")      
      
      # store
      if ((step-burnin)>0 & (step-1-burnin)%%thin==0) {
        idx <- (step-1-burnin)%/%thin + 1
        level_1_draws[idx,,] <- level_1
        level_2_draws[idx,]  <- level_2
      }
      
      # draw individual-level parameters
      level_1["k", ]      <- draw_k(df, level_1, level_2)
      level_1["lambda", ] <- draw_lambda(df, level_1, level_2)
      level_1["mu", ]     <- draw_mu(df, level_1, level_2)
      level_1["tau", ]    <- draw_tau(df, level_1, level_2)
      
      # draw heterogeneity parameters
      level_2["t"]        <- draw_t(level_1, level_2, hyper_prior)
      level_2["gamma"]    <- draw_gamma(level_1, level_2, hyper_prior)
      level_2["r"]        <- draw_r(level_1, level_2, hyper_prior)
      level_2["alpha"]    <- draw_alpha(level_1, level_2, hyper_prior)
      level_2["s"]        <- draw_s(level_1, level_2, hyper_prior)
      level_2["beta"]     <- draw_beta(level_1, level_2, hyper_prior)
    }
    
    # rescale time-dimension back 
#     cat("scale2", scale, "\n")
#     level_1_draws[,"lambda",] <- level_1_draws[,"lambda",] * scale
#     level_1_draws[,"mu",]     <- level_1_draws[,"mu",] * scale
#     level_1_draws[,"tau",]    <- level_1_draws[,"tau",] / scale
#     level_2_draws[,"alpha"]   <- level_2_draws[,"alpha"] / scale
#     level_2_draws[,"beta"]    <- level_2_draws[,"beta"] / scale
    
    # convert MCMC draws into coda::mcmc objects
    return(list(
      level_1 = lapply(1:nr_of_cust, function(i) mcmc(level_1_draws[,,i], start=burnin, thin=thin)),
      level_2 = mcmc(level_2_draws, start=burnin, thin=thin)))
  }

  ## check whether input data meets requirements
  stopifnot(is.data.frame(data))
  stopifnot(all(c("x", "t.x", "T.cal", "litt") %in% names(data)))
  stopifnot(all(is.finite(data$litt)))
  
  # run multiple chains
  draws <- mclapply(1:chains, function(i) run_single_chain(i, data), mc.cores=mc.cores)

  # merge chains into code::mcmc.list objects
  return(list(
    level_1 = lapply(1:1:nrow(data), function(i) mcmc.list(lapply(draws, function(draw) draw$level_1[[i]]))),
    level_2 = mcmc.list(lapply(draws, function(draw) draw$level_2))))
}


#' Calculates P(alive) based on MCMC draws
#'
#' @param data data.frame with column 'T.cal'
#' @param draws MCMC draws returned by \code{pcnbd.mcmc.DrawParameters}
#' @return numeric vector of probabilities
#' @export
pcnbd.mcmc.PAlive <- function(data, draws) {

  nr_of_draws <- niter(draws$level_2) * nchain(draws$level_2)
  nr_of_cust <- length(draws$level_1)
  parameters <- nvar(draws$level_2[[1]])
  if (nr_of_cust != nrow(data))
    stop("mismatch between number of customers in parameters 'data' and 'draws'")
  
  p.alives <- sapply(1:nr_of_cust, function(i) mean(as.matrix(draws$level_1[[i]][, "tau"]) > data$T.cal[i]))
  return(p.alives)
}


#' Samples number of future transactions based on drawn parameters
#'
#' @param data data.frame with column 't.x' and 'T.cal'
#' @param draws MCMC draws returned by \code{pcnbd.mcmc.DrawParameters}
#' @param T.star length of period for which future transactions are counted
#' @return 2-dim array [draw x cust] with sampled future transactions
#' @export
pcnbd.mcmc.DrawFutureTransactions <- function(data, draws, T.star=data$T.star) {

  nr_of_draws <- niter(draws$level_2) * nchain(draws$level_2)
  nr_of_cust <- length(draws$level_1)
  parameters <- varnames(draws$level_1[[1]])
  
  if (nr_of_cust != nrow(data))
    stop("mismatch between number of customers in parameters 'data' and 'draws'")
  if (is.null(T.star))
    stop("T.star is missing")

  x.stars <- array(NA_real_, dim=c(nr_of_draws, nr_of_cust))
  if (length(T.star)==1) T.star <- rep(T.star, nr_of_cust)
  
  draw_left_truncated_gamma <- function(lower, k, lambda) {
    rand <- runif(1, pgamma(lower, k, k*lambda), 1)
    qgamma(rand, k, k*lambda)
  }
  
  for (cust in 1:nrow(data)) {
    Tcal    <- data[cust, "T.cal"]
    Tstar   <- T.star[cust]
    tx      <- data[cust, "t.x"]
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
      if (sum(itts)<minT) itts <- c(itts, rgamma(nr_of_itt_draws * 8, shape=ks[draw], rate=ks[draw]*lambdas[draw]))
      if (sum(itts)<minT)
        stop("not enough inter-transaction times sampled: ", sum(itts), " < ", minT)
      x.stars[draw, cust] <- sum(cumsum(itts)<minT)
    }
    
    # Case: customer churned
    if (any(!alive)) {
      x.stars[!alive, cust] <- 0
    }
  }
  return(x.stars)
}


#' Generate artificial data which follows Pareto/CNBD model assumptions
#'
#' Returns 2-element list
#' * cbs: data.frame with 'cust', 'x', 't.x', 'T.cal', 'T.star', 'x.star' 
#'        this is the summary statistics data.frame which contains all 
#'        needed information for parameter estimation
#' * elog: data.frame with 'cust', 't'
#'
#' @param N number of customers
#' @param T.cal length of calibration period
#' @param T.star length of holdout period
#' @param params list of parameters: {r, alpha, s, beta, t, gamma}
#' @param return.elog if TRUE then the event-log is returned as well; decreases performance
#' 
#' @return 2-elemnt list
#' @export
pcnbd.GenerateData <- function(N, T.cal, T.star, params, return.elog=F) {

  if (length(T.star)==1) T.star <- rep(T.star, N)
  if (length(T.cal)==1) T.cal <- rep(T.cal, N)

  # sample regularity parameter k for each customer
  if (all(c("t", "gamma") %in% names(params))) {
    # Case A: regularity parameter k is gamma-distributed across customers
    ks <- rgamma(N, shape=params$t, rate=params$gamma)
    ks <- pmax(0.01, ks) # ensure that k is not too small, otherwise itt can be 0
    
  } else if ("k" %in% params) {
    # Case B: regularity parameter k is fixed across customers
    ks <- rep(param$k, N)
    
  } else {
    # Case C: k=1 is assumed, i.e. Pareto/NBD
    ks <- rep(1, N)
  }
  
  # sample intertransaction timings parameter lambda for each customer
  lambdas <- rgamma(N, shape=params$r, rate=params$alpha)
  
  # sample lifetime for each customer
  mus <- rgamma(N, shape=params$s, rate=params$beta)
  taus <- rexp(N, rate=mus)
  
  # sample intertransaction timings & churn
  cbs_list <- list()
  elog_list <- list()
  for (i in 1:N) {
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
  row.names(data) <- NULL
  out <- list(cbs=cbs)
  if (return.elog) {
    elog <- do.call(rbind.data.frame, elog_list)
    out$elog <- elog
  }
  return(out)
}
