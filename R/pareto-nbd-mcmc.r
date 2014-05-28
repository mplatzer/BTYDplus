
#' Hierarchical Bayes variant of Pareto/NBD
#'
#' \code{pnbd.mcmc.DrawParameters} samples parameters via MCMC for a given CBS
#' matrix
#' 
#' method 1) If \code{use_data_augmentation==TRUE} then implementation follows chapter
#' 3.2 of Sandeep Conoor's dissertation
#' http://gradworks.umi.com/34/02/3402149.html, i.e. parameter space is expanded
#' with unobserved lifetime tau_i. Note, however, that  we follow the notation
#' of original SMC paper with alpha and beta being the 'rate' parameters of the
#' Gamma distributions, instead of 'scale' parameter.
#' 
#' method 2) If \code{use_data_augmentation==FALSE} then implementation follows
#' Shao-Hui Ma & Jin-Lan Liu paper 
#' http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=4344404, i.e. no data
#' augmentation and draws on individual level need to be done via slice
#' sampling. As such it is 10x slower than method 1)
#' 
#' Estimating parameters via Pareto/NBD MCMC can be 10x slower than Pareto/NBD
#' MLE, which itself can be 10x slower than BG/NBD. Both methods exhibit highly
#' autocorrelated draws of {r, alpha, s, beta} and hence need to be run long, to
#' generate 'enough' draws
#'
#' @param data data.frame with columns 'x', 't.x', 'T.cal'
#' @param mcmc number of MCMC steps
#' @param burnin number of initial MCMC steps which are discarded
#' @param thin only every thin-th MCMC step will be returned
#' @param chains number of MCMC chains to be run
#' @param use_data_augmentation determines MCMC method to be used
#' @param param_init list of 2nd-level parameter start values
#' @param hyper_prior list of hyper parameters for 2nd-level parameters
#' @return 2-element list
#' level_1:  list of coda::mcmc.list objects; one for each customer, containing individual-level draws
#' level_2:  coda::mcmc.list object containing draws of heterogeneity parameters
#' @import coda Rcpp parallel
#' @export
#' @examples
#' #params <- list(r=1.4, alpha=1.3, s=0.7, beta=7)
#' #cbs <- pcnbd.GenerateData(1000, 10, 5, params)$cbs
#' #draws <- pnbd.mcmc.DrawParameters(cbs, mcmc=10000, burnin=10000, thin=10, chains=2)
#' #plot(draws$level_2)
#' #rbind("actual"=unlist(params), "estimated"=summary(draws$level_2, quantiles=0.5)$quantiles)
#' @seealso pcnbd.GenerateData
#' @references http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=4344404
#' @references http://gradworks.umi.com/34/02/3402149.html
pnbd.mcmc.DrawParameters <-
  function(data,
           mcmc = 10000, burnin = 0, thin = 1, chains = 1,
           use_data_augmentation = TRUE,
           param_init = list(r=1, alpha=1, s=1, beta=1),
           hyper_prior = list(r_1=1/1000, r_2=1/1000,
                              alpha_1=1/1000, alpha_2=1/1000,
                              s_1=1/1000, s_2=1/1000,
                              beta_1=1/1000, beta_2=1/1000)) {
    
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
    
    run_single_chain <- function(chain_id=1, df) {
      
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
        draw_lambda <- if (use_data_augmentation) draw_lambda_conoor else draw_lambda_ma_liu
        draw_mu     <- if (use_data_augmentation) draw_mu_conoor else draw_mu_ma_liu
        level_1["lambda", ] <- draw_lambda(df, level_1, level_2)
        level_1["mu", ]     <- draw_mu(df, level_1, level_2)
        level_1["tau", ]    <- draw_tau(df, level_1)
        
        # draw heterogeneity parameters
        lambda_dist      <- draw_gamma_params("lambda", level_1, level_2, hyper_prior)
        level_2["r"]     <- lambda_dist[1]
        level_2["alpha"] <- lambda_dist[2]
        mu_dist          <- draw_gamma_params("mu", level_1, level_2, hyper_prior)
        level_2["s"]     <- mu_dist[1]
        level_2["beta"]  <- mu_dist[2]
      }
      
      # convert MCMC draws into coda::mcmc objects
      return(list(
        level_1 = lapply(1:nr_of_cust, function(i) mcmc(level_1_draws[,,i], start=burnin, thin=thin)),
        level_2 = mcmc(level_2_draws, start=burnin, thin=thin)))
    }
    
    # run multiple chains
    draws <- mclapply(1:chains, function(i) run_single_chain(i, data), mc.cores=1)
    
    # merge chains into code::mcmc.list objects
    return(list(
      level_1 = lapply(1:nrow(data), function(i) mcmc.list(lapply(draws, function(draw) draw$level_1[[i]]))),
      level_2 = mcmc.list(lapply(draws, function(draw) draw$level_2))))
  }
