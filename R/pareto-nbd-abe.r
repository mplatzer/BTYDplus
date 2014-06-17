# 
# # ! not yet ready for use !
# 
# ## Implementation of Pareto/NBD HB as describe in Abe (2009) and Korkmaz (2013)
# ## FIXME: implement passing of covariates
# ## FIXME: check why we gamma_draws are off
# ## TODO: estimate CDnow and compare to Abe/Korkmaz results
# 
# # FIXME docs
# abe.mcmc.DrawParameters <- function(cal.cbs, mcmc = 20000, burnin = mcmc/2, thin = 1) {
# 
#   draw_z <- function(data, theta) {
#     tx     <- data[, "t.x"]
#     Tcal   <- data[, "T.cal"]
#     lambda <- theta["lambda",]
#     mu     <- theta["mu",]
#     prob   <- 1/(1+(mu/(lambda+mu))*(exp((lambda+mu)*(Tcal-tx))-1))
#     prop   <- pmin(prob, 0.95) # impose upper limit; see Note (3)
#     z      <- as.numeric(runif(length(prob)) < prob)
#     return(z)
#   }
#   
#   draw_y <- function(data, theta) {
#     # sample from truncated exponential distribution  
#     lambda <- theta["lambda",]
#     mu     <- theta["mu",]
#     mu_lam <- mu + lambda
#     z      <- theta["z",]
#     tx     <- data[, "t.x"]
#     Tcal   <- data[, "T.cal"]
#     
#     alive  <- z==1
#     y      <- numeric(nrow(data))
#     
#     # Case: still alive - left truncated exponential distribution -> [T.cal, Inf]
#     if (any(alive)) {
#       y[alive]  <- Tcal[alive] + rexp(sum(alive), mu_lam[alive])
#     }
#     
#     # Case: churned     - double truncated exponential distribution -> [tx, T.cal]
#     if (any(!alive)) {
#       mu_lam_tx   <- pmin(700, mu_lam[!alive] * tx[!alive])
#       mu_lam_Tcal <- pmin(700, mu_lam[!alive] * Tcal[!alive])
#       rand        <- runif(n=sum(!alive))
#       y[!alive] <- -log( (1-rand)*exp(-mu_lam_tx) + rand*exp(-mu_lam_Tcal)) / mu_lam[!alive]
#     }
#     return(y)
#   }
#   
#   draw_beta_gamma <- function(data, theta) {
#     # standard multi-variate normal regression update
#     draw <- bayesm::rmultireg(Y=log(t(theta[c("lambda", "mu"),])),
#                               X=covariates, Bbar=beta_0, A=A_0, nu=nu_00, V=gamma_00)
#     return(list(beta=t(draw$B), gamma=draw$Sigma))
#   }
#   
#   draw_theta <- function(data, covariates, theta, beta, gamma) {
#     # sample (lambda, mu) given (z, y, beta, gamma)
#     N      <- nrow(data)
#     x      <- data[, "x"]
#     tx     <- data[, "t.x"]
#     Tcal   <- data[, "T.cal"]
#     z      <- theta["z", ]
#     y      <- theta["y", ]
#     mvmean <- covariates[,] %*% t(beta)
#     inv_gamma <- solve(gamma)
#     
#     cur_lambda <- theta["lambda", ]
#     cur_mu     <- theta["mu", ]
#     
#     log_post <- function(log_theta) {
#       log_lambda  <- log_theta[1,]
#       log_mu      <- log_theta[2,]
#       diff_lambda <- log_lambda - mvmean[,1]
#       diff_mu     <- log_mu     - mvmean[,2]
#       likel <- x*log_lambda + (1-z)*log_mu - (exp(log_lambda)+exp(log_mu))*(z*Tcal+(1-z)*y)
#       prior <- -0.5 * (diff_lambda^2*inv_gamma[1,1] + 2*diff_lambda*diff_mu*inv_gamma[1,2] + diff_mu^2*inv_gamma[2,2])
#       post <- likel + prior
#       post[log_mu>5] <- -Inf # cap !!
#       return(post)
#     }
#     
#     # current state
#     cur_log_theta <- rbind(log(cur_lambda), log(cur_mu))
#     cur_post      <- log_post(cur_log_theta)
#     
#     step <- function(cur_log_theta, cur_post) {
#       # new proposal
#       new_log_theta     <- cur_log_theta + rbind(gamma[1,1]*rt(N, df=3), gamma[2,2]*rt(n=N, df=3))
#       new_log_theta[1,] <- pmax(pmin(new_log_theta[1,], 70), -70)
#       new_log_theta[2,] <- pmax(pmin(new_log_theta[2,], 70), -70)
#       new_post          <- log_post(new_log_theta)
#       
#       # accept/reject new proposal
#       mhratio  <- exp(new_post-cur_post)
#       accepted <- mhratio > runif(n=N)
#       
#       cur_log_theta[,accepted] <- new_log_theta[,accepted]
#       cur_post[accepted] <- new_post[accepted]
#       
#       list(cur_log_theta=cur_log_theta, cur_post=cur_post)
#     }
#     
#     iter <- 1 # how high do we need to set this? 1/5/10/100?
#     for (i in 1:iter) {
#       draw <- step(cur_log_theta, cur_post)
#       cur_log_theta <- draw$cur_log_theta
#       cur_post <- draw$cur_post
#     }
#     cur_theta <- exp(cur_log_theta)
#     
#     return(list(lambda=cur_theta[1,], mu=cur_theta[2,]))
#   }
#   
#   data <- cal.cbs
#   
#   # Setup Regressors (Covariates) for location of 1st-stage prior, i.e. beta = [log(lambda), log(mu)]
#   covariates <- matrix(1, nrow=nrow(data), ncol=1) # ncol=1 # we use a single intercept as regressor here
#   colnames(covariates) <- c("(Intercept)")
#   K <- ncol(covariates) # number of covariates
#   
#   ## initialize arrays for storing draws ##
#   
#   nr_of_draws <- (mcmc-1) %/% thin + 1
#   
#   beta_draws <- array(NA_real_, dim=c(nr_of_draws, 2, K))
#   dimnames(beta_draws)[[2]] <- c("lambda", "mu")
#   dimnames(beta_draws)[[3]] <- colnames(covariates)
#   
#   gamma_draws <- array(NA_real_, dim=c(nr_of_draws, 2, 2))
#   dimnames(gamma_draws)[[2]] <- c("lambda", "mu")
#   dimnames(gamma_draws)[[3]] <- c("lambda", "mu")
#   
#   theta_draws <- array(NA_real_, dim=c(nr_of_draws, 4, nrow(data)))
#   dimnames(theta_draws)[[2]] <- c("lambda", "mu", "y", "z")
#   
#   ## initialize parameters ##
#   
#   theta            <- theta_draws[1,,]
#   theta["lambda",] <- mean(data$x) / mean(ifelse(data$t.x==0, data$T.cal, data$t.x))
#   theta["mu",]     <- 1/mean(data$T.cal)
#   
#   # set hyper-parameters of 2nd-stage prior of beta
#   beta_0 <- matrix(0, nrow=K, ncol=2, dimnames=list(NULL, c("log_lambda", "log_mu")))
#   beta_0[1, "log_lambda"] <- log(mean(theta["lambda",]))
#   beta_0[1, "log_mu"]     <- log(mean(theta["mu",]))
#   A_0                     <- diag(rep(.01, K), ncol=K, nrow=K) # diffuse precision matrix
#   # Note: the Abe paper specifies sigma_0 independent of gamma_0, but the standard prior
#   #  for multivariate regression has in fact a dependent distribution 'beta given sigma_0'
#   #  and models the relationship via a precision matrix A
#   
#   # set diffuse hyper-parameters for 2nd-stage prior of gamma_0
#   # - follows defaults from rmultireg example
#   nu_00    <- 3 + K # 30
#   gamma_00 <- nu_00 * diag(2)  
#   
#   
#   ## run MCMC chain ##
#   
#   for (step in 1:(burnin+mcmc)) {
#     if (step%%1000==0) cat(step, "\n")
#     
#     # draw individual-level parameters
#     theta["z", ] <- draw_z(data, theta)
#     theta["y", ] <- draw_y(data, theta)
#     
#     draw  <- draw_beta_gamma(data, theta)
#     beta  <- draw$beta
#     gamma <- draw$gamma
#     
#     draw <- draw_theta(data, covariates, theta, beta, gamma)
#     theta["lambda", ] <- draw$lambda
#     theta["mu", ]     <- draw$mu
#   
#     # store
#     if ((step-burnin)>0 & (step-1-burnin)%%thin==0) {
#       idx <- (step-1-burnin)%/%thin + 1
#       theta_draws[idx,,] <- theta
#       beta_draws[idx,,]   <- beta
#       gamma_draws[idx,,]  <- gamma
#     }
#   }
# 
#   return(list(theta=theta_draws,
#               beta=beta_draws,
#               gamma=gamma_draws))
# }
# 
# 
# #' Generate artificial data which follows Pareto/CNBD model assumptions
# #'
# #' Returns 2-element list
# #' * cbs: data.frame with 'cust', \code{x}, \code{t.x}, \code{T.cal}, 'T.star', 'x.star' 
# #'        this is the summary statistics data.frame which contains all 
# #'        needed information for parameter estimation
# #' * elog: data.frame with 'cust', \code{t}
# #'
# #' @param n number of customers
# #' @param T.cal length of calibration period
# #' @param T.star length of holdout period
# #' @param params list of parameters: {r, alpha, s, beta, t, gamma}
# #' @param return.elog if \code{TRUE} then the event-log is returned as well; decreases performance
# #' 
# #' @return 2-elemnt list
# #' @export
# #' @import mvtnorm
# abe.GenerateData <- function (n, T.cal, T.star, params, return.elog=FALSE) {
#   
#   if (length(T.cal)==1) T.cal <- rep(T.cal, n)
#   if (length(T.star)==1) T.star <- rep(T.star, n)
#   
#   # sample intertransaction timings parameter lambda for each customer
#   thetas  <- exp(mvtnorm::rmvnorm(n, mean=params$beta, sigma=params$gamma))
#   lambdas <- thetas[,1]
#   mus     <- thetas[,2]
#   
#   # sample lifetime for each customer
#   taus <- rexp(n, rate=mus)
#   
#   # sample intertransaction timings & churn
#   cbs_list <- list()
#   elog_list <- list()
#   for (i in 1:n) {
#     lambda <- lambdas[i]
#     mu <- mus[i]
#     tau <- taus[i]
#     # sample 'sufficiently' large amount of inter-transaction times
#     minT <- min(T.cal[i]+T.star[i], tau)
#     nr_of_itt_draws <- max(10, minT * lambda)
#     itts <- rexp(nr_of_itt_draws * 2, rate=lambda)
#     if (sum(itts)<minT) itts <- c(itts, rexp(nr_of_itt_draws * 4, rate=lambda))
#     if (sum(itts)<minT) itts <- c(itts, rexp(nr_of_itt_draws * 800, rate=lambda))
#     if (sum(itts)<minT)
#       stop("not enough inter-transaction times sampled: ", sum(itts), " < ", tau)
#     times <- cumsum(c(0, itts))
#     times <- times[times<tau]
#     if (return.elog)
#       elog_list[[i]] <- data.frame(cust=i, t=times[times<(T.cal[i]+T.star[i])])
#     # determine frequency, recency, etc.
#     ts.cal   <- times[times<T.cal[i]]
#     ts.star  <- times[times>=T.cal[i] & times<(T.cal[i]+T.star[i])]
#     cbs_list[[i]] <- list(cust   = i,
#                           x      = length(ts.cal)-1,
#                           t.x    = max(ts.cal),
#                           litt   = ifelse(length(ts.cal)-1==0, 0, sum(log(itts[1:(length(ts.cal)-1)]))),
#                           alive  = tau>T.cal[i],
#                           x.star = length(ts.star))
#   }
#   cbs <- do.call(rbind.data.frame, cbs_list)
#   cbs$lambda <- lambdas
#   cbs$mu     <- mus
#   cbs$tau    <- taus
#   cbs$T.cal  <- T.cal
#   cbs$T.star <- T.star
#   rownames(cbs) <- NULL
#   out <- list(cbs=cbs)
#   if (return.elog) {
#     elog <- do.call(rbind.data.frame, elog_list)
#     out$elog <- elog
#   }
#   return(out)
# }
# 
# 
# 
# 
# abe.mcmc.PAlive <- function (data, draws) {
#   p.alives <- colMeans(draws$theta[, "z", ])
#   return(p.alives)
# }
# 
# 
# abe.mcmc.DrawFutureTransactions <- function (data, draws, T.star) {
#   x.stars <- array(NA_real_, dim=dim(draws$theta)[c(1,3)])
#   for (cust in 1:nrow(data)) {
#     T.cal   <- data[cust, "T.cal"]
#     alive   <- as.logical(draws$theta[, "z", cust])
#     ys      <- draws$theta[, "y", cust]
#     lambdas <- draws$theta[, "lambda", cust]
#     # Case: customer alive
#     if (any(alive)) {
#       x.stars[alive, cust] <- rpois(sum(alive), lambdas[alive] * pmin(ys[alive] - T.cal, T.star))
#     }
#     # Case: customer churned
#     if (any(!alive)) {
#       x.stars[!alive, cust] <- 0
#     }
#   }
#   return(x.stars)
# }
# 
# 
# abe.test <- function() {
#   
#   params <- list()
#   params$beta <- c(log(1.2), log(0.08))
#   params$gamma <- matrix(c(0.02, 0, 0, 0.33), ncol=2)
#   set.seed(1)
#   data <- abe.GenerateData(n=1000, T.cal=9, T.star=30, params)$cbs
#   draws <- abe.mcmc.DrawParameters(data)
#   
#   # plot convergence
#   
#   plot(apply(draws$theta[,"mu",],1,median), typ="l")
#   abline(h=exp(params$beta[2]), col="red")
#   
#   plot(apply(draws$theta[,"lambda",],1,median), typ="l")
#   abline(h=exp(params$beta[1]), col="red")
#   
#   plot(draws$beta[,1,1], typ="l")
#   abline(h=params$beta[1], col="red")
#   plot(draws$beta[,2,1], typ="l")
#   abline(h=params$beta[2], col="red")
#   
#   plot(draws$gamma[,1,1], typ="l")
#   abline(h=params$gamma[1,1], col="red")
#   # why do we miss gamma[1,1] completely??
#   
#   plot(draws$gamma[,2,2], typ="l")
#   abline(h=params$gamma[2,2], col="red")
#   
#   plot(draws$gamma[,1,2], typ="l")
#   abline(h=params$gamma[1,2], col="red")
#   # why do we get a positive correlation?
#   
# }
# 
# abe.cdnow <- function() {
#   
#   data(cdnowSummary, package="BTYD")
#   data <- as.data.frame(cdnowSummary$cbs)
#   set.seed(1)
#   draws <- abe.mcmc.DrawParameters(data, burnin=10000, mcmc=50000, thin=100)
#   
#   agg <- function(x) round(c(quantile(x, 0.025), "mean"=mean(x), quantile(x, 0.975)),2)
#   agg(rowMeans(log(draws$theta[, "lambda",])))
#   agg(rowMeans(log(draws$theta[, "mu",])))
#   agg((as.numeric(draws$gamma[, 1,1])))
#   agg((as.numeric(draws$gamma[, 2,2])))
#   agg((as.numeric(draws$gamma[, 1,2])))
#   # numbers match Abe paper
#   
#   mcmc.est.draws <- abe.mcmc.DrawFutureTransactions(data, draws, T.star=39)
#   mcmc.palive    <- abe.mcmc.PAlive(data, draws)
#   mcmc.est       <- colMeans(mcmc.est.draws)
#   #mcmc.est       <- apply(mcmc.est.draws, 2, median)
#   mcmc.est.pos   <- colMeans(mcmc.est.draws>0)
#   
#   cor(data$x.star, mcmc.est)
#   # 0.52 (should be 0.6245)
#   
#   mse <- function(x1, x2) sqrt(mean((x1-x2)^2))
#   mse(data$x.star, mcmc.est)
#   # 2.036 (should be 2.606)
#   
#   mae <- function(x1, x2) mean(abs(x1-x2))
#   mae(data$x.star, mcmc.est)
#   # 0.801 (should be 0.758)
#   
#   # very low estimates? plus numbers don't match!?  -> do we need to have a tighter prior?
#   
#   
#   plot(function(x) dgamma(x, 0.553, 10.580))
#   lines(density(draws$theta[, "lambda", ]), col="red")
#   breaks <- c(seq(0,0.04,len=100), 100)
#   hist(rgamma(100000, 0.553, 10.580), breaks, xlim=c(0,0.04), col=rgb(0,0,1,1/4), prob=T)
#   hist(draws$theta[, "lambda", ], breaks, prob=TRUE, add=TRUE, col=rgb(1,0,0,1/4))
#   
#   plot(function(x) dgamma(x, 0.606, 11.656))
#   lines(density(draws$theta[, "mu", ]), col="red")  
# }
