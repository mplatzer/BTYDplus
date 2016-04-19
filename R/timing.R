
#' Estimate Regularity in Intertransaction Timings
#' 
#' Estimates degree of regularity of intertransaction timings of a customer cohort,
#' 
#' This is done
#'  - by either assuming same regularity across all customers; this then only 
#'      require three transactions per customer
#'      method: wheat
#'  - or by estimating regularity for each customer seperately (as the shape
#'      parameter of a fitted gamma distribution), and then return the 
#'      median across estimates; this requires min. 10 transactions per customer
#'      methods: mle, mle-minka, mle-thom, cv
#' 
#' @param elog data.frame with transaction logs; requires columns customer-id
#'   'cust' and transaction time \code{t} or 'date'
#' @param method \code{wheat}, \code{mle}, \code{mle-minka}, \code{mle-thom}, \code{cv}
#' @param plot if \code{TRUE} then distribution of estimated regularity will be plotted
#' @return estimated real-valued regularity parameter; rounded to an integer,
#'   this can be used as \code{k} for estimating MBG/CNBD-k models
#' @export
#' @seealso \code{\link{mbgcnbd.EstimateParameters}}
#' @example demo/timing.R
estimateRegularity <- function(elog, method = "wheat", plot = FALSE) {
  if (!"cust" %in% names(elog)) 
    stop("Error in estimateRegularity: elog must have a column labelled \"cust\"")
  if (!"date" %in% names(elog) & !"t" %in% names(elog)) 
    stop("Error in estimateRegularity: elog must have a column labelled \"t\" or \"date\"")
  if (!"t" %in% names(elog)) 
    elog$t <- as.numeric(elog$date)
  trans <- split(elog, elog$cust)
  if (method == "wheat") {
    # Wheat, Rita D., and Donald G. Morrison.  'Estimating purchase regularity with two interpurchase times.'
    # Journal of Marketing Research (1990): 87-93.
    M <- unlist(lapply(trans, function(df) {
      itt <- diff(sort(unique(df$t)))
      if (length(itt) > 1) {
        # take last two itt's to capture most recent regularity
        itt2 <- rev(itt)[1:2]
        return(sample(itt2)[1]/sum(itt2))
      }
    }))
    r <- (1 - 4 * var(M))/(8 * var(M))
    if (plot) {
      op <- par(mar = c(1, 2, 1, 2))
      plot(density(M), main = "", sub = "", xlab = "", ylab = "", lwd = 2, frame = FALSE, axes = FALSE)
      polygon(density(M), col = "lightgray", border = 1)
      fn1 <- function(x) dbeta(x, 1, 1)
      fnr <- function(x) dbeta(x, round(r), round(r))
      curve(fn1, add = TRUE, lty = 2, lwd = 2)
      curve(fnr, add = TRUE, lty = 2, lwd = 2)
      par(op)
    }
    return(r)
    
  } else {
    if (method == "mle" | method == "mle-minka") {
      # Maximum Likelihood Estimator http://en.wikipedia.org/wiki/Gamma_distribution#Maximum_likelihood_estimation
      # Approximation for MLE by Minka http://research.microsoft.com/en-us/um/people/minka/papers/minka-gamma.pdf
      ks <- unlist(lapply(trans, function(df) {
        itt <- diff(sort(unique(df$t)))
        if (length(itt) >= 9) {
          s <- log(sum(itt)/length(itt)) - sum(log(itt))/length(itt)
          if (method == "mle") {
          fn <- function(v) {
            return((log(v) - digamma(v) - s)^2)
          }
          k <- optimize(fn, lower = 0.1, upper = 50)$min
          } else if (method == "mle-minka") {
          k <- (3 - s + sqrt((s - 3)^2 + 24 * s))/(12 * s)
          }
          return(k)
        }
      }))
      
    } else if (method == "mle-thom") {
      # Approximation for ML estimator Thom (1968); see Dunn, Richard, Steven Reader, and Neil Wrigley.  'An
      # investigation of the assumptions of the NBD model' Applied Statistics (1983): 249-259.
      ks <- unlist(lapply(trans, function(df) {
        itt <- diff(sort(unique(df$t)))
        if (length(itt) >= 9) {
          hm <- function(v) exp(sum(log(v))/length(v))
          mu <- log(mean(itt)/hm(itt))
          d <- (1/(4 * mu)) * (1 + sqrt(1 + 4 * mu/3))
          return(d)
        }
      }))
      
    } else if (method == "cv") {
      # Estimate regularity by analyzing coefficient of variation Wu, Couchen, and H-L. Chen. 'A consumer purchasing
      # model with learning and departure behaviour.'  Journal of the Operational Research Society (2000): 583-591.
      ks <- unlist(lapply(trans, function(df) {
        itt <- diff(sort(unique(df$t)))
        if (length(itt) >= 9) {
          cv <- sd(itt)/mean(itt)
          k <- 1/cv^2
          return(k)
        }
      }))
    }
    if (length(ks) == 0) 
      stop("No customers with 10 or more transactions.")
    
    if (plot) {
      ymax <- median(ks) * 3
      boxplot(ks, horizontal = TRUE, ylim = c(0, ymax), frame = FALSE, axes = FALSE)
      axis(1, at = 0:ymax)
      axis(3, at = 0:ymax, labels = rep("", 1 + ymax))
      abline(v = 1:ymax, lty = "dotted", col = "lightgray")
      boxplot(ks, horizontal = TRUE, add = TRUE, col = "gray", frame = FALSE, axes = FALSE)
    }
    
    return(median(ks))
  }
} 
