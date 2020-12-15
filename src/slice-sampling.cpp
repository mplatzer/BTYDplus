#include <Rcpp.h>

using namespace Rcpp;

// slice sampling

// This is a stripped down C++ version of diversitree::mcmc thta only returns
// the draw from the last step. Drawing from gamma distribution is now 300x
// faster than compared with diversitree::mcmc.

// The algorithm is described in detail in Neal R.M. 2003. Slice sampling.
// Annals of Statistics 31:705-767. which describes *why* this algorithm works.
// The approach differs from normal Metropolis-Hastings algorithms, and from
// Gibbs samplers, but shares the Gibbs sampler property of every update being
// accepted.  Nothing is required except the ability to evaluate the function at
// every point in continuous parameter space.

// Let x0 be the current (possibly multivariate) position in continuous
// parameter space, and let y0 be the probability at that point.  To update from
// (x0, y0) -> (x1, y1), we update each of the parameters in turn.  For each
// parameter
//
//   1. Draw a random number 'z' on Uniform(0, y0) -- the new point
//      must have at least this probability.
//
//   2. Find a region (x.l, x.r) that contains x0[i], such that x.l
//      and x.r are both smaller than z.
//
//   3. Randomly draw a new position from (x.l, x.r).  If this
//      position is greater than z, this is our new position.
//      Otherwise it becomes a new boundary and we repeat this step
//      (the point x1 becomes x.l if x.l < x0[i], and x.r otherwise so
//      that x0[i] is always contained within the interval).
// Because it is generally more convenient to work with log
// probabilities, step 1 is modified so that we draw 'z' by taking
//   y0 - rexp(1)
// All other steps remain unmodified.

NumericVector slice_sample_cpp(double (*logfn)(NumericVector, NumericVector),
                        NumericVector params,
                        NumericVector x0,
                        int steps = 10,
                        double w = 1,
                        double lower = -INFINITY,
                        double upper = INFINITY) {

  double u, r0, r1, logy, logz, logys;
  NumericVector x, xs, L, R;

  x = clone(x0);
  L = clone(x0);
  R = clone(x0);
  logy = logfn(x, params);

  for (int i = 0; i < steps; i++) {

    for (int j = 0; j < x0.size(); j++) {
      // draw uniformly from [0, y]
      logz = logy - rexp(1)[0];

      // expand search range
      u = runif(1)[0] * w;
      L[j] = x[j] - u;
      R[j] = x[j] + (w-u);
      while ( L[j] > lower && logfn(L, params) > logz )
        L[j] = L[j] - w;
      while ( R[j] < upper && logfn(R, params) > logz )
        R[j] = R[j] + w;

      // sample until draw is within valid range
      r0 = std::max(L[j], lower);
      r1 = std::min(R[j], upper);

      xs = clone(x);
      int cnt = 0;
      do {
        cnt++;
        xs[j] = runif(1, r0, r1)[0];
        logys = logfn(xs, params);
        if ( logys > logz )
          break;
        if ( xs[j] < x[j] )
          r0 = xs[j];
        else
          r1 = xs[j];
      } while (cnt<1e4);
      if (cnt==1e4) ::Rf_error("slice_sample_cpp loop did not finish");

      x = clone(xs);
      logy = logys;
    }
  }

  return x;
}

// // draw from gamma distribution (for test purposes)
//
// double post_gamma(NumericVector x, NumericVector params) {
//   double alpha = params[0];
//   double beta = params[1];
//   return (alpha - 1) * log(x[0]) - beta * x[0];
// }
//
// // [[Rcpp::export]]
// NumericVector slice_sample_gamma(double alpha, double beta, double lower, double upper) {
//   NumericVector params = NumericVector::create(alpha, beta);
//   NumericVector x0 = NumericVector::create(alpha/beta);
//   double steps = 10;
//   double w = 3 * sqrt(alpha) / beta; // approx size of (q95-q05)
//   return slice_sample_cpp(post_gamma, params, x0, steps, w, lower, upper);
// }
//
// // draw from multivariate normal distribution (for test purposes)
//
// double post_mvnorm(NumericVector x, NumericVector sigma) {
//   return -log(2*3.141593) -0.5 * log(sigma[0]*sigma[3]-sigma[1]*sigma[2]) -0.5 * (1/(sigma[0]*sigma[3]-sigma[1]*sigma[2])) *
//     (x[0]*x[0]*sigma[3] - x[0]*x[1]*sigma[2] - x[0]*x[1]*sigma[1] + x[1]*x[1]*sigma[0]);
// }
//
// // [[Rcpp::export]]
// NumericVector slice_sample_mvnorm(NumericVector sigma) {
//   NumericVector x0 = NumericVector::create(0.2, 0.3);
//   double steps = 20;
//   return slice_sample_cpp(post_mvnorm, sigma, x0, steps);
// }

// estimate parameters of gamma distribution

double post_gamma_parameters(NumericVector log_data, NumericVector params) {
  double shape = exp(log_data[0]);
  double rate = exp(log_data[1]);
  double len_x = params[0];
  double sum_x = params[1];
  double sum_log_x = params[2];
  double hyper1 = params[3];
  double hyper2 = params[4];
  double hyper3 = params[5];
  double hyper4 = params[6];
  return len_x * (shape * log(rate) - lgamma(shape)) + (shape-1) * sum_log_x - rate * sum_x +
    (hyper1 - 1) * log(shape) - (shape * hyper2) +
    (hyper3 - 1) * log(rate) - (rate * hyper4);
}

// [[Rcpp::export]]
NumericVector slice_sample_gamma_parameters(NumericVector data, NumericVector init,
                                            NumericVector hyper, double steps = 20, double w = 1) {
  NumericVector params = NumericVector::create(data.size(), sum(data), sum(log(data)),
                                                hyper[0], hyper[1], hyper[2], hyper[3]);
  return exp(slice_sample_cpp(post_gamma_parameters, params, log(init), steps, w, -INFINITY, INFINITY));
}

/*** R

  # unit-test slice sampling for univariate distribution
  alpha <- 2
  beta <- 5
  n <- 1e4
  lower <- 0.3
  upper <- 0.8
  draws1 <- sapply(1:n, function(i) slice_sample_gamma(alpha, beta, lower, upper))
  draws2 <- rgamma(n, alpha, beta)
  draws2 <- draws2[draws2>lower & draws2<upper]
  stopifnot(abs(mean(draws1)-mean(draws2))<.1)

  # unit-test slice sampling for bivariate distribution
  n <- 1e4
  sigma <- c(1, 0.6, 0.6, 1.2)
  draws1 <- t(sapply(1:n, function(i) slice_sample_mvnorm(sigma)))
  draws2 <- MASS::mvrnorm(n, mu=c(0,0), Sigma=matrix(sigma, ncol=2))
  stopifnot(mean(abs(cov(draws2) - cov(draws1)))<.3)

  # unit-test slice_sample_gamma_parameters
  n <- 1e4
  params <- c(1.4, 3.5)
  x <- rgamma(1e4, params[1], params[2])
  draws <- t(replicate(n, slice_sample_gamma_parameters(x, c(1,1), rep(1e-3, 4))))
  stopifnot(max(abs(apply(draws, 2, mean) - params))<.1)
*/


//
// working R implementation of slice sampling
//
//slice_sample <- function(logfn, x0, steps, w, lower=-Inf, upper=Inf) {
//
//  x <- x0
//  logy <- logfn(x0)
//
//  for (i in 1:steps) {
//    # draw uniformly from [0, y]
//    logz <- logy - rexp(1)
//
//    # expand search range
//    u <- runif(1) * w
//    L <- x - u
//    R <- x + (w-u)
//    while ( L > lower && logfn(L) > logz )
//      L <- L - w
//    while ( R < upper && logfn(R) > logz )
//      R <- R + w
//
//    # sample until draw is within valid range
//    r0 <- max(L, lower)
//    r1 <- min(R, upper)
//    repeat {
//      xs <- runif(1, r0, r1)
//      logys <- logfn(xs)
//      if ( logys > logz ) break
//      else if ( xs < x )  r0 <- xs
//      else                r1 <- xs
//    }
//    x <- xs
//    logy <- logys
//  }
//
//  return(x)
//}



// ********* Pareto / GGG **********


// draw of individual-level posterior for Pareto/GGG

double pggg_palive_integrand(double x, double params[6]) {
  double tx = params[1];
  double k = params[3];
  double lambda = params[4];
  double mu = params[5];
  // call pgamma with lower.tail=FALSE and log.p=FALSE; note that 3rd parameter is scale and not rate;
  return (::Rf_pgamma(x-tx, k, 1/(k*lambda), 0, 0) * exp(-mu*x));
}

double simpson38(double (*fn)(double, double[6]), double a, double b, double fn_params[6]) {
  // https://en.wikipedia.org/wiki/Simpson%27s_rule#Simpson.27s_3.2F8_rule_.28for_n_intervals.29
  double n = 12.0;
  double integral = (3.0/8.0) * ((b-a)/n) *
    (fn(a, fn_params) +
    3 * fn(a+(1/n)*(b-a), fn_params) +
    3 * fn(a+(2/n)*(b-a), fn_params) +
    2 * fn(a+(3/n)*(b-a), fn_params) +
    3 * fn(a+(4/n)*(b-a), fn_params) +
    3 * fn(a+(5/n)*(b-a), fn_params) +
    2 * fn(a+(6/n)*(b-a), fn_params) +
    3 * fn(a+(7/n)*(b-a), fn_params) +
    3 * fn(a+(8/n)*(b-a), fn_params) +
    2 * fn(a+(9/n)*(b-a), fn_params) +
    3 * fn(a+(10/n)*(b-a), fn_params) +
    3 * fn(a+(11/n)*(b-a), fn_params) +
    fn(b, fn_params));
  return(integral);
}

// [[Rcpp::export]]
NumericVector pggg_palive(NumericVector x, NumericVector tx, NumericVector Tcal,
                           NumericVector k, NumericVector lambda, NumericVector mu) {
  int N = x.size();
  NumericVector out(N);
  for (int i=0; i<N; i++) {
    // calc numerator
    double one_minus_F = ::Rf_pgamma(Tcal[i]-tx[i], k[i], 1/(k[i]*lambda[i]), 0, 0);
    double numer = one_minus_F * exp(-mu[i]*Tcal[i]);
    // calc denominator by integrating from tx to Tcal
    // - we integrate numerically via Simpson3/8 rule (calling Rdqags crashed under Unix)
    double fn_params[6] = {x[i], tx[i], Tcal[i], k[i], lambda[i], mu[i]};
    double integral = simpson38(pggg_palive_integrand, tx[i], Tcal[i], fn_params);
    double denom = numer + mu[i] * integral;
    out[i] = (numer/denom);
  }
  return(out);
}


double pggg_post_tau(NumericVector data, NumericVector params) {
  double tau_    = data[0];
//  double x      = params[0];
//  double tx     = params[1];
//  double Tcal   = params[2];
//  double litt   = params[3];
  double k      = params[4];
  double lambda = params[5];
  double mu     = params[6];
//  double tau    = params[7];
//  double t      = params[8];
//  double gamma  = params[9];
//  double r      = params[10];
//  double alpha  = params[11];
//  double s      = params[12];
//  double beta   = params[13];

  return(-mu*tau_ + ::Rf_pgamma(tau_, k, 1/(k*lambda), 0, 1));
}


double pggg_post_k(NumericVector data, NumericVector params) {
  double k_     = data[0];
  double x      = params[0];
  double tx     = params[1];
  double Tcal   = params[2];
  double litt   = params[3];
//  double k      = params[4];
  double lambda = params[5];
//  double mu     = params[6];
  double tau    = params[7];
  double t      = params[8];
  double gamma  = params[9];
//  double r      = params[10];
//  double alpha  = params[11];
//  double s      = params[12];
//  double beta   = params[13];

  double log_one_minus_F = ::Rf_pgamma(std::min(Tcal, tau) - tx, k_, 1/(k_*lambda), 0, 1);
  return (t-1) * log(k_) - (k_*gamma) +
    k_ * x * log(k_*lambda) - x * lgamma(k_) - k_ * lambda * tx + (k_-1) * litt +
    log_one_minus_F;
}

double pggg_post_lambda(NumericVector data, NumericVector params) {
  double lambda_ = data[0];
  double x       = params[0];
  double tx      = params[1];
  double Tcal    = params[2];
//  double litt    = params[3];
  double k       = params[4];
//  double lambda  = params[5];
//  double mu      = params[6];
  double tau     = params[7];
//  double t       = params[8];
//  double gamma   = params[9];
  double r       = params[10];
  double alpha   = params[11];
//  double s       = params[12];
//  double beta    = params[13];

  double log_one_minus_F = ::Rf_pgamma(std::min(Tcal, tau) - tx, k, 1/(k*lambda_), 0, 1);
  return (r-1) * log(lambda_) - (lambda_*alpha) +
    k * x * log(lambda_) - k * lambda_ * tx +
    log_one_minus_F;
}


// [[Rcpp::export]]
NumericVector pggg_slice_sample(String what,
                                  NumericVector x, NumericVector tx, NumericVector Tcal, NumericVector litt,
                                  NumericVector k, NumericVector lambda, NumericVector mu, NumericVector tau,
                                  double t, double gamma, double r, double alpha, double s, double beta) {
  int N = x.size();
  NumericVector out(N);

  for (int i=0; i<N; i++) {
    NumericVector params = NumericVector::create(
                              x[i], tx[i], Tcal[i], litt[i],
                              k[i], lambda[i], mu[i], tau[i],
                              t, gamma, r, alpha, s, beta);
    //Rcpp::Rcout << i << " - " << x[i] << " - " << tx[i] << " - " << Tcal[i] << " - " << litt[i] << " - " << k[i] << " - " << lambda[i] << " - " << mu[i] << " - " << tau[i] << " - " << t << " - " << gamma << " - " << r << " - " << alpha << " - " << s << " - " << beta << " - " << std::endl;
    if (what == "k") {
      out[i] = slice_sample_cpp(pggg_post_k, params, NumericVector::create(k[i]), 3, 3 * sqrt(t) / gamma, 1e-1, 1e+3)[0];
    } else if (what == "lambda") {
      out[i] = slice_sample_cpp(pggg_post_lambda, params, NumericVector::create(lambda[i]), 3, 3 * sqrt(r) / alpha, 1e-30, 1e+5)[0];
    } else if (what == "tau") {
      double tau_init = std::min(Tcal[i]-tx[i], ::Rf_rgamma(k[i], 1/(k[i]*lambda[i]))) / 2;
      out[i] = tx[i] + slice_sample_cpp(pggg_post_tau, params, NumericVector::create(tau_init), 6, (Tcal[i]-tx[i])/2, 0, Tcal[i]-tx[i])[0];
    }
  }
  return out;
}


/*** R
  # unit-test slice sampling of pggg_post_tau, by comparing results to pareto/nbd (k=1),
  #   where we can draw directly via https://en.wikipedia.org/wiki/Inverse_transform_sampling
  x <- 0
  tx <- 8
  Tcal <- 14
  litt <- 0
  k <- 1
  lambda <- 1.2
  mu <- 0.01
  n <- 10^4
  draws1 <- pggg_slice_sample("tau", rep(x, n), rep(tx, n), rep(Tcal, n), rep(litt, n),
                               rep(k, n), rep(lambda, n), rep(mu, n), rep(0, n),
                               1,1,1,1,1,1)
  rand <- runif(n)
  draws2 <- -log( (1-rand)*exp(-(mu+lambda) * tx) + rand*exp(-(mu+lambda) * Tcal)) / (mu+lambda)
  err <- abs(mean(draws1)-mean(draws2))
  stopifnot(err<.1)

  # compare results for k!=1 with random draws
  k <- 4
  lambda <- 0.0000001
  mu <- 0.1
  Tcal <- 52
  tx <- 26
  draws1 <- replicate(1e5, pggg_slice_sample("tau", x, tx, Tcal, 0, k, lambda, mu, 0, 0,0,0,0,0,0))
  draws2 <- pmin(rexp(1e7, mu), rgamma(1e7, k, k*lambda))
  draws2 <- draws2[draws2>tx & draws2<Tcal]
  err <- abs(mean(draws1)-mean(draws2))
  stopifnot(err<.1)
*/

/*** R
  # unit-test P(alive) by comparing C++ result matches R results, matches Pareto/NBD result (k=1)
  x <- 1
  tx <- 45
  Tcal <- 52
  k <- 3
  lambda <- 0.0001
  mu <- 0.1

  # C++ implementation
  res1 <- pggg_palive(x, tx, Tcal, k, lambda, mu)
  # R implementation
  ff <- function(y, tx, k, lambda, mu) (1-pgamma(y-tx, k, k*lambda)) * exp(-mu*y)
  res2 <- ff(Tcal, tx, k, lambda, mu) / (ff(Tcal, tx, k, lambda, mu) +
                                   mu * integrate(ff, lower=tx, upper=Tcal, tx=tx, k=k, lambda=lambda, mu=mu)$value)
  # R implementation for k=1
  la_mu <- lambda + mu
  res3 <- exp(-la_mu*Tcal) / (exp(-la_mu*Tcal) + (mu/la_mu)*(exp(-la_mu*tx)-exp(-la_mu*Tcal)))
  # determine P(alive) via sampling
  itts <- rgamma(1e8, k, k*lambda)
  taus <- rexp(1e8, mu)
  res4 <- sum(itts>(Tcal-tx) & taus>Tcal) / (sum(itts>(Tcal-tx) & taus>Tcal) + sum(itts>(taus-tx) & taus<Tcal & taus>tx))

  ape <- function(a, e) abs(a-e)/a
  stopifnot(ape(res1, res2) < 0.01)
  stopifnot(ape(res1, res4) < 0.01)
  stopifnot(ape(res2, res4) < 0.01)
  stopifnot(k!=1 | ape(res1, res2) < 0.01)
*/


// ********* BG/CNBD-k **********

// [[Rcpp::export]]
double xbgcnbd_pmf_cpp(NumericVector params, double t, int x, bool dropout_at_zero = false) {
  if (params.size() != 5) ::Rf_error("params needs to be of size 5 with (k, r, alpha, a, b)");
  if (t == 0) return(0);
  int k        = params[0];
  double r     = params[1];
  double alpha = params[2];
  double a     = params[3];
  double b     = params[4];
  double survivals = x - 1;
  if (dropout_at_zero) survivals = survivals + 1;
  double P1 = exp(lgamma(b+survivals+1)+lgamma(a+b)-lgamma(b)-lgamma(a+b+survivals+1));
  double P2a = 0;
  for (int i = k*x; i <= k*x+k-1; i++) {
    P2a += exp(lgamma(r+i) + r*log(alpha) + i*log(t) - lgamma(i+1) - lgamma(r) - (r+i)*log(alpha+t));
  }
  double P2b;
  if ((dropout_at_zero==FALSE) & (x==0)) {
    P2b = 0;
  } else {
    P2b = a / (b+survivals);
    if (x>0) {
      double cmf = 0;
      for (int i = 0; i <= k*x-1; i++) {
        cmf += exp(lgamma(r+i) + r*log(alpha) + i*log(t) - lgamma(i+1) - lgamma(r) - (r+i)*log(alpha+t));
      }
      P2b = P2b * (1-cmf);
    }
  }
  double res = P1 * (P2a + P2b);
  return res;
}

// [[Rcpp::export]]
NumericVector xbgcnbd_exp_cpp(NumericVector params, NumericVector t, bool dropout_at_zero = false) {
  int N = t.size();
  NumericVector res = rep(0.0, N);

  // heuristic: stop infinite sum at 99.9% quantile of NBD; but at least 100
  int k        = params[0];
  double r     = params[1];
  double alpha = params[2];
  int stop;
  double add;
  for (int j=0; j<N; j++) {
    stop = k * R::qnbinom(0.9999, r, alpha/(alpha+t[j]), TRUE, FALSE);
    //Rcpp::Rcout << " stop:" << stop << std::endl;
    if (stop < 100) stop = 100;
    for (int i=1; i<stop; i++) {
      add = i * xbgcnbd_pmf_cpp(params, t[j], i, dropout_at_zero);
      res[j] += add;
      if ((add < 1e-8) & (i >= 100)) {
        //Rcpp::Rcout << " break:" << i << std::endl;
        break;
      }
    }
  }
  return res;
}


/*** R
  params <- c(3, 0.5, 2, 0.3, 0.6)
  stopifnot(round(xbgcnbd_ll_cpp(params, c(3,4), c(12,12), c(14,14), c(3,3)), 5) == c( -10.42500,-12.22396))
  stopifnot(round(xbgcnbd_pmf_cpp(params, 10, 2), 5) == 0.05669)
  stopifnot(round(xbgcnbd_pmf_cpp(params, 10, 2, TRUE), 5) == 0.04548)
  stopifnot(round(sum(sapply(0:1000, function(i) xbgcnbd_pmf_cpp(params, 10, i))), 5) == 1)
  stopifnot(round(sum(sapply(1:1000, function(i) i * xbgcnbd_pmf_cpp(params, 10, i))), 5) ==
            round(xbgcnbd_exp_cpp(params, 10), 5))
*/
