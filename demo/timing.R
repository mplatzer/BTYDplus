require(data.table)

# generate Erlang-3 with various lambdas
nr_of_customers <- 100
elog <- rbindlist(lapply(1:nr_of_customers, function(i) {
  k <- 3
  lambda <- exp(rnorm(1))
  data.table(cust = i, t = cumsum(rgamma(20, k, k * lambda)))
}))

# estimate regularity parameter k
estimateRegularity(elog, plot = TRUE, method = "wheat")
estimateRegularity(elog, plot = TRUE, method = "mle-minka")
estimateRegularity(elog, plot = TRUE, method = "mle-thom")
estimateRegularity(elog, plot = TRUE, method = "cv") 
