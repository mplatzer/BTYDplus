context("mle")

test_that("NBD", {
    
    # generate artificial NBD data
    set.seed(1)
    params <- c(r = 0.85, alpha = 4.45)
    cbs <- nbd.GenerateData(1000, 32, 32, params)$cbs
    
    # estimate parameters, and compare to true parameters
    est <- nbd.EstimateParameters(cbs[, c("x", "T.cal")])
    
    # require less than 5% deviation in estimated parameters
    ape <- function(act, est) abs(act - est)/act
    expect_true(ape(params[1], est[1]) < 0.05)
    expect_true(ape(params[2], est[2]) < 0.05)
    
    # estimate future transactions in holdout-period with true params
    cbs$x.est <- nbd.ConditionalExpectedTransactions(params, cbs$T.star, cbs$x, cbs$T.cal)
    
    # require less than 5% deviation in estimated transactions
    expect_true(ape(sum(cbs$x.star), sum(cbs$x.est)) < 0.05)
    
    expect_true(min(cbs$x.star) >= 0)
    expect_true(all(cbs$x.star == round(cbs$x.star)))
    
}) 
