
cbgnbd.EstimateParameters <- function (cal.cbs, par.start = c(1, 1, 1, 1), max.param.value = 10000)
{
  return(cbgcnbd.EstimateParameters(cal.cbs, k=1, par.start, max.param.value))
}

cbgnbd.cbs.LL <- function (params, cal.cbs) 
{
  return(cbgcnbd.cbs.LL(params, k=1, cal.cbs))
}

cbgnbd.LL <- function (params, x, t.x, T.cal) 
{
  return(cbgcnbd.LL(params, k=1, x, t.x, T.cal))
}
  
cbgnbd.PAlive <- function (params, x, t.x, T.cal)
{
  return(cbgcnbd.PAlive(params, k=1, x, t.x, T.cal))
}
  
cbgnbd.ConditionalExpectedTransactions <- function (params, T.star, x, t.x, T.cal) 
{
  return(cbgcnbd.ConditionalExpectedTransactions(params, k=1, T.star, x, t.x, T.cal))
}

cbgnbd.GenerateData <- function (n, T.cal, T.star, params, return.elog=F) 
{
  return(cbgcnbd.GenerateData(n, k=1, T.cal, T.star, params, return.elog))
}
