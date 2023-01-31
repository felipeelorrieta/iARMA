

iarmaMLE <- function(iPar = c(phi = 0.5, theta = 0.5, var = 1), serie, time, ...){
  mleEst <- optim(par = iPar,
                  fn = iarmaMinusLogLik,
                  serie = serie,
                  time = time,
                  lower = c(l.phi = 0.01, l.theta = 0.01, l.var = 0.01),
                  upper = c(u.phi = 0.99, u.theta = 0.99, u.var = Inf),
                  method = 'L-BFGS-B',
                  ...)
  
  return(mleEst)
}

# require(FKF)
# source('iarmaGbcf.R')
# source("iarmaSim.R")
# source("iarmaSS.R")
# source("iarmaMinusLogLik.R")
# source("istsa/R/timeSimExp.R")
# 
# sim <- iarmaSim(phi = 0.8, theta = 0.3, var = 15,
#                 time = timeSimExp(N = 1000),
#                 M = 4)
# 
# fit <- iarmaMLE(serie = zoo::coredata(sim[,1]), time = zoo::index(sim), hessian = TRUE)
