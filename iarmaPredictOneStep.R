

iarmaPredictOneStep <- function(phi, theta, var = 1, serie, time){
  # Centered serie
  cSerie <- serie - mean(serie)#
  
  sp <- iarmaSS(phi, theta, var, cSerie, time)#
  
  ans <- FKF::fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
                  Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt,
                  yt = matrix(cSerie, nrow = 1), check.input = TRUE)#
  
  return(zoo::zoo(x = cbind(serie = serie,
                            PredictOneStep = ans$at[-(length(serie)+1)] + mean(serie),#
                            mspeOneStep = ans$Ft),
                  order.by = time))
}

# require(FKF)
# source('iarmaGbcf.R')
# source("iarmaSim.R")
# source("iarmaSS.R")
# source("iarmaMinusLogLik.R")
# source("iarmaMLE.R")
# source("istsa/R/timeSimExp.R")
# 
# sim <- iarmaSim(phi = 0.8, theta = 0.3, var = 15,
#                 time = timeSimExp(N = 1000),
#                 M = 4)
# 
# fit <- iarmaMLE(serie = zoo::coredata(sim[,1]), time = zoo::index(sim), hessian = TRUE)
# 
# predV <- iarmaPredictOneStep(phi = fit$par["phi"], theta = fit$par["theta"],
#                              var = fit$par["var"],
#                              serie = zoo::coredata(sim[,1]),
#                              time = zoo::index(sim))
# 
# plot(predV$serie, xlab = "Time in hours",
#      ylab = "", lwd = 1.5, col = "gray40")
# 
# rug(zoo::index(predV), col = "red")
# 
# lines(predV$PredictOneStep, col = 4, lty = 2, lwd = 2)
# lines(predV$PredictOneStep - (1.96 * sqrt(predV$mspeOneStep)),
#       col = 4, lty = 3)
# lines(predV$PredictOneStep + (1.96 * sqrt(predV$mspeOneStep)),
#       col = 4, lty = 3)
