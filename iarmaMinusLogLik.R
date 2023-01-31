

iarmaMinusLogLik <- function(par, serie, time){
  # Centered serie
  cSerie <- serie - mean(serie)#
  
  sp <- iarmaSS(phi = par[1], theta = par[2], var = par[3], cSerie, time)#
  
  ans <- FKF::fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt,
                  ct = sp$ct, Tt = sp$Tt, Zt = sp$Zt,
                  HHt = sp$HHt, GGt = sp$GGt,
                  matrix(cSerie, nrow = 1), check.input = TRUE)#
  
  return(-ans$logLik)
}

# require(FKF)
# source('iarmaGbcf.R')
# source("iarmaSim.R")
# source("iarmaSS.R")
# source("istsa/R/timeSimExp.R")
# 
# sim <- iarmaSim(phi = 0.8, theta = 0.3, var = 15,
#                 time = timeSimExp(N = 100),
#                 M = 4)
# 
# iarmaMinusLogLik(par = c(0.8, 0.3, 15),
#                  serie = zoo::coredata(sim[,1]),
#                  time = zoo::index(sim))
