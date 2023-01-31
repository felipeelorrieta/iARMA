

iarmaBootSamples <- function(phi, theta, var = 1, serie, time, B = 500){
  N <- length(serie)
  
  # NA, \Delta_{2}, \ldots, \Delta_{N}
  Delta <- c(NA, diff(x = time))
  
  # NA, \phi^{\Delta_{2}}, \ldots, \phi^{\Delta_{N}}
  phiD <- phi**Delta
  
  # NA, \theta^{\Delta_{2}}, \ldots, \theta^{\Delta_{N}}
  thetaD <- theta**Delta
  
  # general backward continued fraction
  c <- zoo::coredata(iarmaGbcf(phi, theta, time))
  
  predV <- zoo::coredata(iarmaPredictOneStep(phi, theta, var, serie, time))
  
  resSt <- ((serie - predV[, 'PredictOneStep']) / sqrt(predV[, 'mspeOneStep']))[-1]
  
  resStCen <- resSt - mean(resSt)
  
  bootSam <- replicate(n = B,
                       expr = {
                         bRes <- sample(resStCen,
                                        size = N,
                                        replace = TRUE)
                         
                         bSam <- sqrt(unname(predV[1, 'mspeOneStep'])) * bRes[1]
                         
                         for(n in 2:N){
                           bSam[n] <- (phiD[n] * bSam[n-1]) + 
                             (sqrt(unname(predV[n, 'mspeOneStep'])) * bRes[n]) +
                             ((thetaD[n] / c[n-1]) * (sqrt(unname(predV[n-1, 'mspeOneStep'])) * bRes[n-1]))
                         }
                         bSam
                         }
                       )
  
  colnames(bootSam) <- colnames(bootSam,
                                do.NULL = FALSE,
                                prefix = "b.s.")
  
  return(zoo::zoo(x = bootSam + mean(serie), order.by = time))#
}

# require(FKF)
# source('iarmaGbcf.R')
# source("iarmaSim.R")
# source("iarmaSS.R")
# source("iarmaMinusLogLik.R")
# source("iarmaMLE.R")
# source("iarmaPredictOneStep.R")
# source("istsa/R/timeSimExp.R")
# 
# sim <- iarmaSim(phi = 0.8, theta = 0.3, var = 15,
#                 time = timeSimExp(N = 1000),
#                 M = 4)
# 
# fit <- iarmaMLE(serie = zoo::coredata(sim[,1]), time = zoo::index(sim), hessian = TRUE)
# 
# bootSam <- iarmaBootSamples(phi = fit$par["phi"], theta = fit$par["theta"], var = fit$par["var"],
#                             serie = zoo::coredata(sim[,1]), time = zoo::index(sim), B = 5)
# 
# bootEst <- apply(X = zoo::coredata(bootSam), MARGIN = 2,
#                  FUN = function(s){
#                    iarmaMLE(serie = s,
#                             time = zoo::index(sim))$par
#                    })
# 
# rowMeans(bootEst)
