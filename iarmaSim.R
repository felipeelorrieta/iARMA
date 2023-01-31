

iarmaSim <- function(phi, theta, var = 1, time, M = 1){
  # trajectory length
  N <- length(x = time)
  
  # NA, \Delta_{2}, \ldots, \Delta_{N}
  Delta <- c(NA, diff(x = time))
  
  # NA, \phi^{\Delta_{2}}, \ldots, \phi^{\Delta_{N}}
  phiD <- phi**Delta
  # NA, \theta^{\Delta_{2}}, \ldots, \theta^{\Delta_{N}}
  thetaD <- theta**Delta
  
  # general backward continued fraction
  c <- zoo::coredata(iarmaGbcf(phi, theta, time))
  
  #simulated errors
  mzeta <- replicate(n = M, rnorm(N))
  
  # transition equation
  malpha <- apply(X = mzeta, MARGIN = 2, FUN = function(zeta){
    alpha <- 0
    for(n in 2:N)
      alpha[n] <- (phiD[n]*alpha[n-1]) + ((phiD[n]+(thetaD[n]/c[n-1]))*sqrt(var*c[n-1])*zeta[n-1])
    alpha
  })
  
  # simulation of M trajectories
  mserie <- sapply(X = 1:M, FUN = function(m){
    serie <- malpha[,m] + sqrt(var*c)*mzeta[,m]
    serie
  })
  
  # putting labels
  colnames(mserie) <- colnames(mserie, do.NULL = FALSE, prefix = "s.")
  
  # creating zoo object
  return(zoo::zoo(x = mserie, order.by = time))
}

# source('iarmaGbcf.R')
# source("istsa/R/timeSimExp.R")
# sim <- iarmaSim(phi = 0.8, theta = 0.3, var = 15,
#                 time = timeSimExp(N = 100),
#                 M = 4)
# plot(sim)
