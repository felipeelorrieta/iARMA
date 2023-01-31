
# create a state space representation out of the IARMA parameters

# we suppose series is centered in the argument

iarmaSS <- function(phi, theta, var, serie, time){
  # trajectory length
  N <- length(x = time)
  
  # \Delta_{2}, \ldots, \Delta_{N}
  Delta <- diff(x = time)
  
  # \phi^{\Delta_{2}}, \ldots, \phi^{\Delta_{N}}
  phiD <- phi**Delta
  # \theta^{\Delta_{2}}, \ldots, \theta^{\Delta_{N}}
  thetaD <- theta**Delta
  
  # general backward continued fraction
  c <- zoo::coredata(iarmaGbcf(phi, theta, time))
  
  # state space representation
  a0 <- 0
  P0 <- matrix(0)
  dn <- matrix((c(phiD, 0) + (c(thetaD, 0)/c)) * serie, nrow = 1)
  Tn <- array(-(c(thetaD, 0)/c), dim = c(1, 1, N))
  HHn <- array(0, dim = c(1, 1, 1))
  cn <- matrix(0)
  Zn <- array(1, dim = c(1, 1, 1))
  GGn <- array(var * c, dim = c(1, 1, N))
  
  return(list(a0 = a0, P0 = P0, dt = dn, ct = cn, Tt = Tn, Zt = Zn, HHt = HHn, GGt = GGn))
}

# source('iarmaGbcf.R')
# source("iarmaSim.R")
# source("istsa/R/timeSimExp.R")
# 
# sim <- iarmaSim(phi = 0.8, theta = 0.3, var = 15,
#                 time = timeSimExp(N = 100),
#                 M = 4)
# 
# serie <- zoo::coredata(sim[,1])
# time <- zoo::index(sim)
# par <- c(0.2, 0.6, 30)
# 
# iarmaSS(phi = par[1], theta = par[2], var = par[3], serie, time)
