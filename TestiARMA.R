iarmaMLE0 <- function(iPar = c(phi = 0.5, var = 1), serie, time, ...){
  mleEst <- optim(par = iPar,
                  fn = iarmaMinusLogLik0,
                  serie = serie,
                  time = time,
                  lower = c(l.phi = 0.01, l.var = 0.01),
                  upper = c(u.phi = 0.99, u.var = Inf),
                  method = 'L-BFGS-B',
                  ...)

  return(mleEst)
}

iarmaMinusLogLik0 <- function(par, serie, time){
  # Centered serie
  cSerie <- serie - mean(serie)#

  sp <- iarmaSS(phi = par[1], theta = 0, var = par[2], cSerie, time)#

  ans <- FKF::fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt,
                  ct = sp$ct, Tt = sp$Tt, Zt = sp$Zt,
                  HHt = sp$HHt, GGt = sp$GGt,
                  matrix(cSerie, nrow = 1), check.input = TRUE)#

  return(-ans$logLik)
}


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

iarmaGbcf <- function(phi, theta, time){
  N <- length(time)

  Delta <- c(NA, diff(x = time))

  phiD <- phi**Delta
  thetaD <- theta**Delta

  c <- (1+(2*phi*theta)+(theta**2))/(1-(phi**2))

  for(n in 2:N) c[n] <- (c[1]*(1-(phiD[n]**2))) - (2*phiD[n]*thetaD[n]) - ((thetaD[n]**2)/c[n-1])

  return(zoo::zoo(x = c,
                  order.by = time))
}

#setwd('/data/felipe/ComplexAR/CIAR2p_2')
#source('ComplexIAR2p_Aux.R')
#setwd("/data/felipe/felipe/IrregularModels")
#source('FunctionsIATS.R')

#set.seed(1234)
#st<-gentime(n=300,lambda1=15,lambda2=2)
#newt=st/min(diff(st))
#newt=1:300
#x=CIAR.sample(n=300,phi.R=0.5,phi.I=0,sT=newt)
#x=IAR.sample(n=300,phi=0.9,sT=newt)
#newt=x$t/min(diff(x$t))
#y1=x$y/sqrt(var(x$y))
#newt=x$times
#y1=x$series
#fit=iarmaMLE0(serie = y1,time = newt,hessian=TRUE)
#fit2=CIAR.kalman(y=y1,t=newt)
#fit3=IAR.loglik(y=y1,sT=newt)