

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


# source("istsa/R/timeSimExp.R")
# time <- timeSimExp(N = 100)
# summary(iarmaGbcf(phi = 0.2, theta = 0.6, time))

