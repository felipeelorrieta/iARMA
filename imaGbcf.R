
imaGbcf <- function(theta, time){
  N <- length(time)

  Delta <- c(NA, diff(x = time))

  c <- 1 + (theta**2)

  for(n in 2:N){
    c[n] <- c[1] - (((sign(theta)*(abs(theta)**Delta[n]))**2) / c[n-1])
  }

  return(list(times = time,
              c = c))
}

# source('/home/cesar/Desktop/iMA/timeSimExp.R')
#
# imaGbcf(theta = -0.8, time = timeSimExp(N = 250, rate1 = 0.1))
