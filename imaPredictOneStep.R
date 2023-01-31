
imaPredictOneStep <- function(theta, var = 1, serie, time){
  N <- length(serie)

  cSerie <- serie - mean(serie)

  c <- imaGbcf(theta, time)$c

  Delta <- c(NA, diff(x = time))

  p1s <- 0

  for(n in 2:N){
    p1s[n] <- ((sign(theta)*(abs(theta)**Delta[n])) / c[n-1]) * (cSerie[n-1] - p1s[n-1])
  }
  
  return(list(times = time,
              pred = cbind(serie = serie,
                           PredictOneStep = p1s + mean(serie),
                           mspeOneStep = var * c)))
}

# load('/home/cesar/Dropbox/PhDThesis/ThesisR/data-cts/asth.rda')
# source('/home/cesar/Desktop/iMA/imaMinusLogLik.R')
# source('/home/cesar/Desktop/iMA/imaMLE.R')
# 
# mleEst <- imaMLE(serie = asth[, 2], time = asth[, 1])
# 
# predV <- imaPredictOneStep(theta = unname(mleEst$par['theta']),
#                            var   = unname(mleEst$par['var']),
#                            serie = asth[, 2],
#                            time  = asth[, 1])
# 
# s <- zoo::zoo(x = predV$pred,
#               order.by = predV$times)
# 
# plot(s$serie,
#      panel = function(...){
#        lines(...)
#        rug(x = zoo::index(s), col = 2)
#       }, xlab = "Time in hours", ylab = "", lwd = 1.5, col = "gray40")
# 
# lines(s$PredictOneStep, col = 4, lty = 2, lwd = 2)
# lines(s$PredictOneStep - (1.96 * sqrt(s$mspeOneStep)),
#       col = 4, lty = 3)
# lines(s$PredictOneStep + (1.96 * sqrt(s$mspeOneStep)),
#       col = 4, lty = 3)
