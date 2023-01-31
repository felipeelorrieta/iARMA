
imaMLE <- function(iPar = c(theta = 0.5, var = 1), serie, time, ...){
  mleEst <- optim(par = iPar,
                  fn = imaMinusLogLik,
                  serie = serie,
                  time = time,
                  lower = c(l.theta = -0.99, l.var = 0.01),
                  upper = c(u.theta = 0.99, u.var = Inf),
                  method = 'L-BFGS-B',
                  ...)

  # an optim object
  return(mleEst)
}

# source('/home/cesar/Desktop/iMA/timeSimExp.R')
# source('/home/cesar/Desktop/iMA/imaSim.R')
# source('/home/cesar/Desktop/iMA/imaMinusLogLik.R')
# 
# set.seed(1234)
# imaSeries <- imaSim(theta = -0.9, time = timeSimExp(N = 250, rate1 = 0.1))
# 
# mleEst <- imaMLE(serie = imaSeries$seriesSim,
#                  time = imaSeries$times)
# 
# # ml estimation
# mleEst$par
# 
# # estimated standard errors
# sqrt(diag(solve(optimHess(mleEst$par,
#                           fn = imaMinusLogLik,
#                           serie = imaSeries$seriesSim,
#                           time = imaSeries$times))))
