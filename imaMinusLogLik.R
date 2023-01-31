
imaMinusLogLik <- function(par, serie, time){
  # serie length
  N <- length(x = serie)

  # successive time differences t_{2}-t_{1},\ldots,t_{N}-t_{N-1}
  Delta <- diff(x = time)

  # IMA covariance matrix
  Gamma <-  par[2] * Matrix::bandSparse(n = N,
                                      k = -0:1,
                                      diagonals = list(rep(1 + (par[1]**2), N), sign(par[1])*(abs(par[1])**Delta)),
                                      symmetric = TRUE)

  # Centered serie
  cSerie <- serie - mean(serie)

  # -log(likelihood) where likelihood = MN(0, Gamma)
  mLogLik <- -sparseMVN::dmvn.sparse(x = cSerie,
                                     mu = rep(0, N),
                                     CH =  Matrix::Cholesky(A = Gamma),
                                     log = TRUE,
                                     prec = FALSE)

  return(mLogLik)
}

# source('/home/cesar/Desktop/iMA/timeSimExp.R')
# source('/home/cesar/Desktop/iMA/imaSim.R')
# 
# set.seed(1234)
# imaSeries <- imaSim(theta = -0.4, time = timeSimExp(N = 250, rate1 = 0.1))
# 
# imaMinusLogLik(par = c(theta = 0.2, var = 2),
#                serie = imaSeries$seriesSim,
#                time = imaSeries$times)
