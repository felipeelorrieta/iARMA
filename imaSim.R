
imaSim <- function(theta, var = 1, time, M = 1){
  # trajectories length
  N <- length(x = time)

  # successive time differences t_{2}-t_{1},\ldots,t_{N}-t_{N-1}
  Delta <- diff(x = time)

  # IMA covariance matrix
  Gamma <-  var * Matrix::bandSparse(n = N,
                                     k = -0:1,
                                     diagonals = list(rep(1 + (theta**2), N), sign(theta)*(abs(theta)**Delta)),
                                     symmetric = TRUE)

  # simulation of M trajectories from MN(0, Gamma)
  series <- sparseMVN::rmvn.sparse(n = M,
                                   mu = rep(0, N),
                                   CH = Matrix::Cholesky(A = Gamma),
                                   prec = FALSE)

  # putting labels
  rownames(series) <- rownames(series,
                               do.NULL = FALSE,
                               prefix = "s.")

  # creating zoo object
  return(list(times = time,
              seriesSim = t(series)))
}

# source('/home/cesar/Desktop/iMA/timeSimExp.R')
# 
# set.seed(1234)
# imaSeries <- imaSim(theta = -0.4,
#                     time = timeSimExp(N = 250, rate1 = 0.1),
#                     M = 4)
# 
# s <- zoo::zoo(x = imaSeries$seriesSim,
#               order.by = imaSeries$times)
# 
# plot(s, nc = 2,
#      panel = function(...){
#        lines(...)
#        rug(x = zoo::index(s), col = 2)
#       })
