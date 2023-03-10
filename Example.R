source('imaMinusLogLik.R')
source('imaMLE.R')
source('imaGbcf.R')
source('imaPredictOneStep.R')
source('imaSim.R')
source('iarmaMinusLogLik.R')
source('iarmaMLE.R')
source('iarmaMLE2.R')
source('iarmaGbcf.R')
source('iarmaPredictOneStep.R')
source('iarmaBootSamples.R')
source('iarmaSS.R')
source('TestiARMA.R')
require(iAR)
st<-gentime(n=300)
newt=st/min(diff(st)) #time gaps should be greater than one
y<-IARsample(phi=0.9,st=newt, n=300)
iarma<-iarmaMLE(serie = y$series,time = newt, hessian = T)

y<-arima.sim(n=300,model=list(ar=0.9))
iarma<-iarmaMLE(serie = y,time = 1:300, hessian = T)