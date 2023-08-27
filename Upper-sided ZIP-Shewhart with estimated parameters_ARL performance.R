### ARL and SDRL performance of the upper-sided ZIP-Shewhart chart with estimated parameters
### Here, we use the approximation for the W random variable
require(VGAM)
#################
ZIP.performance<-function(m,theta,lam0,K){
# m is the size of the preliminary sample
# theta is equal to 1-phi
# lam0 is the in-control lambda of the Poisson distribution
# K is the design parameter of the upper-sided ZIP-Shewhart chart
#######
# Evaluation of the first three moments
  m1<-m*lam0*theta
  m2<-m*lam0*(-theta)*(-1+m*lam0*(-theta)-lam0*(1-theta))
  m3<-m*lam0*(lam0+(1+lam0)^2-3*(-1+m)*lam0*(1+lam0)*(-theta)+(-2+m)*(-1+m)*(lam0^2)*(-theta)^2)*theta
# Evaluation of the 2nd and the 3rd central moment
  mu2<-m2-(m1^2)
  mu3<-m3-3*m1*m2+2*(m1^3)
# Definition of the alpha, beta and gamma, the parameters of the shifted gamma distribution
  alpha1<-(4*(mu2^3))/(mu3^2)
  beta1<-mu3/(2*mu2)
  gamma1<-m1-2*(mu2^2)/(mu3)
  alpha1
  beta1
  gamma1
# CDF and PMF of the shifted gamma distribution
  F1<-function(x){pgamma(x-gamma1,shape=alpha1,scale=beta1,log=FALSE)}
  f1<-function(x){F1(x+0.5)-F1(x-0.5)}
  xmin<-max(0,floor(qgamma(10^(-15),shape=alpha1,scale=beta1,lower.tail=TRUE,log.p=FALSE)+gamma1))
  xmax<-ceiling(qgamma(1-10^(-15),shape=alpha1,scale=beta1,lower.tail=TRUE,log.p=FALSE)+gamma1)
  ###################################################################################################
# Define the UCL
  UCL0<-floor(lam0*theta+K*sqrt(lam0*theta*(1+lam0*(1-theta))))
# The OoC probability in the known parameters case
  beta0<-1-pzipois(UCL0,lambda=lam0,pstr0=(1-theta))
# In-control ARL and SDRL in the known parameters case
  ARL0<-1/(beta0)
  SDRL0<-sqrt(1-beta0)/(beta0)
  ARL0
  SDRL0
  ###################################################################################################
# Evaluation of the unconditional ARL and SDRL
  y<-seq(xmin,xmax,by=1)
  UCL1<-function(x){
    floor((x/(m*theta))*theta+K*sqrt((x/(m*theta))*theta*(1+(x/(m*theta))*(1-theta))))
  }
  Cbeta<-function(x){1-pzipois(UCL1(x),lambda=lam0,pstr0=(1-theta))}
  CARL<-function(x){
    1/Cbeta(x)
  }
  z1<-CARL(y)%*%f1(y)
  UARL<-sum(z1)
  z2<-((2-Cbeta(y))/(Cbeta(y)^2))%*%f1(y)
  CRL2<-sum(z2)
  USDRL<-sqrt(CRL2-UARL^2)
# evaluation of the relative difference between the IC ARL values, in the known and the estimated parameters case
  diff<-abs(sum(z1)-ARL0)/ARL0
cat("m:",m," theta:",theta,"  lam0:",lam0,"  K:",K,"  UARL:",UARL," USDRL:",USDRL," ARL0:",ARL0,"  SDRL0:",SDRL0," diff:",diff)}
###############################
## An example, m=500, phi=0.9, lambda0=1.0, K=6.66 ----> 
ZIP.performance(m=500,theta=0.1,lam0=1,K=6.66)
################# END #############################################################

