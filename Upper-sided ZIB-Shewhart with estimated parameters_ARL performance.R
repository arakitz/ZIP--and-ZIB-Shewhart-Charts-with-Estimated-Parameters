########## The upper-sided ZIB-Shewhart chart with estimated parameters
require(VGAM)
#######################
designmzib<-function(m,theta,p0,n,K){
# m is the size of the preliminary sample
# theta is equal to 1-phi
# p0 is the in-control proportion
# n is the sample size
# K is the design parameter of the upper-sided ZIB-Shewhart chart
###############
# define the first three moments
  m1<-m*n*p0*theta
  m2<-m*n*p0*theta+m*(-1+n)*n*(p0^2)*theta+(-1+m)*m*(n^2)*(p0^2)*(theta^2)
  m3<-m*n*p0*theta*(1+p0*(-3+2*p0-3*n*(-1+p0)*(1+(-1+m)*theta)+(n^2)*p0*(1+(-1+m)*theta*(3+(-2+m)*theta))))
  m4<-m*n*p0*theta*(1+p0*(-7-6*(-2+p0)*p0+n*(-1+p0)*(-7+11*p0)*(1+(-1+m)*theta)-6*(n^2)*(-1+p0)*p0*(1+(-1+m)*theta*(3+(-2+m)*theta))+(n^3)*(p0^2)*(1+(-1+m)*theta*(7+(-2+m)*theta*(6+(-3+m)*theta)))))
 # m1
 # m2
 # m3
 # m4
# define the 2nd and the 3rd central moments 
 mu2<-m2-(m1^2)
  mu3<-m3-3*m1*m2+2*(m1^3)
# define the parameters alpha, beta and gamma of the shifted gamma distribution
  alpha1<-(4*(mu2^3))/(mu3^2)
  beta1<-mu3/(2*mu2)
  gamma1<-m1-2*(mu2^2)/(mu3)
  #alpha1
  #beta1
  #gamma1
# CDF and PMF of the shifted gamma distribution 
 F1<-function(x){pgamma(x-gamma1,shape=alpha1,scale=beta1,log=FALSE)}
  f1<-function(x){F1(x+0.5)-F1(x-0.5)}
  xmin<-max(0,floor(qgamma(0.000001,shape=alpha1,scale=beta1,lower.tail=TRUE,log.p=FALSE)+gamma1))
  xmax<-ceiling(qgamma(0.999999,shape=alpha1,scale=beta1,lower.tail=TRUE,log.p=FALSE)+gamma1)
  ###################################################################################################
# UCL control limit of the upper-sided ZIB-Shewhart chart in the known parameters case
  UCL0<-floor(n*p0*theta+K*sqrt(n*p0*theta*(1-p0+n*p0*(1-theta))))
# the OoC probability
  beta0<-1-pzibinom(UCL0,size=n,prob=p0,pstr0=(1-theta),lower.tail=TRUE,log.p=FALSE)
# IC ARL and SDRL
  ARL0<-1/(beta0)
  SDRL0<-sqrt(1-beta0)/(beta0)
  #ARL0
  #SDRL0
  ###################################################################################################
# Evaluation of the unconditional ARL and SDRL
  y<-seq(xmin,xmax,by=1)
  UCL1<-function(x){
    floor(n*(x/(n*m*theta))*theta+K*sqrt(n*(x/(n*m*theta))*theta*(1-(x/(n*m*theta))+n*(x/(n*m*theta))*(1-theta))))
  }
  Cbeta<-function(x){1-pzibinom(UCL1(x),size=n,prob=p0,pstr0=(1-theta),lower.tail=TRUE,log.p=FALSE)}
  CARL<-function(x){
    1/Cbeta(x)
  }
  z1<-CARL(y)%*%f1(y)
  UARL<-sum(z1)
  z2<-((2-Cbeta(y))/(Cbeta(y)^2))%*%f1(y)
  CRL2<-sum(z2)
  USDRL<-sqrt(CRL2-UARL^2)
# Relative difference between the ARL in the known and the estimated parameters case
  diffs<-abs(sum(z1)-ARL0)/ARL0
  cat("theta",theta,"   n:",n,"  p0:",p0,"xmin:",xmin,"  xmax:",xmax,"  m:",m,"  diffs",diffs,"  K:",K,"  UARL",UARL,"  USDRL",USDRL,"  ARL0",ARL0,"  SDRL0",SDRL0,"\n")}
####################
## An example, m=1000, phi=0.8, p0=0.002, n=200, K=5.92 --->
designmzib(m=1000,theta=0.2,p0=0.002,n=200,K=5.92)
##################### END ##############################################################