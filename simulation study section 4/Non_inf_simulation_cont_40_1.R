#nCONTAMINATED DATA,
#SIMUALTION STUDY WITH N=20 data points per group, N=40 in total

rm(list=ls())


#install.packages("metRology")
library(metRology)
library(MASS)
library(TeachingDemos)


#proposal distributions
proposal=function(R){
  c(runif(R,110,130), (runif(R,-3,9)), runif(R,1,8))
}


# function that generated the data (20 datapoints per group)
ygen=function(n=20, theta1, theta2,sigma ){
  y1=rnorm(n,  theta1,sigma )
  y2=rnorm(n,  theta1+theta2,sigma )
  dat=data.frame(x=Treatment, y=c( y1, y2))
  return(dat=dat)
}

# function that transforms the output of a Simulation-based CD into a confidence density
to_cd=function(values=val[[i]], ndens=500, from=NULL, to=NULL, bw =0.35){
  if (length(values)<=2) {
    values=c(-999,-999)
    return(list(theta=values))
  }
  values=values[!is.na(values)]
  F1=density(values, n=ndens, from=from,bw=bw,  to=to,kernel="rectangular" )
  for(i in 2:length(F1$y)){    if(F1$y[i]<F1$y[i-1]) F1$y[i]=F1$y[i-1]   }
  #plot(F1)
  Fcheck=  (F1$y-min(F1$y))
  Fcheck=Fcheck /max(Fcheck)
  # plot(Fcheck)
  quantiles= seq(0,1, length=round(ndens*0.5))
  qq=sapply(quantiles, function (x) which.min(abs(Fcheck-x)))
  return(list( theta=F1$x[qq]  ))
}

#number of simulations
R=2000
#matricx for Confidence intervals (type hpd and quantile)
CD_hpd95_median=CD_q95_median=CD_hpd90_median=CD_q90_median=
  ABC_hpd95_median=ABC_q95_median=ABC_hpd90_median=ABC_q90_median=matrix(NA, ncol=2,nrow=R)
CD_hpd95_m=CD_q95_m=CD_hpd90_m=CD_q90_m=
  ABC_hpd95_m=ABC_q95_m=ABC_hpd90_m=ABC_q90_m=matrix(NA, ncol=2,nrow=R)
CD_hpd95_ree=CD_q95_ree=CD_hpd90_ree=CD_q90_ree=
  ABC_hpd95_ree=ABC_q95_ree=ABC_hpd90_ree=ABC_q90_ree=matrix(NA, ncol=2,nrow=R)
CD_hpd95_asy_lm=CD_q95_asy_lm=CD_hpd90_asy_lm=CD_q90_asy_lm=matrix(NA, ncol=2,nrow=R)
CD_hpd95_asy_m=CD_q95_asy_m=CD_hpd90_asy_m=CD_q90_asy_m=matrix(NA, ncol=2,nrow=R)

# confidence point estimators
CDmedian_median=NULL
CDmedian_m=NULL
CDmedian_ree=NULL
CDmedian_asy_lm=NULL
CDmedian_asy_m=NULL
ABCmedian_median=NULL
ABCmedian_m=NULL
ABCmedian_ree=NULL

#I type error
theta0=1
ree_theta0=0
m_theta0=0
median_theta0=0
abc_ree_theta0=0
abc_m_theta0=0
abc_median_theta0=0
lmasy_theta0=0
masy_theta0=0
 

############################################################################
############################################################################
# Robust estimating function with 
# y :  simulated data
# beta: coefficient of the model on the real data (parameter estimated)
# lambda is the simulated nuisance parameter (intercept of the model)
# sigma is the simulated nuisance parameter (sd)
REE=function (M=M,newdata=dati, lambda=lambda, sigma=sigma){
  y=newdata$y
  X=model.matrix(M)
  beta=M$coefficients
  beta[1]=lambda
  J <- integrate(f = Vectorize(FUN = function(x) dnorm(x) * 
                sum(MASS:::psi.huber(1/sigma*X[,2]*(x-X%*%beta), k = 1.345)^(2))), 
                lower = -10, upper = 10)$value
  b2=robustBLME:::vpsi_huber(1/sigma*X[,2]*(y-X%*%beta),1.345,n*2)
  sum(b2)*1/sqrt(J)
}

r=1
R=2000
n=20 #for each group
theta0=1 #true value in this scenario
DELTA=0.1 # ABC- appriximation tolerance

############################################################################
############################################################################
#simulations
set.seed(1)
for (r in 1:2000) {

  # generate the data (real data for each simulation)
  y1=rnorm(n, mean =  115.5,4)
  y1=sort(y1,decreasing = F)
  y1[1:4]=115.5- abs(rcauchy(4,location=0,scale = 4)) # add CONTAMINATION
  y2=rnorm(n,mean = 115.5+theta0, 4)
  
  # create a dataset
  Treatment=c(rep("New",n), rep("Standard", n))
  dati=data.frame(x=Treatment, y=c( y1, y2))
  LM=lm(y~x, data=dati )
  
  # generate the proposals and the simulate the fake data given the proposal
  # using the central model
  # the same fake data are used for each method
  gendata=function (R){
    thetai=matrix(ncol=3, nrow=R)
    thetai=t(apply(thetai, 1, function (c) proposal(1)))
    Y=apply(thetai,1,function (x) ygen(n=n,theta1=x[1],theta2=x[2],sigma =x[3] ))
    return(list(Y=Y, thetai=thetai[,2], lambdai=thetai[,1],sigmai=thetai[,3]))
  }
  ris=gendata(R =R)
  Y=ris$Y
  
  ### MEDIAN
  stat_median=sapply(1:length(Y),  function (x) -median(Y[[x]]$y[1:(n)] )+
                       median(Y[[x]]$y[(n+1):(2*n)] ))
  s= median(y2)-median(y1)
  #CD
  accok=  sapply(stat_median,  function (x) (s<=x))
  ris_median=list();   ris_median$theta=ris$thetai[accok]
  #ABC
  accok=  sapply(stat_median,  function (x)  abs(s-x)<=DELTA)
  ris_ABCmedian=list();   ris_ABCmedian$theta=ris$thetai[accok]
  
  
  ### Robust estimating equation
  statREE=function(M=M,lambda=lambda,sigma=sigma, dati=dati){
    REE(newdata = dati,M = M,lambda=lambda, sigma = sigma)
  }
  M=rlm(y~x, data=dati )
  stat_ree=sapply(1:length(Y),  function (x)  statREE(M=M,ris$Y[[x]],
                                        lambda = ris$lambdai[x], sigma = ris$sigmai[x]))
  s=sapply(1:1, function (x) statREE(M=M,dati= dati,
                                 lambda = M$coefficients[1], sigma = M$s ))
  # CD
  accok=  sapply(1:length(stat_ree),  function (x) (s <=stat_ree[x]))
  ris_ree=list(); ris_ree$theta=ris$thetai[accok]
  #ABC
  accok=  sapply(1:length(stat_ree),  function (x)  abs(s -stat_ree[x])<=DELTA)
  ris_ABCree=list();   ris_ABCree$theta=ris$thetai[accok]
  
  # M estimator
  Mest=function (X) { 
    dat=data.frame(x=Treatment, y=c(X$y))
    MM=rlm(y ~x, data=dat) 
    s=MM$coefficients[2]
    return(s)
  }
  s=M$coefficients[2]
  stat_M=sapply(1:length(Y),  function(x) Mest(Y[[x]]))
  accok=  sapply(stat_M,  function (x) (s<=x))
  ris_M=list();   ris_M$theta=ris$thetai[accok]
  accok=  sapply(stat_M,  function (x)  abs(s-x)<=DELTA)
  ris_ABCM=list();   ris_ABCM$theta=ris$thetai[accok]
  
  # asymptotic WALD
  est_LM=LM$coefficients[2]; sd_est_LM=summary(LM)$coefficients[2,2]
  # asymptotic WALD/ROBUST M EST
  est_M=M$coefficients[2]; sd_est_M=summary(M)$coefficients[2,2]
  
  asy_LM=rt.scaled(10000,(38),mean = est_LM, sd_est_LM)
  asy_M=rnorm(10000,mean = est_M, sd_est_M)
  
  
  # CONVERT CDistributions to cdensity
  ris_cd_median=to_cd(ris_median$theta,ndens = R, from=-3,to=8, bw=0.2)
  ris_cd_M=to_cd(ris_M$theta,ndens = R, from=-3,to=8, bw=0.2)
  ris_cd_ree=to_cd(ris_ree$theta,ndens = R, from=-3,to=8, bw=0.2)
  
  par(mfrow=c(2,2))
  
  
  # store the results
  # median
  CD_hpd95_median[r,]=emp.hpd(ris_cd_median$theta,conf = 0.95)
  CD_q95_median[r,]=c(quantile(ris_cd_median$theta,0.025),
                      quantile(ris_cd_median$theta,0.975))
  CD_hpd90_median[r,]=emp.hpd(ris_cd_median$theta,conf = 0.90)
  CD_q90_median[r,]=c(quantile(ris_cd_median$theta,0.05),
                      quantile(ris_cd_median$theta,0.95))
  ABC_hpd95_median[r,]=emp.hpd(ris_ABCmedian$theta,conf = 0.95)
  ABC_q95_median[r,]=c(quantile(ris_ABCmedian$theta,0.025),
                       quantile(ris_ABCmedian$theta,0.975))
  ABC_hpd90_median[r,]=emp.hpd(ris_ABCmedian$theta,conf = 0.90)
  ABC_q90_median[r,]=c(quantile(ris_ABCmedian$theta,0.05),
                       quantile(ris_ABCmedian$theta,0.95))
  
  median_theta0=median_theta0+as.numeric(quantile(ris_cd_median$theta,0.95)>theta0)
  ree_theta0=ree_theta0+as.numeric(quantile(ris_cd_ree$theta,0.95)>theta0)
  m_theta0=m_theta0+as.numeric(quantile(ris_cd_M$theta,0.95)>theta0)
  abc_median_theta0=median_theta0+as.numeric(quantile(ris_cd_median$theta,0.95)>theta0)
  abc_ree_theta0=ree_theta0+as.numeric(quantile(ris_cd_ree$theta,0.95)>theta0)
  abc_m_theta0=m_theta0+as.numeric(quantile(ris_cd_M$theta,0.95)>theta0)
  masy_theta0=masy_theta0+as.numeric(quantile(asy_M,0.95)>theta0)
  lmasy_theta0=lmasy_theta0+as.numeric(metRology::qt.scaled(0.95,df = 38,mean = est_LM ,sd=sd_est_LM)>theta0)
  
  # M estimators
  CD_hpd95_m[r,]=emp.hpd(ris_cd_M$theta,conf = 0.95)
  CD_q95_m[r,]=c(quantile(ris_cd_M$theta,0.025),
                 quantile(ris_cd_M$theta,0.975))
  CD_hpd90_m[r,]=emp.hpd(ris_cd_M$theta,conf = 0.90)
  CD_q90_m[r,]=c(quantile(ris_cd_M$theta,0.05),
                 quantile(ris_cd_M$theta,0.95))
  ABC_hpd95_m[r,]=emp.hpd(ris_ABCM$theta,conf = 0.95)
  ABC_q95_m[r,]=c(quantile(ris_ABCM$theta,0.025),
                  quantile(ris_ABCM$theta,0.975))
  ABC_hpd90_m[r,]=emp.hpd(ris_ABCM$theta,conf = 0.90)
  ABC_q90_m[r,]=c(quantile(ris_ABCM$theta,0.05),
                  quantile(ris_ABCM$theta,0.95))
  #Robust estimating equations
  ABC_hpd95_ree[r,]=emp.hpd(ris_ABCree$theta,conf = 0.95)
  ABC_q95_ree[r,]=c(quantile(ris_ABCree$theta,0.025 ),
                    quantile(ris_ABCree$theta,0.975 ))
  ABC_hpd90_ree[r,]=emp.hpd(ris_ABCree$theta,conf = 0.90)
  ABC_q90_ree[r,]=c(quantile(ris_ABCree$theta,0.05 ),
                    quantile(ris_ABCree$theta,0.95))
  
  CD_hpd95_ree[r,]=emp.hpd(ris_cd_ree$theta,conf = 0.95)
  CD_q95_ree[r,]=c(quantile(ris_cd_ree$theta,0.025),
                   quantile(ris_cd_ree$theta,0.975))
  CD_hpd90_ree[r,]=emp.hpd(ris_cd_ree$theta,conf = 0.90)
  CD_q90_ree[r,]=c(quantile(ris_cd_ree$theta,0.05),
                   quantile(ris_cd_ree$theta,0.95))
  
  CD_hpd95_asy_lm[r,]=emp.hpd(asy_LM,conf = 0.95)
  CD_q95_asy_lm[r,]=c(quantile(asy_LM,0.025),
                      quantile(asy_LM,0.975))
  CD_hpd90_asy_lm[r,]=emp.hpd(asy_LM,conf = 0.90)
  CD_q90_asy_lm[r,]=c(quantile(asy_LM,0.05),
                      quantile(asy_LM ,0.95))
  
  CD_hpd95_asy_m[r,]=emp.hpd(asy_M,conf = 0.95)
  CD_q95_asy_m[r,]=c(quantile(asy_M,0.025),
                     quantile(asy_M,0.975))
  CD_hpd90_asy_m[r,]=emp.hpd(asy_M,conf = 0.90)
  CD_q90_asy_m[r,]=c(quantile(asy_M,0.05),
                     quantile(asy_M ,0.95))
  
  CDmedian_median[r]=median(ris_cd_median$theta)
  CDmedian_m[r]=median(ris_cd_M$theta)
  CDmedian_ree[r]=median(ris_cd_ree$theta)
  
  CDmedian_asy_lm[r]=median(asy_LM)
  CDmedian_asy_m[r]=median(asy_M)
  ABCmedian_median[r]=median(ris_ABCmedian$theta)
  ABCmedian_m[r]=median(ris_ABCM$theta)
  ABCmedian_ree[r]=median(ris_ABCree$theta)
  
  cat(r);  
  print(c("cd median ",
          mean(apply(cbind(CD_hpd95_median,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
          mean(apply(cbind(CD_hpd90_median,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
          mean(apply(cbind(CD_q95_median,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
          mean(apply(cbind(CD_q90_median,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T)))
  
  print(c("abc median ",
          mean(apply(cbind(ABC_hpd95_median,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
          mean(apply(cbind(ABC_hpd90_median,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
          mean(apply(cbind(ABC_q95_median,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
          mean(apply(cbind(ABC_q90_median,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T)))
  
  
  print(c("cd m",
          mean(apply(cbind(CD_hpd95_m,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
          mean(apply(cbind(CD_hpd90_m,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
          mean(apply(cbind(CD_q95_m,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
          mean(apply(cbind(CD_q90_m,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T)))
  
  print(c("abc m",
          mean(apply(cbind(ABC_hpd95_m,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
          mean(apply(cbind(ABC_hpd90_m,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
          mean(apply(cbind(ABC_q95_m,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
          mean(apply(cbind(ABC_q90_m,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T)))
  
  
  print(c("cd ree",
          mean(apply(cbind(CD_hpd95_ree,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
          mean(apply(cbind(CD_hpd90_ree,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
          mean(apply(cbind(CD_q95_ree,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
          mean(apply(cbind(CD_q90_ree,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T)))
  
  print(c("abc ree",
          
          mean(apply(cbind(ABC_hpd95_ree,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
          mean(apply(cbind(ABC_hpd90_ree,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
          mean(apply(cbind(ABC_q95_ree,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
          mean(apply(cbind(ABC_q90_ree,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T)))
  cat(r); cat("\n")
  cat((ree_theta0)/r)
  cat((lmasy_theta0)/r)
  save.image("cont_40_1.RData")
}
 