#NONINF
#SIMULATION STUDY 10 CONT
#load("cont_80.RData")

library(MASS)
library(TeachingDemos)

proposal=function(R){
  c(runif(R,110,130), (runif(R,-3,9)), runif(R,1,8))
}


ygen=function(n=40, theta1, theta2,sigma ){
  y1=rnorm(n,  theta1,sigma )
  y2=rnorm(n,  theta1+theta2,sigma )
  dat=data.frame(x=Treatment, y=c( y1, y2))
  return(dat=dat)
}


to_cd=function(values=val[[i]], ndens=500, from=NULL, to=NULL, bw=0.35 ){
  if (length(values)<=2) {
    values=c(-999,-999)
    return(list(theta=values))
  }
  values=values[!is.na(values)]
  F1=density(values, n=ndens, from=from,bw=bw,  to=to)
  for(i in 2:length(F1$y)){    if(F1$y[i]<F1$y[i-1]) F1$y[i]=F1$y[i-1]   }
  plot(F1)
  Fcheck=  (F1$y-min(F1$y))
  Fcheck=Fcheck /max(Fcheck)
  # plot(Fcheck)
  quantiles= seq(0,1, length=round(ndens*0.5))
  qq=sapply(quantiles, function (x) which.min(abs(Fcheck-x)))
  
  theta=F1$x[qq]
  #plot(density(theta))
  return(list( theta=theta  ))
}
 
R=2000

CD_hpd95_median=CD_q95_median=CD_hpd90_median=CD_q90_median=
  ABC_hpd95_median=ABC_q95_median=ABC_hpd90_median=ABC_q90_median=matrix(NA, ncol=2,nrow=R)
CD_hpd95_m=CD_q95_m=CD_hpd90_m=CD_q90_m=
  ABC_hpd95_m=ABC_q95_m=ABC_hpd90_m=ABC_q90_m=matrix(NA, ncol=2,nrow=R)
CD_hpd95_ree=CD_q95_ree=CD_hpd90_ree=CD_q90_ree=
  ABC_hpd95_ree=ABC_q95_ree=ABC_hpd90_ree=ABC_q90_ree=matrix(NA, ncol=2,nrow=R)
CD_hpd95_asy_lm=CD_q95_asy_lm=CD_hpd90_asy_lm=CD_q90_asy_lm=matrix(NA, ncol=2,nrow=R)
CD_hpd95_asy_m=CD_q95_asy_m=CD_hpd90_asy_m=CD_q90_asy_m=matrix(NA, ncol=2,nrow=R)


CDmedian_median=NULL
CDmedian_m=NULL
CDmedian_ree=NULL
CDmedian_asy_lm=NULL
CDmedian_asy_m=NULL
ABCmedian_median=NULL
ABCmedian_m=NULL
ABCmedian_ree=NULL


theta0=1
ree_margintheta0=0
m_margintheta0=0
median_margintheta0=0
abc_ree_margintheta0=0
abc_m_margintheta0=0
abc_median_margintheta0=0
lmasy_margintheta0=0
masy_margintheta0=0
 

REE=function (M=M,newdata=dati, lambda=lambda, sigma=sigma){
  y=newdata$y
  X=model.matrix(M)
  beta=M$coefficients
  beta[1]=lambda
  J  <- integrate(f = Vectorize(FUN = function(x) dnorm(x) * 
                                     sum(MASS:::psi.huber(1/sigma*X[,2]*(x-X%*%beta), k = 1.345)^(2))), 
                     lower = -10, upper = 10)$value
  
  ##integrate(f = Vectorize(FUN = function(x) dnorm(x) * 
  #                             sum(MASS:::psi.huber(1/sigma*X[,2]*(x-X%*%beta), k = 1.345)^(2))), 
  #             lower = -10, upper = 10)$value
  
  b2=robustBLME:::vpsi_huber(1/sigma*X[,2]*(y-X%*%beta),1.345,n*2)
  sum(b2)*1/sqrt(J)
}

r=1
R=2000
n=40
set.seed(2907)
DELTA=0.1

theta0=1
for (r in r:2000) {
  save.image("cont_80.RData")
  #cat(ree_margintheta0/r)
  y1=rnorm(n, mean =  115.5,4)
  y1=sort(y1,decreasing = F)
  y1[1:4]=115.5- abs(rcauchy(4,location=0,scale = 4)) #CONTAMINATION
  y2=rnorm(n,mean = 115.5+theta0, 4)
  
  Treatment=c(rep("New",n), rep("Srandard", n))
  dati=data.frame(x=Treatment, y=c( y1, y2))
  # dati=datireg
  gendata=function (R){
    thetai=matrix(ncol=3, nrow=R)
    thetai=t(apply(thetai, 1, function (c) proposal(1)))
    Y=apply(thetai,1,function (x) ygen(n=40,theta1=x[1],theta2=x[2],sigma =x[3] ))
    return(list(Y=Y, thetai=thetai[,2], lambdai=thetai[,1],sigmai=thetai[,3]))
  }
  ris=gendata(R =R)
  Y=ris$Y
  
  
  stat_median=sapply(1:length(Y),  function (x) -median(Y[[x]]$y[1:40] )+
                       median(Y[[x]]$y[41:80] ))
  s= median(y2)-median(y1)
  accok=  sapply(stat_median,  function (x) (s<=x))
  ris_median=list();   ris_median$theta=ris$thetai[accok]
  accok=  sapply(stat_median,  function (x)  abs(s-x)<=DELTA)
  ris_ABCmedian=list();   ris_ABCmedian$theta=ris$thetai[accok]
  
  #ree
  statREE=function(M=M,lambda=lambda,sigma=sigma, dati=dati){
    REE(newdata = dati,M = M,lambda=lambda, sigma = sigma)
  }
  M=rlm(y~x, data=dati )
  LM=lm(y~x, data=dati )
  stat_ree=sapply(1:length(Y),  
                  function (x)  statREE(M=M,ris$Y[[x]],
                                        lambda = ris$lambdai[x], sigma = ris$sigmai[x]))
  # s=   statREE(dati = dati,M = M,lambda = M$coefficients[1], sigma = M$s)
  s=sapply(1:1,  
           function (x)  statREE(M=M,dati= dati,
                                 lambda = M$coefficients[1], sigma = M$s ))
  accok=  sapply(1:length(stat_ree),  function (x) (s <=stat_ree[x]))
  ris_ree=list();   ris_ree$theta=ris$thetai[accok]
  accok=  sapply(1:length(stat_ree),  function (x)  abs(s -stat_ree[x])<=DELTA)
  ris_ABCree=list();   ris_ABCree$theta=ris$thetai[accok]
  
  #MEST
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
  
  
  est_LM=LM$coefficients[2]; sd_est_LM=summary(LM)$coefficients[2,2]
  est_M=M$coefficients[2]; sd_est_M=summary(M)$coefficients[2,2]
  #tt=t.test(dati$y[41:80],dati$y[1:40], data=dati, alternative="less",var.equal=T)
  #tt$statistic  #tt$p.value   #tt$conf.int
  # metRology::qt.scaled(0.95,df = 78,mean = est_LM ,sd=sd_est_LM)
  asy_LM=rt.scaled(10000,78,mean = est_LM, sd_est_LM)
  asy_M=rnorm(10000,mean = est_M, sd_est_M)
  # ris_cd_mean=to_cd(ris_mean$theta,ndens = 1500, from=-2,to=6,bw=0.35)
  ris_cd_median=to_cd(ris_median$theta,ndens = 1500, from=-2,to=6, bw=0.15)
  ris_cd_M=to_cd(ris_M$theta,ndens = 1500, from=-2,to=6, bw=0.15)
  ris_cd_ree=to_cd(ris_ree$theta,ndens = 1500, from=-2,to=6, bw=0.15)
  
  
  par(mfrow=c(2,2))
  
  #median
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
  
  median_margintheta0=median_margintheta0+as.numeric(quantile(ris_cd_median$theta,0.95)>theta0)
  ree_margintheta0=ree_margintheta0+as.numeric(quantile(ris_cd_ree$theta,0.95)>theta0)
  m_margintheta0=m_margintheta0+as.numeric(quantile(ris_cd_M$theta,0.95)>theta0)
  abc_median_margintheta0=median_margintheta0+as.numeric(quantile(ris_cd_median$theta,0.95)>theta0)
  abc_ree_margintheta0=ree_margintheta0+as.numeric(quantile(ris_cd_ree$theta,0.95)>theta0)
  abc_m_margintheta0=m_margintheta0+as.numeric(quantile(ris_cd_M$theta,0.95)>theta0)
  masy_margintheta0=masy_margintheta0+as.numeric(quantile(asy_M,0.95)>theta0)
  lmasy_margintheta0=lmasy_margintheta0+as.numeric(metRology::qt.scaled(0.95,df = 78,mean = est_LM ,sd=sd_est_LM)>theta0)
  
  #MEST
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
  #ree
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
  
  cat(r);   beta0=1
  print(c("cd median ",
          mean(apply(cbind(CD_hpd95_median,beta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
          mean(apply(cbind(CD_hpd90_median,beta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
          mean(apply(cbind(CD_q95_median,beta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
          mean(apply(cbind(CD_q90_median,beta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T)))
  
  print(c("abc median ",
          mean(apply(cbind(ABC_hpd95_median,beta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
          mean(apply(cbind(ABC_hpd90_median,beta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
          mean(apply(cbind(ABC_q95_median,beta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
          mean(apply(cbind(ABC_q90_median,beta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T)))
  
  
  print(c("cd m",
          mean(apply(cbind(CD_hpd95_m,beta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
          mean(apply(cbind(CD_hpd90_m,beta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
          mean(apply(cbind(CD_q95_m,beta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
          mean(apply(cbind(CD_q90_m,beta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T)))
  
  print(c("abc m",
          mean(apply(cbind(ABC_hpd95_m,beta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
          mean(apply(cbind(ABC_hpd90_m,beta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
          mean(apply(cbind(ABC_q95_m,beta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
          mean(apply(cbind(ABC_q90_m,beta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T)))
  
  
  print(c("cd ree",
          mean(apply(cbind(CD_hpd95_ree,beta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
          mean(apply(cbind(CD_hpd90_ree,beta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
          mean(apply(cbind(CD_q95_ree,beta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
          mean(apply(cbind(CD_q90_ree,beta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T)))
  
  print(c("abc ree",
          
          mean(apply(cbind(ABC_hpd95_ree,beta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
          mean(apply(cbind(ABC_hpd90_ree,beta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
          mean(apply(cbind(ABC_q95_ree,beta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
          mean(apply(cbind(ABC_q90_ree,beta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T)))
  cat(r); cat("\n")
  cat((ree_margintheta0)/r)
  cat((lmasy_margintheta0)/r)
}

save.image("cont_80_1.RData")












 