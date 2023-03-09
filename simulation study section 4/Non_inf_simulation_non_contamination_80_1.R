#noninf NON CONTAMINATED SIMUALTION STUDY


###################################################################
rm(list=ls())

#SIMULATION STUDY non contaminated

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


to_cd=function(values=val[[i]], ndens=500, from=NULL, to=NULL, bw =0.3){
  if (length(values)<=2) {
    values=c(-999,-999)
    return(list(theta=values))
  }
  values=values[!is.na(values)]
  F1=density(values, n=ndens, from=from,bw=bw,  to=to,kernel="rectangular")
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
ree_theta0=0
abc_ree_theta0=0
m_theta0=0
abc_m_theta0=0
median_theta0=0
abc_median_theta0=0
lmasy_theta0=0
masy_theta0=0



REE=function (M=M,newdata=dati, lambda=lambda, sigma=sigma){
  y=newdata$y
  X=model.matrix(M)
  beta=M$coefficients
  beta[1]=lambda
  J <- integrate(f = Vectorize(FUN = function(x) dnorm(x) * 
        sum(MASS:::psi.huber(1/sigma*X[,2]*(x-X%*%beta), k = 1.345)^(2))), 
        lower = -10, upper = 10)$value
  
 # J <-1
  b2=robustBLME:::vpsi_huber(1/sigma*X[,2]*(y-X%*%beta),1.345,n*2)
  sum(b2)*1/sqrt(J)
}

r=1
R=2000
n=40 #for each group
set.seed(2907)
DELTA=0.1

for (r in r:2000) {
  y1=rnorm(n, mean =  115.5,4)
  y2=rnorm(n,mean = 115.5+theta0, 4)
  
  Treatment=c(rep("New",n), rep("Standard", n))
  dati=data.frame(x=Treatment, y=c( y1, y2))
  
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
  ris_median=list()
  ris_median$theta=ris$thetai[accok]
  accok=  sapply(stat_median,  function (x)  abs(s-x)<=DELTA)
  ris_ABCmedian=list()
  ris_ABCmedian$theta=ris$thetai[accok]
  
  
  
  #ree
  statREE=function(M=M,lambda=lambda,sigma=sigma, dati=dati){
    REE(newdata = dati,M = M,lambda=lambda, sigma = sigma)
  }
  M=rlm(y~x, data=dati )
  stat_ree=sapply(1:length(Y),  
                  function (x)  statREE(M=M,ris$Y[[x]],
                                        lambda = ris$lambdai[x], sigma = ris$sigmai[x]))
  # s=   statREE(dati = dati,M = M,lambda = M$coefficients[1], sigma = M$s)
  s=sapply(1:1,  
           function (x)  statREE(M=M,dati= dati,
                                 lambda = M$coefficients[1], sigma = M$s ))
  accok=  sapply(1:length(stat_ree),  function (x) (s <=stat_ree[x]))
  ris_ree=list()
  ris_ree$theta=ris$thetai[accok]
  accok=  sapply(1:length(stat_ree),  function (x)  abs(s -stat_ree[x])<=DELTA)
  ris_ABCree=list()
  ris_ABCree$theta=ris$thetai[accok]
  
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
  ris_M=list()
  ris_M$theta=ris$thetai[accok]
  accok=  sapply(stat_M,  function (x)  abs(s-x)<=DELTA)
  ris_ABCM=list()
  ris_ABCM$theta=ris$thetai[accok]
  
  LM=lm(y~x, data=dati); 
  est_LM=LM$coefficients[2]; sd_est_LM=summary(LM)$coefficients[2,2]
  est_M=M$coefficients[2]; sd_est_M=summary(M)$coefficients[2,2]
  
  asy_LM=rnorm(10000,mean = est_LM, sd_est_LM)
  asy_M=rnorm(10000,mean = est_M, sd_est_M)
   ris_cd_median=to_cd(ris_median$theta,ndens = 1500, from=-2,to=6, bw=0.35)
  ris_cd_M=to_cd(ris_M$theta,ndens = 1500, from=-2,to=6, bw=0.35)
  ris_cd_ree=to_cd(ris_ree$theta,ndens = 1500, from=-2,to=6, bw=0.35)
  
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
  
  #median_margin2=median_margin2+as.numeric(quantile(ris_cd_median$theta,0.975)>margin2)
  median_theta0=median_theta0+as.numeric(quantile(ris_cd_median$theta,0.95)>theta0)
   
  #ree_margin2=ree_margin2+as.numeric(quantile(ris_cd_ree$theta,0.975)>margin2)
  ree_theta0=ree_theta0+as.numeric(quantile(ris_cd_ree$theta,0.95)>theta0)
  #ree_margin3=ree_margin3+as.numeric(quantile(ris_cd_ree$theta,0.975)>margin3)
  
  
  #m_margin2=m_margin2+as.numeric(quantile(ris_cd_M$theta,0.975)>margin2)
  m_theta0=m_theta0+as.numeric(quantile(ris_cd_M$theta,0.95)>theta0)
  masy_theta0=masy_theta0+as.numeric(quantile(asy_M,0.95)>theta0)
  abc_median_theta0=abc_median_theta0+as.numeric(quantile(ris_ABCmedian,0.95)>theta0)
  abc_m_theta0=abc_m_theta0+as.numeric(quantile(ris_ABCM,0.95)>theta0)
  abc_ree_theta0=abc_ree_theta0+as.numeric(quantile(ris_ABCree,0.95)>theta0)
  lmasy_theta0=lmasy_theta0+as.numeric(quantile(asy_LM,0.95)>theta0)
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
  #ABCmedian_mean[r]=median(ris_ABCmean$theta)
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
  
  cat((ree_theta0)/r)
}
 


save.image("nocont_80_1.RData")

############ 
############


# FIGURE 1

#install.packages("tidyverse")
library(tidyverse)
#install.packages("hrbrthemes")
library(hrbrthemes)
#install.packages("kableExtra")
library(kableExtra)
options(knitr.table.format = "html")
library(viridis)
library(ggplot2)
# Load dataset from github

REP=10^5*5
data <- data.frame(type=c(rep("Higher is better",REP), rep("Lower is better",REP),
                          rep("Higher is better",REP), rep("Lower is better",REP),
                          rep("Higher is better",REP), rep("Lower is better",REP),
                          rep("Higher is better",REP), rep("Lower is better",REP)),
                   result=c(rep("Superiority not established",REP),
                            rep("Superiority not established",REP),
                            rep("Superiority established",REP)
                            ,rep("Superiority established",REP),
                            rep("Non inferiority established",REP),
                            rep("Non inferiority established",REP),
                            rep("Non inferiority not established",REP)
                            ,rep("Non inferiority not established",REP)),
                   delta=c(rnorm(REP,0.5,1),rnorm(REP,-0.5,1),
                           rnorm(REP,2.75,1),rnorm(REP,-2.75,1),
                           rnorm(REP,-1.8,1),rnorm(REP,1.84,1),
                           rnorm(REP,-5.5,1),rnorm(REP,5.5,1)),
                   Test=c(rep("Traditional",REP*4),
                          rep("Non inferiority",REP*4)),
                   pos=c(rep(-5,REP  ),rep(5,REP )),
                   
                   w=c(rep("A",REP),rep("B",REP),
                       rep("A",REP),rep("B",REP),
                       rep("C",REP),rep("D",REP),
                       rep("C",REP),rep("D",REP)))
data$type=as.factor(data$type)
levels(data$type)=c("Higher is better", "Lower is better")


dataA=data[data$Test=="Traditional",]
#colnames(dataA)[2]="result"
dataB=data[data$Test=="Non inferiority",]
dataC=data[data$w=="C",]
dataD=data[data$w=="D",]
#install.packages("remotes")
#remotes::install_github("coolbutuseless/ggpattern")
library(ggpattern)
pA=ggplot( data = dataA , aes(x=delta )) +
  scale_fill_viridis(discrete=TRUE) +
  # geom_col_pattern(aes(x=delta,pattern_color=result))+
  geom_density_pattern(aes(pattern_density=(result)) ) +
  scale_color_viridis(discrete=TRUE) +
  theme_classic() +
  geom_vline( aes(xintercept=0))+
  theme(legend.title = element_blank(),
        #legend.position="none",
        panel.spacing = unit(0.1, "lines"),
        strip.text.x = element_text(size = 8)
  ) +
  xlab("Treatment difference (new-standard)") +
  facet_grid(cols=vars(type))+
  #  facet_grid(~result )+ult
  ylab("Traditional")
pA
pB=ggplot( data = dataB,aes(  x=delta )) +
  
  scale_fill_viridis(discrete=TRUE) +
  geom_density_pattern(aes(pattern_density=(result)) ) +
  scale_color_viridis(discrete=TRUE) +
  theme_classic() +
  geom_vline( aes(xintercept=pos))+
  theme(legend.title = element_blank(),
        #theme(legend.position = "top",
        #legend.position="none",
        panel.spacing = unit(0.1, "lines"),
        strip.text.x = element_text(size = 8)
  ) +
  
  xlab("Treatment difference (new-standard)") +
  #  facet_grid(~result )+ult
  facet_wrap(~ type ,shrink = T,scales = "free_x" ,drop = T)+
  ylab("Non inferiority")
pB


#install.packages("devtools")
#library(devtools)
#install_github("easyGgplot2", "kassambara")
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
#multiplot(pA, pB)
#install.packages("ggpubr")
ggpubr::ggarrange(pA,pB, nrow = 2)
 