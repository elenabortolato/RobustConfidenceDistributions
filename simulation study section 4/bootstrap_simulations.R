# simulation bootstrap
# model :

#y= alpha + theta0 group_standard + error
 
# first scenario
theta0=2.6 

# second scenario
# theta0=1 

######################################
# 40 per group (80), no contamination
n=40
r=1
library(MASS)
library(boot)

# computes the parameter of inerest
test.fun <- function(data) {
  M=rlm(y~x, data=data)
  M$coefficients[2]
}

Treatment=c(rep("New",n), rep("Standard", n))

# data generator
test.rg <- function( data, mle=est0) {
  y1=rnorm(n, mean =  mle[1],mle[3])
  y2=rnorm(n,mean = mle[1]+mle[2],mle[3])
  dati=data.frame(x=Treatment, y=c( y1, y2))
  dati
}

basic=norm=perc=matrix(nrow=2000,ncol=2)
basic95=norm95=perc95=matrix(nrow=2000,ncol=2)
M_basic=M_norm=M_perc=matrix(nrow=2000,ncol=2)

set.seed(90810)
for (r in r:2000) {
  y1=rnorm(n, mean =  115.5,4)
  y2=rnorm(n,mean = 115.5+theta0, 4)
  
  Treatment=c(rep("New",n), rep("Standard", n)) # OPTIMIZED 
  dati=data.frame(x=Treatment, y=c( y1, y2))
  M0=rlm(y~x, data=dati)
  est0=c(M0$coefficients[1:2],M0$s)
 
  b=boot(data=dati, sim="parametric",R = 2000,ncpus=4,
         statistic = test.fun,ran.gen = test.rg,mle = est0,parallel = "multicore")
  ci=boot.ci(b, type = c("basic","norm", "perc"),conf = c(0.9,0.95,0))

  basic[r,]=c(ci$basic[1,4], ci$basic[1,5])
  norm[r,]=c(ci$normal[1,2], ci$normal[1,3])
  perc[r,]=c(ci$perc[1,4], ci$perc[1,5])
  
  basic95[r,]=c(ci$basic[2,4], ci$basic[2,5])
  norm95[r,]=c(ci$normal[2,2], ci$normal[2,3])
  perc95[r,]=c(ci$perc[2,4], ci$perc[2,5])
  
  M_basic[r,]=c(ci$basic[3,4] )
  M_norm[r,]=c(ci$normal[3,2])
  M_perc[r,]=c(ci$perc[3,4])
  cat("\n")
  cat("\n")
  
  cat(mean( apply(basic, 1,  function (x) x[1]<=theta0 && x[2]>=theta0), na.rm=T))
  cat(mean( apply(norm, 1,  function (x) x[1]<=theta0 && x[2]>=theta0), na.rm=T))
  cat(mean( apply(perc, 1,  function (x) x[1]<=theta0 && x[2]>=theta0), na.rm=T))
  cat(r)
}


intboot80=xtable::xtable(100*cbind(rbind(
(mean( apply(basic95, 1,  function (x) x[1]<=theta0 && x[2]>=theta0), na.rm=T)),
(mean( apply(norm95, 1,  function (x) x[1]<=theta0 && x[2]>=theta0), na.rm=T)),
(mean( apply(perc95, 1,  function (x) x[1]<=theta0 && x[2]>=theta0), na.rm=T))),
rbind(
(mean( apply(basic, 1,  function (x) x[1]<=theta0 && x[2]>=theta0), na.rm=T)),
(mean( apply(norm, 1,  function (x) x[1]<=theta0 && x[2]>=theta0), na.rm=T)),
(mean( apply(perc, 1,  function (x) x[1]<=theta0 && x[2]>=theta0), na.rm=T)))))

errboot80=xtable::xtable(
cbind(
rbind(  

c(mean( apply(M_basic, 1,  function (x) x[1]<theta0), na.rm=T)),
c(mean( apply(M_norm, 1,  function (x) x[1]<theta0), na.rm=T)),
c(mean( apply(M_perc, 1,  function (x) x[1]<theta0), na.rm=T))),


rbind(c(mean( apply(M_basic, 1,  function (x) x[1]-theta0), na.rm=T)),
c(mean( apply(M_norm, 1,  function (x) x[1]-theta0), na.rm=T)),
c(mean( apply(M_perc, 1,  function (x) x[1]-theta0), na.rm=T))),


rbind(c(mean( apply(basic, 1,  function (x)  x[2]<theta0), na.rm=T)),
c(mean( apply(norm, 1,  function (x)  x[2]<theta0), na.rm=T)),
c(mean( apply(perc, 1,  function (x)  x[2]<theta0), na.rm=T)))))
 




######################################
# 40 per group (80),  contamination
n=40
r=1

basic=norm=perc=matrix(nrow=2000,ncol=2)
basic95=norm95=perc95=matrix(nrow=2000,ncol=2)
set.seed(9876)
for (r in r:2000) {
  y1=rnorm(n, mean =  115.5,4)
  y1=sort(y1,decreasing = F)
  y1[1:4]=115.5- abs(rcauchy(4,location=0,scale = 4)) #CONTAMINATION
  y2=rnorm(n,mean = 115.5+theta0, 4)
  
  Treatment=c(rep("New",n), rep("Standard", n))
  dati=data.frame(x=Treatment, y=c( y1, y2))
  
  M0=rlm(y~x, data=dati)
  est0=c(M0$coefficients[1:2],M0$s)
  est0
  test.fun <- function(data) {
    M=rlm(y~x, data=data)
    M$coefficients[2]
  }
  
  test.rg <- function( data, mle=est0) {
    y1=rnorm(n, mean =  mle[1],mle[3])
    y2=rnorm(n,mean = mle[1]+mle[2],mle[3])
    Treatment=c(rep("New",n), rep("Standard", n))
    dati=data.frame(x=Treatment, y=c( y1, y2))
    dati
  }
  b=boot(data=dati, sim="parametric",R = 2000,statistic = test.fun,ran.gen = test.rg,mle = est0,parallel = "multicore")
   
  ci=boot.ci(b, type = c("basic","norm", "perc"),conf = c(0.9,0.95,0.0))
  #  ciM=boot.ci(b, type = c("basic","norm", "perc"),conf = c(0.5))
  
  
  basic[r,]=c(ci$basic[1,4], ci$basic[1,5])
  norm[r,]=c(ci$normal[1,2], ci$normal[1,3])
  perc[r,]=c(ci$perc[1,4], ci$perc[1,5])
  
  basic95[r,]=c(ci$basic[2,4], ci$basic[2,5])
  norm95[r,]=c(ci$normal[2,2], ci$normal[2,3])
  perc95[r,]=c(ci$perc[2,4], ci$perc[2,5])
  M_basic[r,]=c(ci$basic[3,4] )
  M_norm[r,]=c(ci$normal[3,2])
  M_perc[r,]=c(ci$perc[3,4])
  cat("\n")
  cat("\n")
  
  cat(mean( apply(basic, 1,  function (x) x[1]<theta0 && x[2]>theta0), na.rm=T))
  cat(mean( apply(norm, 1,  function (x) x[1]<theta0 && x[2]>theta0), na.rm=T))
  cat(mean( apply(perc, 1,  function (x) x[1]<theta0 && x[2]>theta0), na.rm=T))
  cat(r)
}



intboot80cont=xtable::xtable(100*cbind(rbind(
  (mean( apply(basic95, 1,  function (x) x[1]<=theta0 && x[2]>=theta0), na.rm=T)),
  (mean( apply(norm95, 1,  function (x) x[1]<=theta0 && x[2]>=theta0), na.rm=T)),
  (mean( apply(perc95, 1,  function (x) x[1]<=theta0 && x[2]>=theta0), na.rm=T))),
  rbind(
    (mean( apply(basic, 1,  function (x) x[1]<=theta0 && x[2]>=theta0), na.rm=T)),
    (mean( apply(norm, 1,  function (x) x[1]<=theta0 && x[2]>=theta0), na.rm=T)),
    (mean( apply(perc, 1,  function (x) x[1]<=theta0 && x[2]>=theta0), na.rm=T)))))

errboot80cont=xtable::xtable(
  cbind(
    rbind(  
      
      c(mean( apply(M_basic, 1,  function (x) x[1]<theta0), na.rm=T)),
      c(mean( apply(M_norm, 1,  function (x) x[1]<theta0), na.rm=T)),
      c(mean( apply(M_perc, 1,  function (x) x[1]<theta0), na.rm=T))),
    
    
    rbind(c(mean( apply(M_basic, 1,  function (x) x[1]-theta0), na.rm=T)),
          c(mean( apply(M_norm, 1,  function (x) x[1]-theta0), na.rm=T)),
          c(mean( apply(M_perc, 1,  function (x) x[1]-theta0), na.rm=T))),
    
    
    rbind(c(mean( apply(basic, 1,  function (x)  x[2]<theta0), na.rm=T)),
          c(mean( apply(norm, 1,  function (x)  x[2]<theta0), na.rm=T)),
          c(mean( apply(perc, 1,  function (x)  x[2]<theta0), na.rm=T)))))
 
######################################
# 20 per group (40), no contamination
n=20
r=1
test.fun <- function(data) {
  M=rlm(y~x, data=data)
  M$coefficients[2]
}
Treatment=c(rep("New",n), rep("Standard", n))
test.rg <- function( data, mle=est0) {
  y1=rnorm(n, mean =  mle[1],mle[3])
  y2=rnorm(n,mean = mle[1]+mle[2],mle[3])
  dati=data.frame(x=Treatment, y=c( y1, y2))
  dati
}

basic=norm=perc=matrix(nrow=2000,ncol=2)
basic95=norm95=perc95=matrix(nrow=2000,ncol=2)
M_basic=M_norm=M_perc=matrix(nrow=2000,ncol=2)


r=1
set.seed(90810)
for (r in r:2000) {
  y1=rnorm(n, mean =  115.5,4)
  y2=rnorm(n,mean = 115.5+theta0, 4)
  
  Treatment=c(rep("New",n), rep("Standard", n)) # OPTIMIZED 
  dati=data.frame(x=Treatment, y=c( y1, y2))
  #dati$y=c( y1, y2)
  M0=rlm(y~x, data=dati)
  est0=c(M0$coefficients[1:2],M0$s)
  
  b=boot(data=dati, sim="parametric",R = 2000,ncpus=4,
         statistic = test.fun,ran.gen = test.rg,mle = est0,parallel = "multicore")
  ci=boot.ci(b, type = c("basic","norm", "perc"),conf = c(0.9,0.95,0.0))
  #  ciM=boot.ci(b, type = c("basic","norm", "perc"),conf = c(0.5))
  
  
  basic[r,]=c(ci$basic[1,4], ci$basic[1,5])
  norm[r,]=c(ci$normal[1,2], ci$normal[1,3])
  perc[r,]=c(ci$perc[1,4], ci$perc[1,5])
  
  basic95[r,]=c(ci$basic[2,4], ci$basic[2,5])
  norm95[r,]=c(ci$normal[2,2], ci$normal[2,3])
  perc95[r,]=c(ci$perc[2,4], ci$perc[2,5])
  
  M_basic[r,]=c(ci$basic[3,4] )
  M_norm[r,]=c(ci$normal[3,2])
  M_perc[r,]=c(ci$perc[3,4])
  cat("\n")
  cat("\n")
  
  cat(mean( apply(basic, 1,  function (x) x[1]<=theta0 && x[2]>=theta0), na.rm=T))
  cat(mean( apply(norm, 1,  function (x) x[1]<=theta0 && x[2]>=theta0), na.rm=T))
  cat(mean( apply(perc, 1,  function (x) x[1]<=theta0 && x[2]>=theta0), na.rm=T))
  cat(r)
}


intboot40=xtable::xtable(100*cbind(rbind(
  (mean( apply(basic95, 1,  function (x) x[1]<=theta0 && x[2]>=theta0), na.rm=T)),
  (mean( apply(norm95, 1,  function (x) x[1]<=theta0 && x[2]>=theta0), na.rm=T)),
  (mean( apply(perc95, 1,  function (x) x[1]<=theta0 && x[2]>=theta0), na.rm=T))),
  rbind(
    (mean( apply(basic, 1,  function (x) x[1]<=theta0 && x[2]>=theta0), na.rm=T)),
    (mean( apply(norm, 1,  function (x) x[1]<=theta0 && x[2]>=theta0), na.rm=T)),
    (mean( apply(perc, 1,  function (x) x[1]<=theta0 && x[2]>=theta0), na.rm=T)))))

errboot40=xtable::xtable(
  cbind(
    rbind(  
      
      c(mean( apply(M_basic, 1,  function (x) x[1]<theta0), na.rm=T)),
      c(mean( apply(M_norm, 1,  function (x) x[1]<theta0), na.rm=T)),
      c(mean( apply(M_perc, 1,  function (x) x[1]<theta0), na.rm=T))),
    
    
    rbind(c(mean( apply(M_basic, 1,  function (x) x[1]-theta0), na.rm=T)),
          c(mean( apply(M_norm, 1,  function (x) x[1]-theta0), na.rm=T)),
          c(mean( apply(M_perc, 1,  function (x) x[1]-theta0), na.rm=T))),
    
    
    rbind(c(mean( apply(basic, 1,  function (x)  x[2]<theta0), na.rm=T)),
          c(mean( apply(norm, 1,  function (x)  x[2]<theta0), na.rm=T)),
          c(mean( apply(perc, 1,  function (x)  x[2]<theta0), na.rm=T)))))

######################################
# 20 per group (40),  contamination
n=20
r=1 

basic=norm=perc=matrix(nrow=2000,ncol=2)
basic95=norm95=perc95=matrix(nrow=2000,ncol=2)
set.seed(9876)
for (r in r:2000) {
   
  y1=rnorm(n, mean =  115.5,4)
  y1=sort(y1,decreasing = F)
  y1[1:4]=115.5- abs(rcauchy(4,location=0,scale = 4)) #CONTAMINATION
  y2=rnorm(n,mean = 115.5+theta0, 4)
  
  Treatment=c(rep("New",n), rep("Standard", n))
  dati=data.frame(x=Treatment, y=c( y1, y2))
  
  M0=rlm(y~x, data=dati)
  est0=c(M0$coefficients[1:2],M0$s)
  est0
  test.fun <- function(data) {
    M=rlm(y~x, data=data)
    M$coefficients[2]
  }
  
  test.rg <- function( data, mle=est0) {
    y1=rnorm(n, mean =  mle[1],mle[3])
    y2=rnorm(n,mean = mle[1]+mle[2],mle[3])
    Treatment=c(rep("New",n), rep("Standard", n))
    dati=data.frame(x=Treatment, y=c( y1, y2))
    dati
  }
  b=boot(data=dati, sim="parametric",R = 2000,statistic = test.fun,ran.gen = test.rg,mle = est0,parallel = "multicore")
  ci=boot.ci(b, type = c("basic","norm", "perc"),conf = c(0.9,0.95,0.))
  #  ciM=boot.ci(b, type = c("basic","norm", "perc"),conf = c(0.5))
  
  
  basic[r,]=c(ci$basic[1,4], ci$basic[1,5])
  norm[r,]=c(ci$normal[1,2], ci$normal[1,3])
  perc[r,]=c(ci$perc[1,4], ci$perc[1,5])
  
  basic95[r,]=c(ci$basic[2,4], ci$basic[2,5])
  norm95[r,]=c(ci$normal[2,2], ci$normal[2,3])
  perc95[r,]=c(ci$perc[2,4], ci$perc[2,5])
  
  M_basic[r,]=c(ci$basic[3,4] )
  M_norm[r,]=c(ci$normal[3,2])
  M_perc[r,]=c(ci$perc[3,4])
  cat("\n")
  cat("\n")
  
  cat(mean( apply(basic, 1,  function (x) x[1]<=theta0 && x[2]>=theta0), na.rm=T))
  cat(mean( apply(norm, 1,  function (x) x[1]<=theta0 && x[2]>=theta0), na.rm=T))
  cat(mean( apply(perc, 1,  function (x) x[1]<=theta0 && x[2]>=theta0), na.rm=T))
  cat(r)
}

 
intboot40cont=xtable::xtable(100*cbind(rbind(
  (mean( apply(basic95, 1,  function (x) x[1]<=theta0 && x[2]>=theta0), na.rm=T)),
  (mean( apply(norm95, 1,  function (x) x[1]<=theta0 && x[2]>=theta0), na.rm=T)),
  (mean( apply(perc95, 1,  function (x) x[1]<=theta0 && x[2]>=theta0), na.rm=T))),
  rbind(
    (mean( apply(basic, 1,  function (x) x[1]<=theta0 && x[2]>=theta0), na.rm=T)),
    (mean( apply(norm, 1,  function (x) x[1]<=theta0 && x[2]>=theta0), na.rm=T)),
    (mean( apply(perc, 1,  function (x) x[1]<=theta0 && x[2]>=theta0), na.rm=T)))))

errboot40cont=xtable::xtable(
  cbind(
    rbind(  
      
      c(mean( apply(M_basic, 1,  function (x) x[1]<theta0), na.rm=T)),
      c(mean( apply(M_norm, 1,  function (x) x[1]<theta0), na.rm=T)),
      c(mean( apply(M_perc, 1,  function (x) x[1]<theta0), na.rm=T))),
    
    
    rbind(c(mean( apply(M_basic, 1,  function (x) x[1]-theta0), na.rm=T)),
          c(mean( apply(M_norm, 1,  function (x) x[1]-theta0), na.rm=T)),
          c(mean( apply(M_perc, 1,  function (x) x[1]-theta0), na.rm=T))),
    
    
    rbind(c(mean( apply(basic, 1,  function (x)  x[2]<theta0), na.rm=T)),
          c(mean( apply(norm, 1,  function (x)  x[2]<theta0), na.rm=T)),
          c(mean( apply(perc, 1,  function (x)  x[2]<theta0), na.rm=T)))))


 