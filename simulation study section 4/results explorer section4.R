# results simulations M estimator CDs
load()
'filename="nocont_40.RData"

if(filename=="nocont_40.RData" ||
   filename=="nocont_80.RData" ||
   filename=="cont_40.RData" ||
   filename=="cont_80.RData"
   ){
  load(filename)
  theta0=2.6
}

if(filename=="nocont_40_1.RData"||
   filename=="nocont_80_1.RData" ||
   filename=="cont_40_1.RData" ||
   filename=="cont_80_1.RData"){
  load(filename)
  theta0=1
}
'
###############################################à
#bias
b=(cbind(c("Wald/Mean", "wald/M-test", "ABC/Median", "ABC/M-EE", "ABC/M-est", 
           "CDensity/Median", "CDensity/M-EE", "CDensity/M-est"),
                       abs(rbind(
                       summary(CDmedian_asy_lm, digits=2),
                       summary(CDmedian_asy_m, digits=2),
                       summary(ABCmedian_median, digits=2),
                       summary(ABCmedian_ree, digits=2),
                       summary(ABCmedian_m, digits=2),
                       summary(CDmedian_median, digits=2),
                       summary(CDmedian_ree, digits=2),
                       summary(CDmedian_m, digits=2))[,4]-theta0)))

#prob underestimation
pu=(cbind(c("Wald/Mean", "wald/M-test", "ABC/Median", "ABC/M-EE", "ABC/M-est", 
            "CDensity/Median", "CDensity/M-EE", "CDensity/M-est"),
          (rbind(
                       mean((CDmedian_asy_lm)>theta0),
                       mean((CDmedian_asy_m)>theta0),
                       mean((ABCmedian_median)>theta0),
                       mean((ABCmedian_ree)>theta0),
                       mean((ABCmedian_m) >theta0),
                       mean((CDmedian_median)>theta0),
                       mean((CDmedian_ree)>theta0),
                       mean((CDmedian_m) >theta0)))))



error=round(matrix(c(1-(lmasy_theta0)/r,
1-(masy_theta0)/r,
1-(abc_median_theta0)/r,
1-(abc_m_theta0)/r,
1-(abc_ree_theta0)/r,
1-(median_theta0)/r,
1-(m_theta0)/r,
1-(ree_theta0)/r
)),2)

xtable::xtable(cbind(b,pu, error)[,c(1,2,4,5)],digits=2)
#coverage of intervals
INT=rbind(
rbind(c(mean(apply(cbind(CD_hpd95_asy_lm,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
mean(apply(cbind(CD_hpd90_asy_lm,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T)),
c(mean(apply(cbind(CD_q95_asy_lm,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
mean(apply(cbind(CD_q90_asy_lm,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T))),
rbind(c(mean(apply(cbind(CD_hpd95_asy_m,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
mean(apply(cbind(CD_hpd90_asy_m,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T)),
c(mean(apply(cbind(CD_q95_asy_m,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
mean(apply(cbind(CD_q90_asy_m,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T))),
rbind(c(mean(apply(cbind(ABC_hpd95_median,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
mean(apply(cbind(ABC_hpd90_median,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T)),
c(mean(apply(cbind(ABC_q95_median,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
mean(apply(cbind(ABC_q90_median,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T))),
rbind(c(mean(apply(cbind(CD_hpd95_ree,theta0),1, function (x)   x[3]>x[1] && x[3]<x[2]),na.rm = T),
        mean(apply(cbind(CD_hpd90_ree,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T)),
      c(mean(apply(cbind(CD_q95_asy_lm,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
        mean(apply(cbind(CD_q90_asy_lm,theta0),1, function (x)   x[3]>x[1] && x[3]<x[2]),na.rm = T))),
rbind(c(mean(apply(cbind(ABC_hpd95_m,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
mean(apply(cbind(ABC_hpd90_m,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T)),
c(mean(apply(cbind(ABC_q95_m,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
mean(apply(cbind(ABC_q90_m,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T))),
rbind(c(mean(apply(cbind(CD_hpd95_median,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
        mean(apply(cbind(CD_hpd90_median,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T)),
      c(mean(apply(cbind(CD_q95_median,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
      mean(apply(cbind(CD_q90_median,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T))),
rbind(c(
mean(apply(cbind(CD_hpd95_ree,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
mean(apply(cbind(CD_hpd90_ree,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T)),
c(mean(apply(cbind(CD_q95_ree,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
mean(apply(cbind(CD_q90_ree,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T))),
rbind(c(
  mean(apply(cbind(CD_hpd95_m,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
  mean(apply(cbind(CD_hpd90_m,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T)),
  c(mean(apply(cbind(CD_q95_m,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T),
    mean(apply(cbind(CD_q90_m,theta0),1, function (x)  x[3]>x[1] && x[3]<x[2]),na.rm = T))))
 
INT=INT[c(2,4,6,8,10,12,14,16),]
rownames(INT)=(c("Wald/Mean", "wald/M-test", "ABC/Median", "ABC/M-EE", "ABC/M-est", 
                 "CDensity/Median", "CDensity/M-EE", "CDensity/M-est"))


xtable::xtable(INT)
###############################################################################