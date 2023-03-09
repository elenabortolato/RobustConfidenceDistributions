set.seed(1)
yyy=rgamma(9000000,2.92,2)
med=median(yyy)
meann=mean(yyy)
par(mar=c(0.,0,0,0))
end=quantile(yyy, 0.95)
start=quantile(yyy, 0.05)

library(TeachingDemos)
hpd=emp.hpd(yyy,0.5)
start1=hpd[1]
end1=hpd[2]
poly_range <- density(yyy)$x > start & density(yyy)$x < end# polygon x-range
poly_range1 <- density(yyy)$x > start1 & density(yyy)$x < end1# polygon x-range

dd=density(
  yyy)
#dd2=density(yyy2,kernel = "rectangular",bw=0.001,n=10000)
plot(dd,type="n",xlim=c(0,6.4),  xlab="",ylab=" ",
     main="")# non inferiority  NON CONTAMINATION


plot(dd,lwd=1,xlim=c(0,6.4), ylim=c(0,0.6) ,xlab="",ylab=" ", main="", zero.line = F,box=F, type="l")# non inferiority  NON CONTAMINATION
#
polygon(c(start, density(yyy)$x[poly_range], end),                # X-Coordinates of polygon range
        c(0, density(yyy)$y[poly_range], 0),                 # Y-Coordinates of polygon range
        lwd=0.1,border=1, 
        col  ="#B8DE29FF") 

polygon(c(start1, density(yyy)$x[poly_range1], end1),                # X-Coordinates of polygon range
        c(0, density(yyy)$y[poly_range1], 0), density = 12,                 # Y-Coordinates of polygon range
        border=1,col="#29AF7FFF",lwd=0.1 ) 

poly_range2 <- density(yyy)$x > quantile(yyy,0.975) # polygon x-range
polygon(c( quantile(yyy,0.975), density(yyy)$x[poly_range2], 7),                # X-Coordinates of polygon range
        c(0, density(yyy)$y[poly_range2], 0),                  # Y-Coordinates of polygon range
        col = "#FDE927FF", lwd=0.1) 
#start=quantile(yyy, 0.05)
#end=quantile(yyy, 0.95)

#seqq=seq(from=(start), to=end, length=20)
#x2=yyy2
#polygon(c(min(dd2$x[dd2$x>start]),dd2$x[dd2$x>start]),                # X-Coordinates of polygon
#       c(0, dd2$y[dd2$x>start]),border = 0, lty = 1,                                # Y-Coordinates of polygon
#       angle = 0,  col="#1b90e2")                                     # Color of polygon    
#polygon(border = c( 10.4,2), 1,col=3,density = 3.0)

text(x = c(4.20), y=0.074814,c("one sided p-value"),cex = 0.78)

#rm(list=ls())
# abline(v=quantile(yyy,0.975), lty=2)
abline(v=med, lty=3, col=9)
abline(v=meann, lty=3, col=9)
abline(v=0.950, lty=3, col=9)
text(med,0.5795, "point estimators", cex = 0.78)
text(med-0.039,-0.012, "med", cex = 0.78)
text((meann+0.07 ), -0.0132, "mean", cex=0.78)
text(0.97,-0.012, "mode", cex=0.78)
# abline(v=0.5035, lty=1, col=2)
text(x = c(1.77), y=0.187,paste("1-  ","% CI"),cex = 0.78)
text(x = c(1.73), y=0.187,expression(alpha),cex = 0.78)


text(0.9,0.1319843, expression(psi[1]), cex = 0.78)
text(0.98,0.1319843, "<", cex = 0.78)
text(1.05,0.1319843, expression(psi), cex = 0.78)
text(1.12,0.1319843, "<", cex = 0.78)
text(1.225,0.1319843, expression(psi[2]), cex = 0.78)
text(3.45,-0.012, expression(psi), cex = 0.78)

text(hpd[1]+0.08,0.0125, expression(psi[1]), cex = 0.78)
text(hpd[2]+0.08,0.0125, expression(psi[2]), cex = 0.78)
 