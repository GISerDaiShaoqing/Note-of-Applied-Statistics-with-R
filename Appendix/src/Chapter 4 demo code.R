# Name: Note of Applied Statistics with R(Chapter 4 demo code) Normal Distribution and visualization
# Purpose: Generation of Normal Distribution and visualization
# Author:      Dai shaoqing
#
# Created:     05/06/2017
# Copyright:   (c) Dai shaoqing <dsq1993qingge@163.com> 2017
#------------------------------------------------------------
#Set Fonts
windowsFonts(RT=windowsFont("Times New Roman"))

#Generation of Normal Distribution(Different mean and sd, same samples)
d<-rnorm(5000,mean=5,sd=1)
e<-rnorm(5000,mean=5,sd=2)
f<-rnorm(5000,mean=15,sd=1)

#Generation of range
m<-range(f)
n<-range(d)
xrange<-c(n[1],m[2])

#plot
jpeg(file = "ndplot.jpg",width = 5000,height = 4000,units = "px",res = 1000,family="RT")
plot(density(d),col="blue",xlim = xrange,family="RT")
lines(density(e),col="blue",family="RT")
lines(density(f),col="blue",family="RT")
text(10,0.1,labels = expression(sigma^2==4),family="RT")
text(8,0.3,labels = expression(sigma^2==1),family="RT")
text(17,0.3,labels = expression(sigma^2==1),family="RT")
dev.off()

#Plot function of Normal Distribution
x <- seq(-5,5,length.out=100)
y <- dnorm(x,0,1)

jpeg(file = "ndplot2.jpg",width = 500,height = 400,units = "px",res = 1000,family="RT")
plot(x,y,col="red",xlim=c(-5,5),ylim=c(0,1),type='l',
     xaxs="i", yaxs="i",ylab='density',xlab='',
     main="The Normal Density Distribution",family="RT")

lines(x,dnorm(x,0,0.5),col="green",family="RT")
lines(x,dnorm(x,0,2),col="blue",family="RT")
lines(x,dnorm(x,-2,1),col="orange",family="RT")

legend("topright",legend=paste("m=",c(0,0,0,-2)," sd=", c(1,0.5,2,1)), lwd=1, col=c("red", "green","blue","orange"))
dev.off()

#Generate Standard Normal Distribution and visualization
a<-rnorm(5000,mean=0,sd=1)
c<-1+log10(length(a))/log10(2)
c
jpeg(file = "ndplot3.jpg",width = 500,height = 400,units = "px",res = 1000,family="RT")
hist(a,breaks = 15,freq=F,col="lightblue")
lines(density(a),col="red")
dev.off()

#Calculate the Qd
Qd<-function(x) { 
  q=fivenum(x)
  Qd=q[4]-q[2]
  cat("The quartile deviation is ", Qd)
}

#Normal Test based on Qd and sd
Normaltestindex<-function(x) { 
  q=fivenum(x)
  Qd=q[4]-q[2]
  s=sd(x)
  Normaltestindex=Qd/s
  cat("The Qd/s :", Normaltestindex)
}

#Compare Three Disrtibution
#Population Distribution
X<-c(18,20,22,24)
jpeg(file="ndplot4.jpg",width = 500,height = 400,units = "px",res = 1000,family="RT")
hist(X,breaks = 5,freq = F)
dev.off()
MX<-mean(X)
MSD<-sqrt(sum((X-mean(X))^2)/(length(X)))

##Population sd or var
Populationsd<-function(x){
  n=length(x)
  m=mean(x)
  Psd=sqrt(sum((x-m)^2)/n)
  cat("The Standard deviation of Population : ",Psd)
}
Populationsd(X)

PopulationVar<-function(x){
  n=length(x)
  m=mean(x)
  Psd=sum((x-m)^2)/n
  cat("The Variance of Population : ",Psd)
}
PopulationVar(X)

#Construct Sample Distribution
sdm<-c(18,19,20,21,19,20,21,22,20,21,22,23,21,22,23,24)
sdmc<-1+log10(length(sdm))/log10(2)
sdmc
jpeg(file="ndplot5.jpg",width = 500,height = 400,units = "px",res = 1000,family="RT")
hist(sdm,breaks = 9,freq = F,col="lightblue")
dev.off()

#Plot Chi-Square Distribution
cx<-seq(1,20,1)
cy<-dchisq(cx,1)

jpeg(file = "ndplot6.jpg",width = 500,height = 400,units = "px",res = 1000,family="RT")
plot(cx,cy,col="red",xlim=c(1,20),ylim=c(0,0.25),type='l',
     xaxs="i", yaxs="i",ylab='density',xlab='',
     main="The Chi-Square Distribution",family="RT")

lines(cx,dchisq(cx,4),col="green",family="RT")
lines(cx,dchisq(cx,10),col="blue",family="RT")
lines(cx,dchisq(cx,20),col="orange",family="RT")

legend("topright",legend=paste("df=",c(1,4,10,20)), lwd=1, col=c("red", "green","blue","orange"))
dev.off()

#plot t Distribution
tx<-seq(-5,5,length.out = 100)
ty<-dnorm(tx,0,1)

jpeg(file = "ndplot7.jpg",width = 500,height = 400,units = "px",res = 1000,family="RT")
plot(tx,ty,col="red",xlim=c(-5,5),ylim=c(0,0.5),type='l',
     xaxs="i", yaxs="i",ylab='density',xlab='',
     main="The t Distribution",family="RT")

lines(tx,dt(tx,df=5),col="green",family="RT")
lines(tx,dt(tx,df=13),col="blue",family="RT")
text(0,0.45,labels = "red:Standard Normal Distribution")
text(2,0.2,labels = "green: df=5")
text(1.5,0.4,labels = "blue: df=13")

dev.off()

#plot F Distribution
fx<-seq(1,10,1)
fy<-df(fx,df1=1,df2=10)

jpeg(file="ndplot8.jpg",width = 500,height = 400,units = "px",res = 1000,family="RT")
plot(fx,fy,col="red",xlim = c(1,10),ylim=c(0,0.7),type="l",
     xaxs="i", yaxs="i",ylab='density',xlab='',
     main="The F Distribution",family="RT")

lines(fx,df(fx,df1=5,df2=10),col="green",family="RT")
lines(fx,df(fx,df1=10,df2=10),col="blue",family="RT")

legend("topright",legend=paste("df1=",c(1,4,10,20)," df2=", c(1,4,10,20)), lwd=1, col=c("red", "green","blue","orange"))
dev.off()