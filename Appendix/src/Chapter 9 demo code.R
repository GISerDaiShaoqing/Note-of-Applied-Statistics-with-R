# Name: Note of Applied Statistics with R(Chapter 9 demo code) Linear Regression
# Purpose: Linear Regression
# Author:      Dai shaoqing
#
# Created:     06/11/2017
# Copyright:   (c) Dai shaoqing <dsq1993qingge@163.com> 2017
#------------------------------------------------------------
#Certain vs Uncertain
a<-seq(1,8,1)
b<-seq(4,32,4)
c<-seq(32,4,-4)
d<-c(3,7,10,18,25,20,16,12)
e<-c(3,7,11,17,21,23,29,33)
f<-c(33,29,23,21,17,11,7,3)
g<-runif(8,min=1,max=80)

modelb<-lm(b~a)
modelc<-lm(c~a)
modeld<-lm(d~I(a^2)+a)
modele<-lm(e~a)
modelf<-lm(f~a)

#different correlationships
jpeg("plot.jpeg",width=12,height=8,units="in",res=300)
layout(matrix(c(1,2,3,4,5,6),nrow=2,byrow=T))
plot(a,b,pch=16,col="red",main="完全正线性相关")
abline(modelb,col="blue")
plot(a,c,pch=16,col="red",main="完全负线性相关")
abline(modelc,col="blue")
plot(a,d,pch=16,col="red",main="非线性相关")
lines(a,fitted(modeld),col="blue")
plot(a,e,pch=16,col="red",main="正线性相关")
abline(modele,col="blue")
plot(a,f,pch=16,col="red",main="负线性相关")
abline(modelf,col="blue")
plot(a,g,pch=16,col="red",main="不相关")
dev.off()

#confidencel interval
#new data, predict and confidence inteval estimate 
h<-data.frame(h=seq(1,8,1))
con<-predict.lm(modele,h,interval="confidence",level=0.95)
pre<-predict.lm(modele,h,interval="prediction",level=0.95)

#plot
summary(modele)
jpeg("plot2.jpeg",width=8,height=6,units="in",res=300)
#layout(matrix(c(1),nrow=1,byrow=T))
plot(a,e,type="n",xlab="", ylab="", xlim=c(1,8),  ylim=c(0, 40))
polygon(c(h[,1], rev(h[,1])), c(con[,3], rev(con[,2])),border="red",lwd=1,lty = c("dashed", "solid"))
polygon(c(h[,1], rev(h[,1])), c(pre[,3], rev(pre[,2])),border="blue",lwd=1,lty = c("dashed", "solid"))
points(a,e,pch=16,col="blue")
abline(modele,col="black")
text(3,30,labels=expression(paste("y=4.2857x-1.2857,","Adj. ", italic(R)^2, "=0.9941, ", italic(P), "<0.01")), cex=1, adj=0)
dev.off()

#residual plot
jpeg("plot3.jpeg",width=8,height=6,units="in",res=300)
layout(matrix(c(1,2,3,4),nrow=2,byrow=T))
plot(modele)
dev.off()

#plot 3d
x<-seq(1,8,1)
y<-seq(1,24,3)
f<-function(x,y){
  r<-3*x+y
}
z<-outer(x,y,f)
jpeg("plot4.jpeg",width=8,height=6,units="in",res=300)
#layout(matrix(c(1),nrow=1,byrow=T))
persp(x,y,z,theta=210,phi=35,col="red3")
dev.off()

#plot quadratic regression model
model.qr<-function(x,beta1=2,beta2=3){
  r<-4.5+beta1*x+beta2*x^2
  e<-runif(length(x),-1,1)
  pre<-r+e
  return(pre)
}

x<-seq(1,8,1)
y<-model.qr(x,-16,1)

jpeg("plot5.jpeg",width=10,height=4,units="in",res=300)
layout(matrix(c(1,2,3,4),nrow=1,byrow=T))
plot(x,y,type="l",col="red",main="β1<0,β2<0")

y<-model.qr(x,4,2)
plot(x,y,type="l",col="green",main="β1>0,β2>0")

y<-model.qr(x,-4,-2)
plot(x,y,type="l",col="blue",main="β1<0,β2<0")

y<-model.qr(x,16,-1)
plot(x,y,type="l",col="yellow",main="β1>0,β2<0")
dev.off()

#no-linear regression model
#x/ax+b
model.1<-function(x,a=2,b=3){
  pre<-x/(a*x+b)
  return(pre)
}

jpeg("plot6.jpeg",width=10,height=4,units="in",res=300)
layout(matrix(c(1,2),nrow=1,byrow=T))
x<-seq(1,8,1)
y<-model.1(x,1,-1)
plot(x,y,type="l",col="red",main="β<0")

y<-model.1(x,1,1)
plot(x,y,type="l",col="red",main="β>0")
dev.off()

#ax^b
model.2<-function(x,a=2,b=3){
  pre<-a*x^b
  return(pre)
}

jpeg("plot7.jpeg",width=10,height=4,units="in",res=300)
layout(matrix(c(1,2),nrow=1,byrow=T))
x<-seq(1,8,1)
y<-model.2(x,2,1)
plot(x,y,type="l",col="red")
y<-model.2(x,2,2)
lines(x,y,col="red")
y<-model.2(x,2,0.5)
lines(x,y,col="red")

y<-model.2(x,2,-1)
plot(x,y,type="l",col="red",ylim=c(0,2))
y<-model.2(x,2,-2)
lines(x,y,col="red")
y<-model.2(x,2,-0.5)
lines(x,y,col="red")
dev.off()

#logx
model.3<-function(x,a=2,b=3){
  pre<-a+b*log(x)
  return(pre)
}

jpeg("plot8.jpeg",width=10,height=4,units="in",res=300)
layout(matrix(c(1,2),nrow=1,byrow=T))
x<-seq(1,8,1)
y<-model.3(x,2,1)
plot(x,y,type="l",col="red",main="β<0")

y<-model.3(x,2,-1)
plot(x,y,type="l",col="red",main="β>0")
dev.off()

#ab^x
model.4<-function(x,a=2,b=3){
  pre<-a^(b*x)
  return(pre)
}

jpeg("plot9.jpeg",width=10,height=4,units="in",res=300)
layout(matrix(c(1,2),nrow=1,byrow=T))
x<-seq(1,8,1)
y<-model.4(x,2,1)
plot(x,y,type="l",col="red",main="β<0")

y<-model.4(x,2,-1)
plot(x,y,type="l",col="red",main="β>0")
dev.off()

#S
model.5<-function(x,a=2,b=3){
  pre<-1/(a+b*exp(-x))
  return(pre)
}

jpeg("plot10.jpeg",width=5,height=4,units="in",res=300)
x<-seq(-8,8,1)
y<-model.5(x,2,1)
#layout(matrix(c(1),nrow=1,byrow=T))
plot(x,y,type="l",col="red",ylim=c(0,0.5))
dev.off()

#logistic regression
x<-seq(1,25,1)
y<-seq(1,25,1)
y[1:12]<-1
y[13:25]<-0
dataxy<-data.frame(x,y)

jpeg("plot11.jpeg",width=10,height=4,units="in",res=300)
layout(matrix(c(1,2),nrow=1,byrow=T))
plot(x,y,pch=16,col="blue")

m.sl<-lm(y~x,data=dataxy)
plot(x,y,pch=16,col="blue")
abline(m.sl,col="red",lty=2)
dev.off()

x<-seq(from=-10,to=10,by = 0.01)
y=exp(x)/(1+exp(x))
jpeg("plot12.jpeg",width=5,height=4,units="in",res=300)
#layout(matrix(c(1),nrow=1,byrow=T))
plot(x,y,type="l",col="blue")
dev.off()