# Name: Note of Applied Statistics with R(Chapter 12 demo code) Priciple Component Analysis
# Purpose: Priciple Component Analysis
# Author:      Dai shaoqing
#
# Created:     09/21/2017
# Copyright:   (c) Dai shaoqing <dsq1993qingge@163.com> 2017
#------------------------------------------------------------
#Load packages
library(psych)
library(corrplot)

#Input data
a<-read.csv('.../city2006.csv',header=T)
aclean<-na.omit(a)

#Plot the correlation
acor<-corr.test(aclean[,-1])
acorp<-acor$p
acorp[upper.tri(acorp)]=0
jpeg("1.jpeg",width=6,height=6,units="in",res=300)
par(fig=c(0,1,0.1,1))
corrplot.mixed(acor$r,upper="square",lowe="number",
               diag="u",tl.cex=1,tl.col="black",tl.pos="lt",number.cex=0.8,
               cl.cex=1,p.mat=acorp, sig.level=0.05,insig=c("blank"))
dev.off()

#Function1 of PCA
#Base covirance matrix
apca1<-princomp(aclean[,-c(1:2,9)])

#Screenplot
jpeg("2.jpeg",width=6,height=6,units="in",res=300)
screeplot(apca1,type='line')
dev.off()

#Result of PCA
summary(apca1,loadings=T)

#Biplot of PCA
jpeg("3.jpeg",width=6,height=6,units="in",res=300)
biplot(apca1)
dev.off()

#Base correlation matrix
apca1<-princomp(aclean[,-c(1:2,9)],cor=T)

#Screenplot
jpeg("4.jpeg",width=6,height=6,units="in",res=300)
screeplot(apca1,type='line')
dev.off()

#Result of PCA
summary(apca1,loadings=T)

#Biplot of PCA
jpeg("5.jpeg",width=6,height=6,units="in",res=300)
biplot(apca1)
dev.off()

#Function2 of PCA
#Screenplot
jpeg("6.jpeg",width=6,height=6,units="in",res=300)
fa.parallel(aclean[,-c(1:2,9)],fa="pc",n.iter = 100,show.legend = T,main="Scree plot with Parallel analysis")
dev.off()

#PCA function
apca2<-principal(aclean[,c(3:8,10,11,12)],nfactors=2,rotate="varimax",scores=T)

#Diagram of PCA
jpeg("7.jpeg",width=6,height=6,units="in",res=300)
fa.diagram(apca2)
dev.off()

#Biplot of PCA
jpeg("8.jpeg",width=6,height=6,units="in",res=300)
biplot(apca2)
dev.off()