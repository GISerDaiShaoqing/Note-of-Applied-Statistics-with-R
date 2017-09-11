# Name: Note of Applied Statistics with R(Chapter 11 demo code) Discriminant Analysis
# Purpose: Discriminant Analysis
# Author:      Dai shaoqing
#
# Created:     09/11/2017
# Copyright:   (c) Dai shaoqing <dsq1993qingge@163.com> 2017
#------------------------------------------------------------
#Load packages
library(raster)
library(WMDB)
library(MASS)
library(maptools)
library(klaR)

#Generate false color image
aplot<-brick(".../GF1test.tif")
jpeg("rawImage.jpeg",width=6,height=7.2,units="in",res=300)
plotRGB(aplot,r=4,g=3,b=2)
dev.off()

#Plot each band in image
jpeg("singleband.jpeg",width=6,height=7.2,units="in",res=300)
spplot(aplot)
dev.off()

#Plot train data
traindatasp<-readShapePoints(".../trainGF.shp")
jpeg("traindata.jpeg",width=6,height=7.2,units="in",res=300)
par(new=T)
layout(matrix(seq(1,4,1),nrow=2,byrow=T))
plotRGB(aplot,r=4,g=3,b=2)
plot(traindatasp,col=traindatasp$type,pch=16,add=TRUE)
dev.off()

#Plot validation data
validationsp<-readShapePoints(".../validation.shp")
jpeg("validation.jpeg",width=6,height=7.2,units="in",res=300)
plotRGB(aplot,r=4,g=3,b=2)
plot(validationsp,col=validationsp$type,pch=16,add=TRUE)
dev.off()

#Pre-processing
traindata<-as.data.frame(traindatasp)
traindata<-traindata[,-c(6,7)]
classification<-as.factor(traindata[,1])
bandcl<-traindata[,-1]
validation<-as.data.frame(validationsp)
validation<-validation[,-c(6,7)]
validationcl<-as.factor(validation[,1])
validationdata<-validation[,-1]
imagedata<-as.data.frame(aplot)
imagedata[is.na(imagedata)]<-0

#Distance
wmd(bandcl,classification)

#Fisher
traindatanew<-data.frame(traindata[,-1],classification)
z<-lda(classification~Band_1+Band_2+Band_3+Band_4,traindatanew)
z
new=predict(z,bandcl)$class
table(new,classification)
49/55

#Bayes
dbayes(bandcl,classification)

#Raster classification(Fisher)
imagepre<-data.frame(Band_1=imagedata$GF1test.1,Band_2=imagedata$GF1test.2,Band_3=imagedata$GF1test.3,Band_4=imagedata$GF1test.4)
classresultimage<-predict(z,imagepre)$class

refc<-crs(a)
extc<-extent(a)
colc<-ncol(a)
rowc<-nrow(a)
resc<-res(a)
maskc<-a

cnewraster<-raster(ncol=colc,nrow=rowc,ext=extc,crs=refc,resolution=resc)
cnewraster[]<-classresultimage
craster<-mask(cnewraster,maskc)

#Raster classification(Bayes)
Bayesmodel<-NaiveBayes(classification~Band_1+Band_2+Band_3+Band_4,data=traindatanew)
bayesimage<-predict(Bayesmodel,imagepre)

bnewraster<-raster(ncol=colc,nrow=rowc,ext=extc,crs=refc,resolution=resc)
bnewraster[]<-bayesimage$class
braster<-mask(bnewraster,maskc)

jpeg('twoclassification.jpeg',width=10,height=10,units="in",res=300)
layout(matrix(seq(1,4,1),nrow=2,byrow=T))
plotRGB(aplot,r=4,g=3,b=2)
plot(craster,xaxt="n",yaxt="n")
plotRGB(aplot,r=4,g=3,b=2)
plot(braster,xaxt="n",yaxt="n")
dev.off()

#Raster Classsification output
writeRaster(craster,'Fisher.tif',format="GTiff")
writeRaster(braster,'Bayes.tif',format="GTiff")

#Accuracy
accuracysp<-readShapePoints(".../validationnew.shp")
accuracy<-as.data.frame(accuracysp)
table(accuracy$type,accuracy$Bayes)
table(accuracy$type,accuracy$Fisher)
