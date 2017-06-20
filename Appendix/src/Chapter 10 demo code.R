# Name: Note of Applied Statistics with R(Chapter 10 demo code) Cluster Analysis
# Purpose: Cluster Analysis
# Author:      Dai shaoqing
#
# Created:     06/20/2017
# Copyright:   (c) Dai shaoqing <dsq1993qingge@163.com> 2017
#------------------------------------------------------------
#load packages
library(reshape)
library(fpc)
library(cluster)
library(sp)
library(maptools)
library(leaflet)

#preprocessing data
aqidata<-read.csv(".../air20170607.csv",header=T)
airdata<-melt(aqidata[,-c(1,3,12)],id=c("stationname","time"))
airdata<-airdata[c(289:576),]
airnew<-cast(airdata,stationname~time+variable)

#hierarchical clustering method
#Caculate Distance and visualization
stationname<-airnew[,1]
airclu<-airnew[,-1]
rownames(airclu)<-stationname
dist.pm25<-dist(airnew[,-1],method='euclidean')
jpeg("plot1.jpeg",width=12,height=8,units="in",res=300)
heatmap(as.matrix(dist.pm25),labRow=stationname,labcol=F)
dev.off()

#Cluster
model1=hclust(dist.pm25,method="ward")

#Visualization
#Choose one from the two lines below
jpeg("plot2.jpeg",width=12,height=8,units="in",res=300)
plot(model1,labels=stationname,hang=-1,las=1)
plclust(model1,labels=stationname,hang=-1)
dev.off()

#result output
result=cutree(model1,k=3)

#result visualization
jpeg("plot3.jpeg",width=12,height=8,units="in",res=300)
plot(airnew[,2],airnew[,3],col=result,pch=as.integer(result))
dev.off()

#K-means clustering method
kres<-kmeans(airnew[,-1],centers=3,nstart=10)

#result visualization
jpeg("plot4.jpeg",width=12,height=8,units="in",res=300)
plotcluster(airnew[,-1],kres$cluster)
dev.off()

jpeg("plot5.jpeg",width=12,height=8,units="in",res=300)
layout(matrix(seq(1,2,1),nrow=1,byrow=T))
plot(airnew[,2],airnew[,3],col=result,pch=as.integer(result))
plot(airnew[,2],airnew[,3],col=kres$cluster,pch=as.integer(kres$cluster))
dev.off()

#station
#data pre-processing
airclres<-data.frame(stationname=airnew[,1],cluster1=result,cluster2=kres$cluster)
stationlocation<-read.table(".../station.txt",header=T,sep=",")
airsta<-stationlocation
airsta<-merge(airsta,airclres,by="stationname")
coordinates(airsta)=~lng+lat

#Construct Interaccion Map
leaflet(airsta)%>%addTiles()%>%addProviderTiles("OpenStreetMap.Mapnik")%>%
  addMarkers(airsta,lng=airsta$lng,lat=airsta$lat,label=airsta$stationname)

spplot(airsta,zcol="cluster1")