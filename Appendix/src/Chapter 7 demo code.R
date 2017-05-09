# Name: Note of Applied Statistics with R(Chapter 7 demo code) Goodness of Fit
# Purpose: Goodness of Fit
# Author:      Dai shaoqing
#
# Created:     05/10/2017
# Copyright:   (c) Dai shaoqing <dsq1993qingge@163.com> 2017
#------------------------------------------------------------
#If you have alreadly installed the packages, please comment.
#install.packages(plyr) 

#Goodness of Fit about Multinomial distribution
goftestmp<-function(x,e,alpha=0.05){
  library(plyr)
  cn=colnames(x)
  countx=count(x,var=cn[1])
  k=nrow(countx)
  chi=sum((countx$freq-e)^2/e)
  chit=qchisq(1-alpha,k-1)
  cat("Hypothesis is :",ifelse(abs(chi)>chit,"Rejected","Accepted"))
}

#QQ plot
#generation of random number that fall in normal distribution
a<-rnorm(200,0,1)

#plot
jpeg("plot1.jpg",width = 5000,height = 4000,units = "px",res = 1000)
qqnorm(a)
qqline(a,col="red")
dev.off()