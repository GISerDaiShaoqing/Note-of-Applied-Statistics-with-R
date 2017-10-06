# Name: Note of Applied Statistics with R(Chapter 13 demo code) Factor Analysis
# Purpose: Factor Analysis
# Author:      Dai shaoqing
#
# Created:     10/06/2017
# Copyright:   (c) Dai shaoqing <dsq1993qingge@163.com> 2017
#------------------------------------------------------------
#Load packages
library(psych)

#Input data
a<-read.csv('.../city2006.csv',header=T)
aclean<-na.omit(a)

#Function1 of Factor Analysis
factanal(aclean[,-c(1:2,9)],3,data=aclean[,-c(1:2,9)])
afa<-factanal(aclean[,-c(1:2,9)],1,data=aclean[,-c(1:2,9)])

#Function2 of Factor Analysis
factor.analysis=function(x,m){
  p=nrow(x);x.diag=diag(x);sum.rank=sum(x.diag)
  rowname=paste("X",1:p,sep="")
  colname=paste("Factor",1:m,sep="")
  
  #Construct loading matrix
  A=matrix(0,nrow=p,ncol=m,dimnames=list(rowname,colname))
  eig=eigen(x)#The 'eig' includes two elements,'values' is characteristic root, 'vectors' is characteristic vector
  for (i in 1:m){
    A[,i]=sqrt(eig$values[i])*eig$vectors[,i]#Fill A matrix
  }
  var.A=diag(A%*%t(A))#Variance of Factor
  rowname1=c('SS loadings','ProportationVar','Cumulative Var')
  #Construct result matrix
  result=matrix(0,nrow=3,ncol=m,dimnames=list(rowname1,colname))
  for (i in 1:m) {
    result[1,i]=sum(A[,i]^2)#Variance of each factor
    result[2,i]=result[1,i]/sum.rank#variance contribution
    result[3,i]=sum(result[1,1:i])/sum.rank#Cumulative variacnce contribution
  }
  method=c("Principal Component Method")
  #Print result
  list(method=method,loadings=A,var=cbind(common=var.A,specific=x.diag-var.A),result=result)
}
R=cor(aclean[,-c(1,2,9)])
factor.analysis(R,3)

#Function3 of Factor Analysis
afa<-fa(aclean[,-c(1:2,9)],2)
fa(aclean[,-c(1:2,9)],3,fm="mle")
afa<-fa(aclean[,-c(1:2,9)],3,fm="mle")

jpeg("1.jpeg",width=6,height=6,units="in",res=300)
fa.diagram(afa)
dev.off()

jpeg("2.jpeg",width=6,height=6,units="in",res=300)
biplot(afa)
dev.off()