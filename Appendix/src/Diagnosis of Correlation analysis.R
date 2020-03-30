# Name: Diagnosis of Correlation analysis
# Purpose: Correlation Analysis
# Author:      Dai shaoqing
#
# Created:     03/30/2020
# Copyright:   (c) Dai shaoqing <dsq1993qingge@163.com> 2020
#------------------------------------------------------------
#Load packages
library(psych)
library(corrplot)

#Generate the random dataset
at <- cbind(runif(41, -1, 1), runif(41, -2, 2), runif(41, -3, 3), runif(41, -5, 5))

#Convert the matrix to data frame
at <- data.frame(at)
colnames(at) <- c("a", "b", "c", "d")

#Set one element to NA
at[6, 4] <- NA
at

#Correlation analysis under the different parameters 
corr.test(at[, c(2,3)], use = "pairwise", adjust = 'none')
corr.test(at[, c(2:4)], use = "pairwise", adjust = 'none')

corr.test(at[, c(2,3)], use = "complete", adjust = 'none')
corr.test(at[, c(2:4)], use = "complete", adjust = 'none')

corr.test(at[, c(2,3)], use = "complete.obs", adjust = 'none')
corr.test(at[, c(2:4)], use = "complete.obs", adjust = 'none')

#------------------------------------------------------------
#The Lewis's Case
x <- matrix(c(-2,-1,0,1,2,1.5,2,0,1,2,NA,NA,0,1,2),5)

#Convert the matrix to data frame
xu <- data.frame(x)
colnames(xu) <- c("a", "b", "c")

#Correlation analysis under the different parameters 
cor(x, use = "everything")
cor(x, use = "pairwise")
cor(x, use = "complete.obs")
cor(x, use = "pairwise.complete.obs")

corr.test(xu, use = "everything", , adjust = 'none')
corr.test(xu, use = "pairwise", adjust = 'none')
corr.test(xu, use = "complete.obs", adjust = 'none')
corr.test(xu, use = "pairwise.complete.obs", adjust = 'none')

#Plot the data and calculate the covriance
plot(xu[,1], xu[,2], xlab = "a", ylab = "b", col = 'red', pch = 16)
cov(xu[,1], xu[,2])

#Used corrplot package visualize the correlation analysis under the different parameters.
c1m <- corr.test(xu, use = "everything", , adjust = 'none')
c2m <- corr.test(xu, use = "pairwise", adjust = 'none')
c3m <- corr.test(xu, use = "complete.obs", adjust = 'none')

layout(mat = matrix(c(1,2,3), nrow = 1, byrow = T))
corrplot.mixed(c1m$r, upper = "number", lower = "circle")
corrplot.mixed(c2m$r, upper = "number", lower = "circle")
corrplot.mixed(c3m$r, upper = "number", lower = "circle")
