# Name: Odds ratio calculation
# Purpose: OR
# Author:      Dai shaoqing
#
# Created:     03/16/2021
# Copyright:   (c) Dai shaoqing <dsq1993qingge@163.com> 2021
#------------------------------------------------------------
## Install the r package namely aplore3 which provides the datasets.
install.packages('aplore3')
install.packages('dplyr')

## Load package
library(aplore3)
library(dplyr)

data(icu)
summary(icu)
head(icu)

plot(icu$age, icu$sta, pch = 16, col = 'red', xlab = 'age', ylab = 'outcome')

## Pre-processing
icu$agegroup <- ifelse(icu$age > 0 & icu$age <= 20, 0,
                       ifelse(icu$age > 20 & icu$age <= 30, 1,
                              ifelse(icu$age > 30 & icu$age <= 40, 2,
                                     ifelse(icu$age > 40 & icu$age <= 50, 3,
                                            ifelse(icu$age > 50 & icu$age <= 60, 4,
                                                   ifelse(icu$age > 60 & icu$age <= 70, 5,
                                                          ifelse(icu$age > 70 & icu$age <= 80, 6,
                                                                 ifelse(icu$age > 80 & icu$age <= 90, 7, 8))))))))

icu$deathc <- ifelse(icu$sta=='Lived', 0, 1)

icuagegroup <- icu %>%
  group_by(agegroup) %>%
  summarise(pr = sum(deathc))

groupn <- icu %>%
  group_by(agegroup) %>%
  count()

icuagegroup$pr <- icuagegroup$pr/groupn$n

plot(icuagegroup, pch = 16, col = 'red', xlab = 'age group', ylab = 'pr(death)')

## Display the p-logit(p) plot
p <- runif(10000, min = 0, 1)

logitf <- function(p) {
  logp <- log(p/(1-p))
  return(logp)
}

logitp <- logitf(p)

plpdf <- data.frame(p, logitp)
plpdf <- plpdf[order(plpdf$p),]

plot(plpdf$p, plpdf$logitp, type = 'l', col = 'red', xlab = 'p', ylab = 'logit(p)')

## OR calculate
## 1 without any packages
modellogit <- glm(sta~type, data = icu, family = binomial)
ORDF <- data.frame(exp(cbind(OR = coef(modellogit), confint(modellogit))))
ORDF

## 2 Using epiDisplay
install.packages('epiDisplay')
library(epiDisplay)
modellogit <- glm(sta~type, data = icu, family = binomial)
ORDF <- logistic.display(modellogit)
ORDF

## 3 Using questionr
install.packages('questionr')
library(questionr)
modellogit <- glm(sta~type, data = icu, family = binomial)
ORDF <- odds.ratio(modellogit)
ORDF

## Set reference group
icu$typen <- relevel(icu$type, ref = "Emergency")
modellogit <- glm(sta~typen, data = icu, family = binomial)
ORDF <- odds.ratio(modellogit)
ORDF
