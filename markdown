---
title: "bayesian mixed model"
author: "Zach"
date: "July 26, 2018"
output: html_document
---

#########################IMPORTING DATA

library(readxl)
ThesisSP18Data <- read_excel("~/Graduate School/URI Statistics MS/STA 599 Thesis/Thesis Analyses/ThesisSP18Data.xlsx")

ThesisSP18Data <-as.data.frame(ThesisSP18Data)

#Need to omit final two rows of S18 (they're blank, it's leftover code from a quick analysis in the excel file)

ThesisSP18Data <- ThesisSP18Data[-c(457, 458),]

##Subsetting
S18Var <- c("MA1B", "MA2B", "MA3B", "MA4B", "MA5B", "MA6B", "MA7B", "MA8B","QANX1B", "QANX2B", "QANX3B", "QANX4B", "QUSE1B", "QUSE2B", "QUSE3B", "QUSE4B", "QINFL1B", "QINFL2B", "QINFL3B", "QINFL4B", "QINFL5B", "QINFL6B", "QINFL7B", "QSF1B", "QSF2B", "QSF3B", "QSF4B", "QHIND1B", "QHIND2B", "QHIND3B", "QHIND4B", "QHIND5B", "QSC1B", "QSC2B", "QSC3B", "QSC4B", "QSE1B", "QSE2B", "QSE3B", "QSE3B", "QSE5B", "QSE6B", "SATS1B", "SATS2B", "SATS3B", "SATS4B", "SATS5B", "SATS6B", "SATS7B", "SATS8B", "SATS9B", "SATS10B", "SATS11B", "SATS12B", "SATS13B", "SATS14B", "SATS15B", "SATS16B", "SATS17B", "SATS18B", "SATS19B", "SATS20B", "SATS21B", "SATS22B", "SATS23B", "SATS24B", "SATS25B", "SATS26B", "SATS27B", "SATS28B", "SATS29B", "SATS30B", "SATS31B", "SATS32B", "SATS33B", "SATS34B", "SATS35B", "SATS36B", "MA1M", "MA2M", "MA3M", "MA4M", "MA5M", "MA6M", "MA7M", "MA8M", "QANX1M", "QANX2M", "QANX3M", "QANX4M", "QUSE1M", "QUSE2M", "QUSE3M", "QUSE4M", "QINFL1M", "QINFL2M", "QINFL3M", "QINFL4M", "QINFL5M", "QINFL6M", "QINFL7M", "QSF1M", "QSF2M", "QSF3M", "QSF4M", "QHIND1M", "QHIND2M", "QHIND3M", "QHIND4M", "QHIND5M", "QSC1M", "QSC2M", "QSC3M", "QSC4M", "QSE1M", "QSE2M", "QSE3M", "QSE3M", "QSE5M", "QSE6M", "SATS1M", "SATS2M", "SATS3M", "SATS4M", "SATS5M", "SATS6M", "SATS7M", "SATS8M", "SATS9M", "SATS10M", "SATS11M", "SATS12M", "SATS13M", "SATS14M", "SATS15M", "SATS16M", "SATS17M", "SATS18M", "SATS19M", "SATS20M", "SATS21M", "SATS22M", "SATS23M", "SATS24M", "SATS25M", "SATS26M", "SATS27M", "SATS28M", "SATS29M", "SATS30M", "SATS31M", "SATS32M", "SATS33M", "SATS34M", "SATS35M", "SATS36M", "MA1E", "MA2E", "MA3E", "MA4E", "MA5E", "MA6E", "MA7E", "MA8E","QANX1E", "QANX2E", "QANX3E", "QANX4E", "QUSE1E", "QUSE2E", "QUSE3E", "QUSE4E", "QINFL1E", "QINFL2E", "QINFL3E", "QINFL4E", "QINFL5E", "QINFL6E", "QINFL7E", "QSF1E", "QSF2E", "QSF3E", "QSF4E", "QHIND1E", "QHIND2E", "QHIND3E", "QHIND4E", "QHIND5E", "QSC1E", "QSC2E", "QSC3E", "QSC4E", "QSE1E", "QSE2E", "QSE3E", "QSE3E", "QSE5E", "QSE6E", "SATS1E", "SATS2E", "SATS3E", "SATS4E", "SATS5E", "SATS6E", "SATS7E", "SATS8E", "SATS9E", "SATS10E", "SATS11E", "SATS12E", "SATS13E", "SATS14E", "SATS15E", "SATS16E", "SATS17E", "SATS18E", "SATS19E", "SATS20E", "SATS21E", "SATS22E", "SATS23E", "SATS24E", "SATS25E", "SATS26E", "SATS27E", "SATS28E", "SATS29E", "SATS30E", "SATS31E", "SATS32E", "SATS33E", "SATS34E", "SATS35E", "SATS36E", "Resil1B", "Resil2B", "Resil3B", "Resil4B", "Resil5B", "Resil6B", "Grit1B", "Grit2B", "Grit3B", "Grit4B", "Grit5B", "Grit6B", "Grit7B", "Grit8B", "Resil1M", "Resil2M", "Resil3M", "Resil4M", "Resil5M", "Resil6M", "Grit1M", "Grit2M", "Grit3M", "Grit4M", "Grit5M", "Grit6M", "Grit7M", "Grit8M", "Resil1E", "Resil2E", "Resil3E", "Resil4E", "Resil5E", "Resil6E", "Grit1E", "Grit2E", "Grit3E", "Grit4E", "Grit5E", "Grit6E", "Grit7E", "Grit8E")

S18<-ThesisSP18Data[S18Var]

library(MissMech)
TestMCARNormality(S18)  
###Error in solve.default(sigoo) : system is computationally singular: reciprocal condition number = 1.03251e-17
##Usually these types of error come from linear dependence, which makes sense since the data are longitudinal.
##Could proceed by doing data reduction (factor scores/averages), then testing for MCAR/MVN on the reduced dataset?

#################CODE COPY/PASTED FROM YUAN & MACKINNON (2009)

model
{
# specify the multilevel model
# N1 and N2 are the number of first-level units and second-level units,
respectively.
for( j in 1:N2){
for(i in 1:N1){
# Specify the first-level model: Mij = β2j + αjXij + e2ij

m[j, i] ~ dnorm(mean.m[j, i], prec.m)
mean.m[j, i] < − Beta2[j] + AlphaBeta[j, 1]*x[j, i]
# Specify the first-level model: Yij = β3j + βjMij + τ′jXij + e3ij
y[j, i] ~ dnorm(mean.y[j, i], prec.y)
mean.y[j, i] < − Beta3[j] + AlphaBeta[j, 2]*m[j,i] + Tau.p[j]*x[j,
i]
}
# Specify the second-level models
Beta2[j] ~ dnorm(beta2, prec.beta2)
Beta3[j] ~ dnorm(beta3, prec.beta3)
Tau.p[j] ~ dnorm(tau.p, prec.taup)
# bivariate normal distribution for αi and β i
AlphaBeta[j, 1:2] ~ dmnorm(alphabeta[ ], prec.ab[,])
}
# vague bivariate normal prior for α and β
alphabeta[1:2] ~ dmnorm(mean[ ], prec[,])
mean[1] < − 0
mean[2] < − 0
prec[1,1] < − 1.0E-6
prec[1,2] < − 0
prec[2,1] < − 0
prec[2,2] < − 1.0E-6
# vague inverse-Wishart prior for the covariance of αj and β j;
# dwish(·) is the Wishart distribution
prec.ab[1:2 , 1:2] ~ dwish(Omega[,], 2)
Omega[1,1] < − 0.001
Omega[1,2] < − 0.0
Omega[2,1] < − 0.0
Omega[2,2] < − 0.001
# vague normal priors for β2, β3 and τ′
beta2 ~ dnorm(0.0, 1.0E-6)
beta3 ~ dnorm(0.0, 1.0E-6)
tau.p ~ dnorm(0.0, 1.0E-6)
# vague inverse-gamma prior for variances of the first level model
prec.y ~ dgamma(0.001, 0.001)
prec.m ~ dgamma(0.001, 0.001)
# vague uniform priors for standard deviations of the second level model
sigma.beta2 ~ dunif(0, 100)
sigma.taup ~ dunif(0, 100)
sigma.beta3 ~ dunif(0, 100)
prec.beta2 < − 1/(sigma.beta2*sigma.beta2)
prec.taup < − 1/(sigma.taup*sigma.taup)
prec.beta3 < − 1/(sigma.beta3*sigma.beta3)

## define quantities of interest
# covariance matrix of αi and βi
var.ab[1,1] < − prec.ab[2,2]/(prec.ab[1,1]*prec.ab[2,2]-prec.ab[1,2]
*prec.ab[2,1])
var.ab[1,2] < − -prec.ab[1,2]/(prec.ab[1,1]*prec.ab[2,2]-prec.ab[1,2]
*prec.ab[2,1])
var.ab[2,1] < − -prec.ab[2,1]/(prec.ab[1,1]*prec.ab[2,2]-prec.ab[1,2]
*prec.ab[2,1])
var.ab[2,2] < − prec.ab[1,1]/(prec.ab[1,1]*prec.ab[2,2]-prec.ab[1,2]
*prec.ab[2,1])
# average mediated effect ab, total effect c and relative mediated effect
r = ab/c
ab < − alphabeta[1]*alphabeta[2]+ var.ab[1,2]
c < − alphabeta[1]*alphabeta[2]+ var.ab[1,2]+tau.p
r < − ab/c
}
