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

#####################REVERSE SCORING

S18BMM$Resil2B <- 6 - S18BMM$Resil2B
S18BMM$Resil4B <- 6 - S18BMM$Resil4B
S18BMM$Resil6B <- 6 - S18BMM$Resil6B

S18BMM$Resil2M <- 6 - S18BMM$Resil2M
S18BMM$Resil4M <- 6 - S18BMM$Resil4M
S18BMM$Resil6M <- 6 - S18BMM$Resil6M

S18BMM$Resil2E <- 6 - S18BMM$Resil2E
S18BMM$Resil4E <- 6 - S18BMM$Resil4E
S18BMM$Resil6E <- 6 - S18BMM$Resil6E

S18BMM$Grit1B <- 6 - S18BMM$Grit1B
S18BMM$Grit3B <- 6 - S18BMM$Grit3B
S18BMM$Grit5B <- 6 - S18BMM$Grit5B
S18BMM$Grit6B <- 6 - S18BMM$Grit6B

S18BMM$Grit1M <- 6 - S18BMM$Grit1M
S18BMM$Grit3M <- 6 - S18BMM$Grit3M
S18BMM$Grit5M <- 6 - S18BMM$Grit5M
S18BMM$Grit6M <- 6 - S18BMM$Grit6M

S18BMM$Grit1E <- 6 - S18BMM$Grit1E
S18MMM$Grit3E <- 6 - S18BMM$Grit3E
S18BMM$Grit5E <- 6 - S18BMM$Grit5E
S18BMM$Grit6E <- 6 - S18BMM$Grit6E

#####################CREATING FACTOR SCORES

QAnxB <- c("QANX1B", "QANX2B", "QANX3B", "QANX4B")
QAnxBFA <-S18BMM[QAnxB]
QANX_FAB <-fa(QAnxBFA, factors=1,rotation="promax",scores="regression")
QANX_FSB <- as.data.frame(QANX_FAB$scores)
S18BMM <- cbind(S18BMM, QANX_FSB)
S18BMM$QANX_FSB <- as.numeric(unlist(QANX_FSB))

QAnxM <- c("QANX1M", "QANX2M", "QANX3M", "QANX4M")
QAnxMFA <-S18BMM[QAnxM]
QANX_FAM <-fa(QAnxMFA, factors=1,rotation="promax",scores="regression")
QANX_FSM <- as.data.frame(QANX_FAM$scores)
S18BMM <- cbind(S18BMM, QANX_FSM)
S18BMM$QANX_FSM <- as.numeric(unlist(QANX_FSM))

QAnxE <- c("QANX1E", "QANX2E", "QANX3E", "QANX4E")
QAnxEFA <-S18BMM[QAnxE]
QANX_FAE <-fa(QAnxEFA, factors=1,rotation="promax",scores="regression")
QANX_FSE <- as.data.frame(QANX_FAE$scores)
S18BMM <- cbind(S18BMM, QANX_FSE)
S18BMM$QANX_FSE <- as.numeric(unlist(QANX_FSE))

QSEB <- c("QSE1B", "QSE2B", "QSE3B", "QSE4B", "QSE5B", "QSE6B")
QSEBFA <-S18BMM[QSEB]
QSE_FAB <-fa(QSEBFA, factors=1,rotation="promax",scores="regression")
QSE_FSB <- as.data.frame(QSE_FAB$scores)
S18BMM <- cbind(S18BMM, QSE_FSB)
S18BMM$QSE_FSB <- as.numeric(unlist(QSE_FSB))

QSEM <- c("QSE1M", "QSE2M", "QSE3M", "QSE4M")
QSEMFA <-S18BMM[QSEM]
QSE_FAM <-fa(QSEMFA, factors=1,rotation="promax",scores="regression")
QSE_FSM <- as.data.frame(QSE_FAM$scores)
S18BMM <- cbind(S18BMM, QSE_FSM)
S18BMM$QSE_FSM <- as.numeric(unlist(QSE_FSM))

QSEE <- c("QSE1E", "QSE2E", "QSE3E", "QSE4E")
QSEEFA <-S18BMM[QSEE]
QSE_FAE <-fa(QSEEFA, factors=1,rotation="promax",scores="regression")
QSE_FSE <- as.data.frame(QSE_FAE$scores)
S18BMM <- cbind(S18BMM, QSE_FSE)
S18BMM$QSE_FSE <- as.numeric(unlist(QSE_FSE))

ResilB <- c("Resil1B", "Resil2B", "Resil3B", "Resil4B", "Resil5B", "Resil6B")
ResilBFA <- S18BMM[ResilB]
R_FAB <-fa(ResilBFA, factors=1, scores="regression")
R_FSB <-as.data.frame(R_FAB$scores)
S18BMM <- cbind(S18BMM, R_FSB)
S18BMM$R_FSB <- as.numeric(unlist(R_FSB))

ResilM <- c("Resil1M", "Resil2M", "Resil3M", "Resil4M", "Resil5M", "Resil6M")
ResilMFA <- S18BMM[ResilM]
R_FAM <-fa(ResilMFA, factors=1, scores="regression")
R_FSM <-as.data.frame(R_FAM$scores)
S18BMM <- cbind(S18BMM, R_FSM)
S18BMM$R_FSM <- as.numeric(unlist(R_FSM))

ResilE <- c("Resil1E", "Resil2E", "Resil3E", "Resil4E", "Resil5E", "Resil6E")
ResilEFA <- S18BMM[ResilE]
R_FAE <-fa(ResilEFA, factors=1, scores="regression")
R_FSE <-as.data.frame(R_FAE$scores)
S18BMM <- cbind(S18BMM, R_FSE)
S18BMM$R_FSE <- as.numeric(unlist(R_FSE))

GritBF <- c("Grit1B", "Grit2B", "Grit3B", "Grit4B", "Grit5B", "Grit6B", "Grit7B", "Grit8B")
GritBFA <- S18BMM[GritB]
G_FAB <-fa(GritBFA, nfactors=2, rotate = "promax", scores="tenBerge")
G_FSB <-as.data.frame(G_FAB$scores)
S18BMM <- cbind(S18BMM, G_FSB)
S18BMM$G_FSB <- as.numeric(unlist(G_FSB))

GritM <- c("Grit1M", "Grit2M", "Grit3M", "Grit4M", "Grit5M", "Grit6M", "Grit7M", "Grit8M")
GritMFA <- S18BMM[GritM]
G_FAM <-fa(GritMFA, nfactors=2, rotate = "promax", scores="tenBerge")
G_FSM <-as.data.frame(G_FAM$scores)
S18BMM <- cbind(S18BMM, G_FSM)
S18BMM$G_FSM <- as.numeric(unlist(G_FSM))

GritE <- c("Grit1E", "Grit2E", "Grit3E", "Grit4E", "Grit5E", "Grit6E", "Grit7E", "Grit8E")
GritEFA <- S18BMM[GritE]
G_FAE <-fa(GritEFA, nfactors=2, rotate="promax", scores="tenBerge")
G_FSE <-as.data.frame(G_FAE$scores)
S18BMM <- cbind(S18BMM, G_FSE)
S18BMM$G_FSE <- as.numeric(unlist(G_FSE))

#Grit is a 2-factor structure (F1 = items 2, 4, 7, 8; F2 = items 1, 3, 5, 6)
#was a pain to extract two factors from same analysis
#just did them separaretly

GritBF1 <- c("Grit2B", "Grit4B", "Grit7B", "Grit8B")
GritBF2 <- c("Grit1B", "Grit3B", "Grit5B", "Grit6B")
GritBFA1 <- S18BMM[GritBF1]
GritBFA2 <- S18BMM[GritBF2]
G_FAB1 <- fa(GritBFA1, nfactors=1, scores="regression")
G_FAB2 <- fa(GritBFA2, nfactors=1, scores="regression")
G_FSB1 <-as.data.frame(G_FAB1$scores)
G_FSB2 <-as.data.frame(G_FAB2$scores)
S18BMM <- cbind(S18BMM, G_FSB1)
S18BMM <- cbind(S18BMM, G_FSB2)
S18BMM$G_FSB1 <- as.numeric(unlist(G_FSB1))
S18BMM$G_FSB2 <- as.numeric(unlist(G_FSB2))

GritMF1 <- c("Grit2M", "Grit4M", "Grit7M", "Grit8M")
GritMF2 <- c("Grit1M", "Grit3M", "Grit5M", "Grit6M")
GritMFA1 <- S18BMM[GritMF1]
GritMFA2 <- S18BMM[GritMF2]
G_FAM1 <- fa(GritMFA1, nfactors=1, scores="regression")
G_FAM2 <- fa(GritMFA2, nfactors=1, scores="regression")
G_FSM1 <-as.data.frame(G_FAM1$scores)
G_FSM2 <-as.data.frame(G_FAM2$scores)
S18BMM <- cbind(S18BMM, G_FSM1)
S18BMM <- cbind(S18BMM, G_FSM2)
S18BMM$G_FSM1 <- as.numeric(unlist(G_FSM1))
S18BMM$G_FSM2 <- as.numeric(unlist(G_FSM2))

GritEF1 <- c("Grit2E", "Grit4E", "Grit7E", "Grit8E")
GritEF2 <- c("Grit1E", "Grit3E", "Grit5E", "Grit6E")
GritEFA1 <- S18BMM[GritEF1]
GritEFA2 <- S18BMM[GritEF2]
G_FAE1 <- fa(GritEFA1, nfactors=1, scores="regression")
G_FAE2 <- fa(GritEFA2, nfactors=1, scores="regression")
G_FSE1 <-as.data.frame(G_FAE1$scores)
G_FSE2 <-as.data.frame(G_FAE2$scores)
S18BMM <- cbind(S18BMM, G_FSE1)
S18BMM <- cbind(S18BMM, G_FSE2)
S18BMM$G_FSE1 <- as.numeric(unlist(G_FSE1))
S18BMM$G_FSE2 <- as.numeric(unlist(G_FSE2))


#################Subsetting for Factor Scores

FactorScores <- c("QANX_FSB", "QANX_FSM", "QANX_FSE", "QSE_FSB", "QSE_FSM", "QSE_FSE", "R_FSB", "R_FSM", "R_FSE", "G_FSB1", "G_FSB2", "G_FSM1", "G_FSM2", "G_FSE1", "G_FSE2")
FSBMM <- S18BMM[FactorScores]

TestMCARNormality(FSBMM)

#Data are not MCAR or MVN. Darn.

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
