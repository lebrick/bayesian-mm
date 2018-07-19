#################IMPORTING DATA

library(readxl)

#################REVERSE SCORING

#################SUBSETTING DATA

#################MCAR AND MVN TESTING

library(MissMech)
TestMCARNormality(data)

library(BaylorEdPsych)
LittleMCAR(data)

library(MVN)
mvn(data)

#################FACTOR ANALYSIS, MAP TEST, AND OMEGA (and alpha too)
##probably don't need to report MAP test, all measures should be single factors

library(psych) ##note sample code is below is from MS thesis

F17QAnxB <- c("QANX1B", "QANX2B", "QANX3B", "QANX4B")
F17QAnxBOCC <-F17CC[F17QAnxB]
omega(F17QAnxBOCC, fm="ml") #omega = .91, alpha = .86

VSS(data)

QANX_FAB <-fa(F17QAnxBOCC, factors=1,rotation="promax",scores="regression")
QANX_FSB <- as.data.frame(QANX_FAB$scores)
F17CC <- cbind(F17CC, QANX_FSB)
F17CC$QANX_FSB <- as.numeric(unlist(QANX_FSB))

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
