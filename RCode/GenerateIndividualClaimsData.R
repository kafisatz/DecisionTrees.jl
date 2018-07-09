################################################################
### Individual Claims History Simulation Machine
### Source: https://people.math.ethz.ch/~wmario/simulation.html
### This code is based on the work of M. Wuethrich and A. Gabrielli, see above
### Author: Bernhard Koenig
### Date: July 1, 2018

### Original header of the code:
################################################################
### Demonstration of the Simulation Machine:                 ###
### Generates individual accident insurance claims           ###
### including reporting pattern, cash flow pattern           ###
### and claims closing information.                          ###
### Original Authors: Andrea Gabrielli, Mario V. Wuthrich    ###
### Date: 26.02.2018                                         ###
### Version: 1                                               ###
### Reference: Individual Claims History Simulation Machine. ###
###            SSRN Manuscript ID 3130560.                   ###  
################################################################


##################################
### Load the required packages ###
##################################

library(parallel)
library(foreach)
library(doParallel)
library(data.table)
library(plyr)
library(MASS)
library(ChainLadder)


#######################################################################
### Read in the functions Simulation.Machine and Feature.Generation ###
#######################################################################


### IMPORTANT: Change the working directory to the directory containing the source files for the simulation machine:
### Session --> Set Working Directory --> To Source File Location
setwd("~/ASync/home/Code/Julia/Ressources/Reserving/SimulationWuethrich/Simulation.Machine.V1")
source(file="./Functions.V1.R")


#############################
### Generate the features ###
#############################

V <- 5000000                           # totally expected number of claims (over 12 accounting years)
LoB.dist <- c(0.35,0.20,0.15,0.30)    # categorical distribution for the allocation of the claims to the 4 lines of business
growth <- c(0,0,0,0)                  # growth parameters (per LoB) for the numbers of claims in the 12 accident years
seed1 <- 100                          # setting seed for simulation
features <- Feature.Generation(V = V, LoB.dist = LoB.dist, inflation = growth, seed1 = seed1)

str(features)

###############################################
### Simulate (and store) cash flow patterns ###
###############################################

npb <- nrow(features)                # blocks for parallel computing
seed1 <- 100                         # setting seed for simulation
std1 <- 0.85                         # standard deviation parameter for total claim size simulation
std2 <- 0.85                         # standard deviation parameter for recovery simulation
output <- Simulation.Machine(features = features, npb = npb, seed1 = seed1, std1 = std1, std2 = std2)

str(output)

##########################################################
### CL analysis on aggregated claims of simulated data ###
##########################################################

### We add artificial observations with 0 payments
### to ensure that we have at least one observation per accident year
add.obs <- output[rep(1,12),]
add.obs[,1] <- -1
add.obs[,4] <- 1994:2005
add.obs[,9:20] <- 0
output <- rbind(output,add.obs)

#alreadyReportedby2005<-output[output$AY+output$RepDel<=2005,]

### Cumulative cash flows
#cum_CF <- round((1000)^(-1)*ddply(output, .(AY), summarise, CF00=sum(Pay00),CF01=sum(Pay01),CF02=sum(Pay02),CF03=sum(Pay03),CF04=sum(Pay04),CF05=sum(Pay05),CF06=sum(Pay06),CF07=sum(Pay07),CF08=sum(Pay08),CF09=sum(Pay09),CF10=sum(Pay10),CF11=sum(Pay11))[,2:13])
cum_CF <- ddply(alreadyReportedby2005, .(AY), summarise, CF00=sum(Pay00),CF01=sum(Pay01),CF02=sum(Pay02),CF03=sum(Pay03),CF04=sum(Pay04),CF05=sum(Pay05),CF06=sum(Pay06),CF07=sum(Pay07),CF08=sum(Pay08),CF09=sum(Pay09),CF10=sum(Pay10),CF11=sum(Pay11))[,2:13]
for (j in 2:12){cum_CF[,j] <- cum_CF[,j-1] + cum_CF[,j]}
cum_CF

### Mack chain-ladder analysis
tri_dat <- array(NA, dim(cum_CF))
reserves <- data.frame(array(0, dim=c(12+1,3)))
reserves <- setNames(reserves, c("true Res.","CL Res.","MSEP^(1/2)"))
for (i in 0:11){
  for (j in 0:(11-i)){tri_dat[i+1,j+1] <- cum_CF[i+1,j+1]}
  reserves[i+1,1] <- cum_CF[i+1,12]-cum_CF[i+1,12-i]
}
reserves[13,1] <- sum(reserves[1:12,1])
tri_dat <- as.triangle(as.matrix(tri_dat))
dimnames(tri_dat)=list(origin=1:12, dev=1:12)
Mack <- MackChainLadder(tri_dat, est.sigma="Mack")
for (i in 0:11){reserves[i+1,2] <- round(Mack$FullTriangle[i+1,12]-Mack$FullTriangle[i+1,12-i])}
reserves[13,2] <- sum(reserves[1:12,2])
reserves[1:12,3] <- round(Mack$Mack.S.E[,12])
reserves[13,3] <- round(Mack$Total.Mack.S.E)
reserves                           # true reserves, chain-ladder reserves and square-rooted MSEP

#write triangle to file
write.table(tri_dat, "./triangle.csv", sep=",", row.names=FALSE)
write.table(reserves, "./truth.csv", sep=",", row.names=FALSE)

paidToDate<-diag(tri_dat[nrow(tri_dat):1,])[nrow(tri_dat):1]
ultimate<-paidToDate+reserves[1:nrow(tri_dat),2]
ratioToUltimate<-paidToDate/ultimate
   
#attach ratioToUltimate to data

ayears<-sort(unique(output$AY))
DFratioToUltimte<-data.frame(AY = ayears, ratioToUltimate = ratioToUltimate)
claimsData <- (merge(DFratioToUltimte, output, by = 'AY'))

#cumulative paid amounts
output$PayCum00<-output$Pay00
output$PayCum01<-output$PayCum00+output$Pay01
output$PayCum02<-output$PayCum01+output$Pay02
output$PayCum03<-output$PayCum02+output$Pay03
output$PayCum04<-output$PayCum03+output$Pay04
output$PayCum05<-output$PayCum04+output$Pay05
output$PayCum06<-output$PayCum05+output$Pay06
output$PayCum07<-output$PayCum06+output$Pay07
output$PayCum08<-output$PayCum07+output$Pay08
output$PayCum09<-output$PayCum08+output$Pay09
output$PayCum10<-output$PayCum09+output$Pay10
output$PayCum11<-output$PayCum10+output$Pay11

#discard add.obs rows
output2<-subset(output, ClNr > 0)
#remove certain columns 
output2$Pay00 <- NULL
output2$Pay01 <- NULL
output2$Pay02 <- NULL
output2$Pay03 <- NULL
output2$Pay04 <- NULL
output2$Pay05 <- NULL
output2$Pay06 <- NULL
output2$Pay07 <- NULL
output2$Pay08 <- NULL
output2$Pay09 <- NULL
output2$Pay10 <- NULL
output2$Pay11 <- NULL

output2$Open00 <- NULL
output2$Open01 <- NULL
output2$Open02 <- NULL
output2$Open03 <- NULL
output2$Open04 <- NULL
output2$Open05 <- NULL
output2$Open06 <- NULL
output2$Open07 <- NULL
output2$Open08 <- NULL
output2$Open09 <- NULL
output2$Open10 <- NULL
output2$Open11 <- NULL

write.table(output2, "./Simulated.Cashflows.csv", sep=",", row.names=FALSE)
#the data should be copied into the folder data\reservingAccident of the DecisionTrees.jl Julia package
