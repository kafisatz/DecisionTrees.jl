###############################################
#########  Data Analysis French MTPL2
#########  GLM, Regression Trees and Boosting
#########  Author: Bernhard König
#########  This code is based on the work of Noll, Salzmann and Wuethrich
#########  See also: https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3164764
###############################################

###############################################
#########  load packages and data
###############################################
 
require(MASS)
install.packages("CASdatasets", repos = "http://dutangc.free.fr/pub/RRepos/")
library(CASdatasets) 
require(stats)
library(data.table)
library(plyr)
library(rpart)
library(rpart.plot)
library(Hmisc)
library(compare)

#Read datafile with attached Scores (produced by Boosting Model in Julia)
fi<-'C:\\Users\\bernhard.konig\\Documents\\async\\home\\code\\julia\\DecisionTrees.jl\\data\\freMTPL2\\freMTPL2_withJuliaBoostingScores1000.csv'
dat<-fread(fi)
data(freMTPL2freq)

#amend the original data
datOrig<-freMTPL2freq
datOrig$ClaimNb <- pmin(datOrig$ClaimNb, 4)   # correct for unreasonable observations (that might be data error)
datOrig$Exposure <- pmin(datOrig$Exposure, 1) # correct for unreasonable observations (that might be data error)

#amend dat
dat$VehGas <- factor(dat$VehGas)      # consider VehGas as categorical
dat$n <- 1                            # intercept
dat$Area<-as.factor(dat$Area)

#compare data
comparison<-compare(datOrig,dat,allowAll=TRUE)
comparison$result
summary(comparison)
comparison$detailedResult


###############################################
#########  Poisson deviance statistics
###############################################

Poisson.Deviance <- function(pred, obs){
     2*(sum(pred)-sum(obs)+sum(log((obs/pred)^(obs))))/length(pred)
               }

###############################################
#########  feature pre-processing for GLM
###############################################

dat2 <- dat
dat2$ClaimNb<-as.numeric(dat2$ClaimNb)
dat2$VehBrand <- as.factor(dat2$VehBrand)
dat2$AreaGLM <- as.integer(dat2$Area)
dat2$VehPowerGLM <- as.factor(pmin(dat2$VehPower,9))
VehAgeGLM <- cbind(c(0:110), c(1, rep(2,10), rep(3,100)))
dat2$VehAgeGLM <- as.factor(VehAgeGLM[dat2$VehAge+1,2])
relevel(as.factor(dat2$VehAgeGLM),ref="2")
#dat2[,"VehAgeGLM"] <-relevel(dat2[,"VehAgeGLM"], ref="2")
DrivAgeGLM <- cbind(c(18:100), c(rep(1,21-18), rep(2,26-21), rep(3,31-26), rep(4,41-31), rep(5,51-41), rep(6,71-51), rep(7,101-71)))
dat2$DrivAgeGLM <- as.factor(DrivAgeGLM[dat2$DrivAge-17,2])
relevel(dat2$DrivAgeGLM, ref="5")
#dat2[,"DrivAgeGLM"] <-relevel(dat2[,"DrivAgeGLM"], ref="5")
dat2$BonusMalusGLM <- as.integer(pmin(dat2$BonusMalus, 150))
dat2$DensityGLM <- as.numeric(log(dat2$Density))
dat2$Region<-as.factor(dat2$Region)
relevel(dat2$Region, ref="R24")
#dat2[,"Region"] <-relevel(dat2[,"Region"], ref="R24")
dat2$AreaRT <- as.integer(dat2$Area)
dat2$VehGasRT <- as.integer(dat2$VehGas)
#let us construct N score bands
nscorebands<-50
dat2$scorebands<-as.factor(ceil(dat2$scores/(1000/nscorebands)))
str(dat2)

###############################################
#########  choosing learning and test sample
###############################################

set.seed(100)
ll <- sample(c(1:nrow(dat2)), round(0.9*nrow(dat2)), replace = FALSE)
learn <- dat2[ll,]
test <- dat2[-ll,]
(n_l <- nrow(learn))
(n_t <- nrow(test))

###############################################
#########  GLM analysis
###############################################

### Model GLM1, wihtout the Score
d.glm1 <- glm(ClaimNb ~ VehPowerGLM + VehAgeGLM + DrivAgeGLM + BonusMalusGLM
                        + VehBrand + VehGas + DensityGLM + Region + AreaGLM, 
                        data=learn, offset=log(Exposure), family=poisson())
                   
summary(d.glm1)  
#anova(d.glm1)   # this calculation takes a bit of time                                          

learn$fitGLM <- fitted(d.glm1)
test$fitGLM <- predict(d.glm1, newdata=test, type="response")

options(digits=22)
getOption("digits")

# in-sample and out-of-sample losses (in 10^(-2))
(insampleGLM <- 100*Poisson.Deviance(learn$fitGLM, learn$ClaimNb)) #31.26738
100*Poisson.Deviance(test$fitGLM, test$ClaimNb) #32.17123


#Consider the GLM with the score variable

### Model GLM1, wihtout the Score
d.glm2 <- glm(ClaimNb ~ VehPowerGLM + VehAgeGLM + DrivAgeGLM + BonusMalusGLM
              + VehBrand + VehGas + DensityGLM + Region + AreaGLM +scorebands, 
              data=learn, offset=log(Exposure), family=poisson())

summary(d.glm2)  
#anova(d.glm2)   # this calculation takes a bit of time                                          

#we can see two things
#1) that the score is highly significant (actually it is the most significant factor)
#2) in return the significance of the other variables decrease
#both of the above are expected due to the way the score was defined

learn$fitGLM2 <- fitted(d.glm2)
test$fitGLM2 <- predict(d.glm2, newdata=test, type="response")

# in-sample and out-of-sample losses (in 10^(-2))
(insampleGLM <- 100*Poisson.Deviance(learn$fitGLM2, learn$ClaimNb)) #
100*Poisson.Deviance(test$fitGLM2, test$ClaimNb) 

#VALIDATION ERROR
#when including scores as log linear factor
  #32.005868191510004  (train: 30.13709594018319)
#when including 10 discrete scorebands (as unordered factor)
  #31.835941361638032 (train: 30.044557985296333)
#when including 20 discrete scorebands (as unordered factor)
  #31.691920201655044 (train: 29.802325706494621)
#when including 50 discrete scorebands (as unordered factor)
  #31.586166499590558 (train: 29.66890227788042 )

#todo add single tree