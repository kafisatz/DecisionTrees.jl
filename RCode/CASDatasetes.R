#use R to download the dataset to this folder

install.packages("CASdatasets", repos = "http://dutangc.free.fr/pub/RRepos/")
library(CASdatasets)
?CASdatasets

set.seed (100)
ll <- sample (c (1: nrow ( freMTPL2freq )), round (0.9* nrow ( freMTPL2freq )), replace = FALSE )
learn <- freMTPL2freq [ll ,]
test <- freMTPL2freq [-ll ,]
data(freMTPL2freq)

oneToN<-c(1: nrow ( freMTPL2freq ))

matched<-match(oneToN,ll)
matched[is.na(matched)]=0
matched[matched>0]=1

#matched is the train test identifier with entries 0/1


###############################################
#The following code is copied from code used for the Noll, Salzmann, Wüthrich paper
###############################################
#########  Data Analysis French MTPL2
#########  GLM, Regression Trees and Boosting
#########  Author: Mario Wuthrich
#########  Version May 6, 2018
#########  See also: https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3164764
###############################################

###############################################
#########  load packages and data
###############################################

require(MASS)
library(CASdatasets)
require(stats)
library(data.table)
library(plyr)
library(rpart)
library(rpart.plot)
library(Hmisc)



data(freMTPL2freq)
dat <- freMTPL2freq
dat$VehGas <- factor(dat$VehGas)      # consider VehGas as categorical
dat$n <- 1                            # intercept
dat$ClaimNb <- pmin(dat$ClaimNb, 4)   # correct for unreasonable observations (that might be data error)
dat$Exposure <- pmin(dat$Exposure, 1) # correct for unreasonable observations (that might be data error)
str(dat)


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
dat2$AreaGLM <- as.integer(dat2$Area)
dat2$VehPowerGLM <- as.factor(pmin(dat2$VehPower,9))
VehAgeGLM <- cbind(c(0:110), c(1, rep(2,10), rep(3,100)))
dat2$VehAgeGLM <- as.factor(VehAgeGLM[dat2$VehAge+1,2])
dat2[,"VehAgeGLM"] <-relevel(dat2[,"VehAgeGLM"], ref="2")
DrivAgeGLM <- cbind(c(18:100), c(rep(1,21-18), rep(2,26-21), rep(3,31-26), rep(4,41-31), rep(5,51-41), rep(6,71-51), rep(7,101-71)))
dat2$DrivAgeGLM <- as.factor(DrivAgeGLM[dat2$DrivAge-17,2])
dat2[,"DrivAgeGLM"] <-relevel(dat2[,"DrivAgeGLM"], ref="5")
dat2$BonusMalusGLM <- as.integer(pmin(dat2$BonusMalus, 150))
dat2$DensityGLM <- as.numeric(log(dat2$Density))
dat2[,"Region"] <-relevel(dat2[,"Region"], ref="R24")
dat2$AreaRT <- as.integer(dat2$Area)
dat2$VehGasRT <- as.integer(dat2$VehGas)
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

### Model GLM1
d.glm1 <- glm(ClaimNb ~ VehPowerGLM + VehAgeGLM + DrivAgeGLM + BonusMalusGLM
              + VehBrand + VehGas + DensityGLM + Region + AreaGLM, 
              data=learn, offset=log(Exposure), family=poisson())

summary(d.glm1)  
#anova(d.glm1)   # this calculation takes a bit of time                                          

learn$fitGLM <- fitted(d.glm1)
test$fitGLM <- predict(d.glm1, newdata=test, type="response")

glmfit_alldata<-c(as.vector(learn$fitGLM),as.vector(test$fitGLM))
ids<-c(learn$IDpol,test$IDpol)
srtperm<-order(ids)
ids[srtperm]
is.unsorted(ids[srtperm]) #must be FALSE

options(digits=22)
getOption("digits")
# in-sample and out-of-sample losses (in 10^(-2))
(insampleGLM <- 100*Poisson.Deviance(learn$fitGLM, learn$ClaimNb))
#31.26737970992647 
100*Poisson.Deviance(test$fitGLM, test$ClaimNb)
#32.171234111851014

#add train test index vector and GLM prediction to the data
freMTPL2freq$trnTest=matched

freMTPL2freq$GLMfitted=glmfit_alldata[srtperm]

write.csv(freMTPL2freq,"c:\\temp\\freMTPL2.csv",row.names=FALSE)
