#copied from here
#https://www.hackerearth.com/practice/machine-learning/machine-learning-algorithms/beginners-tutorial-on-xgboost-parameter-tuning-r/tutorial/

#set working directory
path <- "C:\\Users\\bernhard.konig\\Documents\\ASync\\home\\Code\\Julia\\DecisionTrees.jl\\data\\adult\\"
setwd(path)

#load libraries
library(data.table)
library(mlr)
library(xgboost)

#set variable names
setcol <- c("age","workclass","fnlwgt","education","education-num","marital-status","occupation","relationship","race","sex","capital-gain","capital-loss","hours-per-week","native-country","target")

#load data
train <- read.table("adult.data", header = F, sep = ",", col.names = setcol, na.strings = c(" ?"), stringsAsFactors = F)
#test <- read.table("adult.test",header = F,sep = ",",col.names = setcol,skip = 1, na.strings = c(" ?"),stringsAsFactors = F) #if you want to skip row 1
test <- read.table("adult.test",header = F,sep = ",",col.names = setcol,na.strings = c(" ?"),stringsAsFactors = F)
#"74","native.countryHoland-Netherlands"
#BK: train has an observation with native.coutnry="Holand-Netherlands" which does not occur in test
#we remove this to avoid an error when applying the data to test
idx<-train$native.country==" Holand-Netherlands"
idx[is.na(idx)] <- FALSE 
sum(idx)
#remove single observation
train<-train[!idx,]

#convert data frame to data table
setDT(train) 
setDT(test)

#check missing values 
table(is.na(train))
sapply(train, function(x) sum(is.na(x))/length(x))*100

table(is.na(test))
sapply(test, function(x) sum(is.na(x))/length(x))*100

#quick data cleaning
#remove extra character from target variable
library(stringr)
test [,target := substr(target,start = 1,stop = nchar(target)-1)]

#remove leading whitespaces
char_col <- colnames(train)[ sapply (test,is.character)]
for(i in char_col) set(train,j=i,value = str_trim(train[[i]],side = "left"))

for(i in char_col) set(test,j=i,value = str_trim(test[[i]],side = "left"))

#set all missing value as "Missing" 
train[is.na(train)] <- "Missing" 
test[is.na(test)] <- "Missing"

#using one hot encoding 
labels <- train$target 
labels2 <- labels==unique(labels)[1] #convert from 50K to 1 and 0
ts_label <- test$target
ts_label2 <- ts_label==unique(labels)[1]
new_tr <- model.matrix(~.+0,data = train[,-c("target"),with=F]) 
new_ts <- model.matrix(~.+0,data = test[,-c("target"),with=F])

#convert factor to numeric 
labels <- as.numeric(labels2)
ts_label <- as.numeric(ts_label2)

#preparing matrix 
dtrain <- xgb.DMatrix(data = new_tr,label = labels) 
dtest <- xgb.DMatrix(data = new_ts,label=ts_label)

#default parameters
params <- list(booster = "gbtree", objective = "binary:logistic", eta=0.3, gamma=0, max_depth=6, min_child_weight=1, subsample=1, colsample_bytree=1)

xgbcv <- xgb.cv( params = params, data = dtrain, nrounds = 100, nfold = 5, showsd = T, stratified = T, print.every.n = 10, early.stop.round = 20, maximize = F)
##best iteration = 39

xgbcv
#test error is at 0.1264702 for iteration 39

#first default - model training
xgb1 <- xgb.train (params = params, data = dtrain, nrounds = 79, watchlist = list(val=dtest,train=dtrain), print.every.n = 10, early.stop.round = 10, maximize = F , eval_metric = "error")

#check that names match
  n1<-dimnames(dtrain)[[2]]
  n2<-dimnames(dtest)[[2]]
  cpm<-n1==n2
  sum(cpm)==length(n1)

#model prediction
#xgbpred <- predict (xgb1,dtrain)
xgbpred <- predict (xgb1,dtest)
xgbpred <- ifelse (xgbpred > 0.5,1,0)

#confusion matrix
library(caret)
library(e1071)
confusionMatrix (as.factor(xgbpred), as.factor(ts_label))
#Accuracy - 0.8733 

#view variable importance plot
mat <- xgb.importance (feature_names = colnames(new_tr),model = xgb1)
xgb.plot.importance (importance_matrix = mat[1:20])
