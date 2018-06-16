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
myVectors <- sapply(1:length(myList), function(i) any(myList[[i]] == myValue))

matched<-match(oneToN,ll)
matched[is.na(matched)]=0
matched[matched>0]=1

freMTPL2freq$trnTest=matched

write.csv(freMTPL2freq,"c:\\temp\\freMTPL2.csv",row.names=FALSE)
