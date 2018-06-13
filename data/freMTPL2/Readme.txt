#use R to download the dataset to this folder

install.packages("CASdatasets", repos = "http://dutangc.free.fr/pub/RRepos/")
library(CASdatasets)
?CASdatasets
data(freMTPL2freq)
write.csv(freMTPL2freq,"DecisionTrees.jl\\data\\freMTPL2.csv"