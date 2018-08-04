#git clone --depth=50 --branch=master https://github.com/kafisatz/DecisionTrees.jl.git 
#cd kafisatz/DecisionTrees.jl

@assert VERSION>=v"0.7.0-rc2.0"

import Distributed

using Pkg
pkg"activate DecisionTrees.jl"
Distributed.@everywhere using DecisionTrees

dr="C:\\Users\\bernhard.konig\\Documents\\ASync\\home\\Code\\Julia\\DecisionTrees.jl"
isdir(dr)&&cd(dr)
pwd()

Distributed.@everywhere import CSV
Distributed.@everywhere import DataFrames
Distributed.@everywhere import DataFrames: DataFrame

##############################
#Read the data
##############################
@warn("""This assumes pwd() == Pkg.dir("DecisionTrees.jl")""") 
datafile=joinpath("data","freMTPL2","freMTPL2.csv")
@assert isfile(datafile) "You need to unzip the file freMTPL2.zip in the folder data/freMTPL2"
@time fullData=CSV.read(datafile,rows_for_type_detect=100000,allowmissing=:none,categorical=false);

#some data manipulation
    areasSorted=sort(unique(fullData[:Area]))
    AreaInteger=map(x->findall((in)([x]),areasSorted)[1],fullData[:Area])
    fullData[:AreaInteger]=AreaInteger
    #correct for unreasonable observations
    for i=1:size(fullData,1)
        fullData[:ClaimNb][i]=min(4,fullData[:ClaimNb][i])
        fullData[:Exposure][i]=min(1.0,fullData[:Exposure][i])
    end
    
#set independent variables
selected_explanatory_vars=["Area","AreaInteger","VehPower","VehAge","DrivAge","BonusMalus","VehBrand","VehGas","Density","Region"]

##############################
#Prepare the data
##############################
dtmtable,sett,dfprepped=prepare_dataframe_for_dtm!(fullData,keycol="IDpol",trnvalcol="trnTest",numcol="ClaimNb",denomcol="Exposure",weightcol="Exposure",independent_vars=selected_explanatory_vars);

#run a few models

#this runs through as expected
    updateSettings!(sett,minWeight=-0.03,model_type="build_tree",calculatePoissonError=true)
    resultingFiles,resM=dtm(dtmtable,sett)

#this runs through as expected
    updateSettings!(sett,crit="sse",model_type="build_tree")
    resultingFiles,resM=dtm(dtmtable,sett)

#this runs through as expected
    updateSettings!(sett,crit="sse",model_type="build_tree")
    resultingFiles,resM=dtm(dtmtable,sett)

#this one fails!
    updateSettings!(sett,crit="sse",model_type="boosted_tree")
    resultingFiles,resM=dtm(dtmtable,sett)