@assert VERSION>=v"0.7.0-beta" "Your Julia version does not satisfy the requirements for this package. You need to use Julia v0.7"
t0=time_ns()
#cd(string(ENV["HOMEPATH"],"\\Documents\\ASync\\home\\Code\\Julia\\DecisionTrees.jl"))

import Distributed
#addprocs(22)
#Distributed.@everywhere using Revise
Distributed.@everywhere import CSV
Distributed.@everywhere import DataFrames
Distributed.@everywhere import DataFrames: DataFrame

Distributed.@everywhere using DecisionTrees  #this may take some time

tela = (time_ns()-t0)/1e9
@show tela #precompilationcan take considerable time under Julia 0.7beta (up to 300 seconds in some cases (as ofJuly 3, 2018)) 


############################################################
#1. Model for the French MTPL Data
############################################################

##############################
#Read the data
##############################
datafile=joinpath("data","freMTPL2","freMTPL2.csv")
@assert isfile(datafile) "You need to unzip the file freMTPL2.zip in the folder data/freMTPL2"
@time fullData=CSV.read(datafile,rows_for_type_detect=100000,allowmissing=:none,categorical=false);

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
#note: we refrain from using the exposure as explanatory variable, as this is not meaningful in practical applications
#where one aims to predict the claim frequency of a given policy (without knowing how long that policy will remain within a certain portfolio)

dtmtable,sett,dfprepped=prepare_dataframe_for_dtm!(fullData,keycol="IDpol",trnvalcol="trnTest",numcol="ClaimNb",denomcol="Exposure",weightcol="Exposure",independent_vars=selected_explanatory_vars);


updateSettingsMod!(sett,minWeight=-0.03,model_type="build_tree",calculatePoissonError=true)

resultingFiles,resM=dtm(dtmtable,sett)