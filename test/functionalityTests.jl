#smoketest_tree.jl

@testset "Functionality Tests" begin

import Random: rand,randstring

##################################################
#Read the data
##################################################

datafile=joinpath(datadir,"freMTPL2\\freMTPL2.csv")
@test isfile(datafile);
@time fullData=CSV.read(datafile,rows_for_type_detect=100000,allowmissing=:none,categorical=false);

#add AreaInteger as new variable, i.e. A:F -> 1:6
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
updateSettingsMod!(sett,minw=-0.03,model_type="build_tree",boolCalculatePoissonError=true)

##############################
#Run single tree model
##############################

resultingFiles,resM=dtm(dtmtable,sett)
@warn("todo: add tests on lift, poisson error, etc.")

##################################################
#run boosting
##################################################

sett.niter=15
sett.mf=0.1
sett.model_type="boosted_tree"
strs,resm2=dtm(dtmtable,sett)
@test typeof(resm2.modelstats)==DataFrame
@test 1==1

end #testset
