import CSV
import DataFrames
import DataFrames: DataFrame

using DecisionTrees
using Profile
using Pkg

pkgdir=Pkg.dir("DecisionTrees") #how long will this remain supported?
testdir=joinpath(pkgdir,"test")
datadir=joinpath(pkgdir,"data")

#Read the data
datafile=joinpath(datadir,"freMTPL2\\freMTPL2.csv")
@assert isfile(datafile);
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

#Prepare the data
##############################

dtmtable,sett,dfprepped=prepare_dataframe_for_dtm!(fullData,keycol="IDpol",trnvalcol="trnTest",numcol="ClaimNb",denomcol="Exposure",weightcol="Exposure",independent_vars=selected_explanatory_vars);

##############################
#Define model SETTINGS

updateSettingsMod!(sett,minWeight=-0.03,model_type="boosted_tree",iterations=2,learningRate=0.1,subsampling_features_prop=.7,calculatePoissonError=true)

##############################
#Run single Boosting Model
resultingFiles,resM=dtm(dtmtable,sett)

#measure allocation
updateSettingsMod!(sett,
iterations=40,
model_type="boosted_tree",
nScores="1000",
writeStatistics="false",
writeSasCode="false",
writeIterationMatrix="false",
writeResult="false",
writeCsharpCode="false",
writeVbaCode="false",
saveJLDFile="false",
saveResultAsJLDFile="false",

showProgressBar_time="false",
prroduceEstAndLeafMatrices="false",
write_dot_graph="false",

calculatePoissonError="true",
performanceMeasure="Average Poisson Error Val"
)


Profile.clear_malloc_data()
strs,resm2=dtm(dtmtable,sett)
@info("allocation run done.")
quit();
