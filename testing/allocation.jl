
using DecisionTrees
using CSV
using DataFrames
using Profile

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

updateSettingsMod!(sett,minw=-0.03,model_type="boosted_tree",niter=2,mf=0.1,subsampling_features_prop=.7,boolCalculatePoissonError=true)

##############################
#Run single Boosting Model
resultingFiles,resM=dtm(dtmtable,sett)

#measure allocation
updateSettingsMod!(sett,
niter=40,
model_type="boosted_tree",
nscores="1000",
write_statistics="false",
write_sas_code="false",
write_iteration_matrix="false",
write_result="false",
write_csharp_code="false",
write_vba_code="false",
boolSaveJLDFile="false",
boolSaveResultAsJLDFile="false",

showProgressBar_time="false",
boolProduceEstAndLeafMatrices="false",
write_dot_graph="false",

boolCalculatePoissonError="true",
performanceMeasure="Average Poisson Error Val"
)


Profile.clear_malloc_data()
strs,resm2=dtm(dtmtable,sett)
@info("allocation run done.")
quit();
