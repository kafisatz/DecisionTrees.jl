t0=time_ns()
cd(string(ENV["HOMEPATH"],"\\Documents\\ASync\\home\\Code\\Julia\\DecisionTrees.jl"))
@warn("You may need to run 'pkg> instantiate' when you first run this. Use ] to enter the package mode.")

@everywhere using Revise
@everywhere using CSV,DataFrames

@everywhere using DecisionTrees  

tela = (time_ns()-t0)/1e9
@show tela #maybe 500-700 seconds on June 11, 2018

#elt=[Int64,	Float64,	Float64,	Float64,	Int64,	String,	String,	String,	Int64,	String,	String,	String,	String,	String,	String,	String,	String,	String,	String,	Int64,	Int64,	Int64,	Int64,	String,	Int64,	String,	String,	Int64,	String,	String,	Int64,	String,	String,	Int64,	String,	Int64,	String,	String,	Int64,	Int64,	String,	Int64,	Int64,	String,	Int64,	Int64,	String,	Int64,	String,	String,	String,	Int64,	String,	String,	String,	String,	Int64,	Int64,	Int64,	String,	String,	String,	String,	String,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Float64,	Int64,	Float64,	Float64,	Float64,	Float64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64]
#header=["age","workclass","fnlwgt","education","education-num","marital-status","occupation","relationship","race","sex","capital-gain","capital-loss","hours-per-week","native-country","target"]
#@time df_tmp=readtable(string("data\\data1mini.csv"),eltypes=elt); 
@time fullData=CSV.read(string("data\\freMTPL2\\freMTPL2.csv"),rows_for_type_detect=100000,allowmissing=:none,categorical=false);

#add AreaInteger as new variable, being A:F -> 1:6 (as suggested in the paper of Wüthrich)
areasSorted=sort(unique(fullData[:Area]))
AreaInteger=map(x->findall((in)([x]),areasSorted)[1],fullData[:Area])
fullData[:AreaInteger]=AreaInteger

selected_explanatory_vars=["Area","AreaInteger","VehPower","VehAge","DrivAge","BonusMalus","VehBrand","VehGas","Density","Region"]
#view coltypes
for x in selected_explanatory_vars
    println(x)
    println(eltype(fullData[Symbol(x)]))
    println("first ten values")
    print(fullData[1:10,Symbol(x)])
    println("")
end

dtmtable,sett=prepare_dataframe_for_dtm!(fullData,trnvalcol="trnTest",directory="C:\\temp\\",numcol="ClaimNb",denomcol="Exposure",weightcol="Exposure",independent_vars=selected_explanatory_vars);

originalTrnValIndex=deepcopy(fullData[:trnTest])    

#redefine trn and val data sets
resample_trnvalidx!(dtmtable,.7)

#define model settings

sett.minw=-.03
sett.model_type="boosted_tree"
sett.niter=10
sett.mf=0.08
sett.subsampling_features_prop=.7
sett.boolCalculatePoissonError=true

#run model
unused,resM=dtm(dtmtable,sett)
fitted=resM.meanobserved.*resM.rawrelativities

errs=poissonError(dtmtable.numerator,dtmtable.weight,fitted);
sum(errs[dtmtable.trnidx])/length(dtmtable.trnidx)
sum(errs[dtmtable.validx])/length(dtmtable.validx)


#grid search

#sett core settings
updateSettingsMod!(sett,
niter=5,
model_type="boosted_tree",
max_splitting_points_num=250	,
nscores="1000",
write_statistics="true",
write_sas_code="false",
write_iteration_matrix="false",
write_result="false",
write_csharp_code="false",
write_vba_code="false",
boolSaveJLDFile="false",
boolSaveResultAsJLDFile="false",

showProgressBar_time="true",
boolProduceEstAndLeafMatrices="false",
write_dot_graph="false",

boolCalculatePoissonError="true",
performanceMeasure="Average Poisson Error Val"
)

#createGridSearchSettings
settV=createGridSearchSettings(sett,    
    minw=[-0.005,-0.01,-0.02,-0.03,-0.05]
)
    #,mf=[0.01,0.005,0.02,0.04],
    #subsampling_features_prop=[1.0,.75,.5,.25],
    #smoothEstimates=["0","1"]
#
#    randomw=0.0
#	subsampling_prop)

#dtm(dtmtable,settV,fn="C:\\temp\\1\\dtm.CSV")
@show 0

warn("todo:check best tree with poisson error too!")
    #crit::SplittingCriterion # fn version of 4