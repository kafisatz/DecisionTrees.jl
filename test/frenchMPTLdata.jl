t0=time_ns()
cd(string(ENV["HOMEPATH"],"\\Documents\\ASync\\home\\Code\\Julia\\DecisionTrees.jl"))
@warn("You may need to run 'pkg> instantiate' when you first run this. Use ] to enter the package mode.")

using Revise
using CSV,DataFrames

using DecisionTrees  

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

#evaluate confusion matrix
0

#=
#test grid search

sett2=deepcopy(sett)
sett3=deepcopy(sett)
sett4=deepcopy(sett)
sett2.minw=sett.minw/2
sett3.minw=sett.minw/6
sett4.niter=sett.niter*2

setts=vcat(sett,sett2,sett3,sett4)
=#

#=
favorite_settings=deepcopy(sett)
favorite_settings.minw=250
dt=deepcopy(dtmtableTrain)
#derive model from all data
resample_trnvalidx!(dtmtableTrain,.999)
notUsed,resModel=dtm(dt,sett)

predictionsOnHoldOut,leafnrs=predict(resModel,dtmtableTest.features)
preds=ifelse.(predictionsOnHoldOut.<.5,1.0,0.0)
truth=dtmtableTest.numerator

err,accuracy,mat=confusmatBinary(1.+convert(Vector{Int64},truth),1.+convert(Vector{Int64},preds))
=#