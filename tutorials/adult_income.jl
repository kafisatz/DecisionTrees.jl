t0=time_ns()
cd(string(ENV["HOMEPATH"],"\\Documents\\ASync\\home\\Code\\Julia\\DecisionTrees.jl"))
@warn("You may need to run 'pkg> instantiate' when you first run this. Use ] to enter the package mode.")

using Revise
using CSV,DataFrames

using DecisionTrees  

tela = (time_ns()-t0)/1e9
@show tela #maybe 500-700 seconds on June 11, 2018

#elt=[Int,	Float64,	Float64,	Float64,	Int,	String,	String,	String,	Int,	String,	String,	String,	String,	String,	String,	String,	String,	String,	String,	Int,	Int,	Int,	Int,	String,	Int,	String,	String,	Int,	String,	String,	Int,	String,	String,	Int,	String,	Int,	String,	String,	Int,	Int,	String,	Int,	Int,	String,	Int,	Int,	String,	Int,	String,	String,	String,	Int,	String,	String,	String,	String,	Int,	Int,	Int,	String,	String,	String,	String,	String,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Float64,	Int,	Float64,	Float64,	Float64,	Float64,	Int,	Int,	Int,	Int,	Int,	Int]
header=["age","workclass","fnlwgt","education","education-num","marital-status","occupation","relationship","race","sex","capital-gain","capital-loss","hours-per-week","native-country","target"]
#@time df_tmp=readtable(string("data\\data1mini.csv"),eltypes=elt); 
@time dftrain=CSV.read(string("data\\adult\\adult.data"),header=header,rows_for_type_detect=10000,allowmissing=:none,categorical=false,rows=32561);
@time dftest=CSV.read(string("data\\adult\\adult.test"),header=header,rows_for_type_detect=10000,allowmissing=:none,categorical=false,rows=16281);

dftest[:trainOrTest]="Test"
dftrain[:trainOrTest]="Train"

#data prep is done for the aggregated data set in order to ensure that the factors are the same
dfTot=vcat(dftrain,dftest)

#note:the target variable in test and train is not formatted in the same manner! One of the files has a "." at the end
for i=1:size(dfTot,1)
    v=dfTot[:target][i]
    if v==" >50K." || v==" >50K"
        dfTot[:target][i]=" >50K"
    else
        dfTot[:target][i]=" <=50K"
    end
end

#create numerical response variable
dfTot[:response]=ifelse.(dfTot[:target].==" >50K",1,0)
dfTot[:trnVal]=ifelse.(dfTot[:trainOrTest].=="Test",0,1)

selected_explanatory_vars=header[1:end-1]
#view coltypes
for x in selected_explanatory_vars
    @show x,eltype(dfTot[Symbol(x)])
    println("first ten values")
    print(dfTot[1:10,Symbol(x)])
    println("")
end

dtmtableFull,sett=prepare_dataframe_for_dtm!(dfTot,trnvalcol="trnVal",directory="C:\\temp\\",numcol="response",independent_vars=selected_explanatory_vars);

#split the dtmtableFull into two separate tables
dtmtableTrain=deepcopy(dtmtableFull);
dtmtableTrain.key=dtmtableTrain.key[1:size(dftrain,1)];
dtmtableTrain.numerator=dtmtableTrain.numerator[1:size(dftrain,1)];
dtmtableTrain.denominator=dtmtableTrain.denominator[1:size(dftrain,1)];
dtmtableTrain.weight=dtmtableTrain.weight[1:size(dftrain,1)];
dtmtableTrain.features=dtmtableTrain.features[1:size(dftrain,1),:];
#temporarily reassign tnr and val indices. These will be overwritten later on
dtmtableTrain.validx=dtmtableTrain.trnidx[end-1:end];
dtmtableTrain.trnidx=dtmtableTrain.trnidx[1:end-2];
    
dtmtableTest=deepcopy(dtmtableFull);
dtmtableTest.key=dtmtableTest.key[size(dftrain,1)+1:end];
dtmtableTest.numerator=dtmtableTest.numerator[size(dftrain,1)+1:end];
dtmtableTest.denominator=dtmtableTest.denominator[size(dftrain,1)+1:end];
dtmtableTest.weight=dtmtableTest.weight[size(dftrain,1)+1:end];
dtmtableTest.features=dtmtableTest.features[size(dftrain,1)+1:end,:];
#temporarily reassign tnr and val indices. These will be overwritten later on
dtmtableTest.validx=convert(Vector{Int},1:(length(dtmtableTest.weight)-2))
dtmtableTest.trnidx=convert(Vector{Int},(length(dtmtableTest.weight)-2):length(dtmtableTest.weight))
    
#define trn and val data sets
resample_trnvalidx!(dtmtableTrain,.7)

#run model
dtm(dtmtableTrain,sett)

#grid search

#evaluate confusion matrix
0

favorite_settings=deepcopy(sett)
favorite_settings.minw=250
dt=deepcopy(dtmtableTrain)
#derive model from all data
resample_trnvalidx!(dtmtableTrain,.999)
notUsed,resModel=dtm(dt,sett)

predictionsOnHoldOut,leafnrs=predict(resModel,dtmtableTest.features)
preds=ifelse.(predictionsOnHoldOut.<.5,1.0,0.0)
truth=dtmtableTest.numerator

err,accuracy,mat=confusmatBinary(1.+convert(Vector{Int},truth),1.+convert(Vector{Int},preds))
#tp 2281
#tn 621
# R[2,1] =FP = 1565 preds.==0 and truth.==1
#
# dtm(dtmtableTrain,settingsVector)