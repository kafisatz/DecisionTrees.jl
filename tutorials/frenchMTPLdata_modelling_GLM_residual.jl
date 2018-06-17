###########


t0=time_ns()
cd(joinpath(GLOBAL_julia_code_folder,"DecisionTrees.jl"))
#@info("You may want to run 'pkg> instantiate' when you first run this. Use ] to enter the package mode.")

@everywhere using Revise
@everywhere using CSV,DataFrames

@everywhere using DecisionTrees  

tela = (time_ns()-t0)/1e9
@show tela #precompilation can take considerable time under Julia 0.7alpha (up to 15 minutes in some cases (as ofJune 11, 2018)), this should improve soon though 

##############################
#Read the data
##############################
datafile=string("data\\freMTPL2\\freMTPL2.csv")
@assert isfile(datafile);
@time fullData=CSV.read(datafile,rows_for_type_detect=100000,allowmissing=:none,categorical=false);

#add AreaInteger as new variable, i.e. A:F -> 1:6 (as suggested in the paper of Wüthrich, Noll, Salzmann)
areasSorted=sort(unique(fullData[:Area]))
AreaInteger=map(x->findall((in)([x]),areasSorted)[1],fullData[:Area])
fullData[:AreaInteger]=AreaInteger

#correct for unreasonable observations
for i=1:size(fullData,1)
    fullData[:ClaimNb][i]=min(4,fullData[:ClaimNb][i])
    fullData[:Exposure][i]=min(1.0,fullData[:Exposure][i])
end

##############################
#Some Checks
##############################

#check poisson error versus errors reported in R
#this is a good check to ensure that the data is the same (train/test split, fitted GLM values are the same as in R)

glmFit=fullData[:GLMfitted]
errs=poissonError(fullData[:ClaimNb],fullData[:GLMfitted])

inSampleErrorR=31.26737970992647  #see R Code
outOfSampleErrorR=32.171234111851014 #see R Code
valBitArrayIdx=fullData[:trnTest].==0;
trnBitArrayIdx=fullData[:trnTest].==1;
trnAvgError=sum(errs[trnBitArrayIdx])/sum(trnBitArrayIdx)
valAvgError=sum(errs[valBitArrayIdx])/sum(valBitArrayIdx)
trnAvgError-inSampleErrorR/100 #should be zero
valAvgError-outOfSampleErrorR/100 #should be zero

originalTrnValIndex=deepcopy(fullData[:trnTest])    

##############################
#Prepare the data
##############################

#set independent variables
selected_explanatory_vars=["Area","AreaInteger","VehPower","VehAge","DrivAge","BonusMalus","VehBrand","VehGas","Density","Region"]

#Let us now model the RESIDUAL of the GLM.
#thus we model observed/glmfitted, if we can model that ratio successfully, our model will improve.
#dfprepped is an intermediary dataframe which is generally not needed
dtmtable,sett,dfprepped=prepare_dataframe_for_dtm!(fullData,keycol="IDpol",trnvalcol="trnTest",numcol="ClaimNb",denomcol="GLMfitted",weightcol="Exposure",independent_vars=selected_explanatory_vars);

##############################
#Define model SETTINGS
##############################

#let us consider a larger validation data (say 30%) set for the following
resample_trnvalidx!(dtmtable,.7)

updateSettingsMod!(sett,minw=-0.01,model_type="boosted_tree",niter=100,mf=0.01,subsampling_features_prop=.7,boolCalculatePoissonError=false)
sett.scorebandsstartingpoints=[1,50,100,150,200,300,400,500,600,700,850,900,950]

#resM is the resulting Model
resultingFiles,resM=dtm(dtmtable,sett)

