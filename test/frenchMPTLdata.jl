############################################################
#Boosting Model for the French MTPL Data
############################################################

##############################
#Load packages
##############################
t0=time_ns()
cd(string(ENV["HOMEPATH"],"\\Documents\\ASync\\home\\Code\\Julia\\DecisionTrees.jl"))
#@info("You may want to run 'pkg> instantiate' when you first run this. Use ] to enter the package mode.")

@everywhere using Revise
@everywhere using CSV,DataFrames

@everywhere using DecisionTrees  

tela = (time_ns()-t0)/1e9
@show tela #precompilationcan take considerable time under Julia 0.7alpha (up to 15 minutes in some cases (as ofJune 11, 2018)), this should improve soon though 


##############################
#Read the data
##############################
datafile=string("data\\freMTPL2\\freMTPL2.csv")
@assert isfile(datafile);
@time fullData=CSV.read(datafile,rows_for_type_detect=100000,allowmissing=:none,categorical=false);

#the data is described in the R package CASdatasets
#also the Noll, Salzmann, Wüthrich Paper has some descriptive graphs of the data.

#add AreaInteger as new variable, i.e. A:F -> 1:6 (as suggested in the paper of Wüthrich, Noll, Salzmann)
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

#check type of each column
for x in selected_explanatory_vars
    println(x)
    println(eltype(fullData[Symbol(x)]))
    println("first ten values")
    print(fullData[1:10,Symbol(x)])
    println("")
end

#consider unique elements for each variable
for x in 1:size(fullData,2)
    println(x)
    @show size(unique(fullData[:,x]))
end

##############################
#Prepare the data
##############################

#dtmtable::DTMTable is a custom data format for DecisionTrees
#sett are the model settings
#dfprepped is an intermediary dataframe which is generally not needed
dtmtable,sett,dfprepped=prepare_dataframe_for_dtm!(fullData,trnvalcol="trnTest",numcol="ClaimNb",denomcol="Exposure",weightcol="Exposure",independent_vars=selected_explanatory_vars);
#consider 
fieldnames(dtmtable)
dtmtable.key #an identifier (String type), ideally it is a unique identifier for each row
dtmtable.numerator #a vector
dtmtable.denominator #a vector
dtmtable.weights #a vector
#each of the above have the same length

dtmtable.trnidx #an index vector for the training data. This is an integer vector
#therefore we generally have
length(dtmtable.trnidx)<length(dtmtable.numerator)
#dtmtable.validx is the corresponding validation data
#we trnidx and validx should generally be a partition of 1:n (where n==length(dtmtable.weights))

dtmtable.features #explanatory variables
#note that the data is generally stored in a compressed format
#for numerical variables, the algorithm considers sett.max_splitting_points_num  (default 250) splitting points which are uniformly spaced
#you can increase this number and re-prepare the data if you want
#=
    sett.max_splitting_points_num=100
    dtmtable=prep_data_from_df(df_prepped,sett,fn)
=#
#we realize that it may be sutiable to define the splitting in a different manner (than uniformly spaced).
#this feature might be added at a later time. You can consider the function add_coded_numdata!(...) to see how splitting points are chosen

originalTrnValIndex=deepcopy(fullData[:trnTest])    

#Redefine trn and val data sets
#if you prefer train on X% of the data, you can use this function to adjust the training set. By default it is sampled randomly
#resample_trnvalidx!(dtmtable,.85)

##############################
#Define model SETTINGS
##############################

#the minimum weight of each individual tree shall be 3% of the exposure per leaf
#You can set minw to any positive value, e.g. sett.minw=20000
#if minw is set to a value in [-1,0] its absolute value is interpreted as a percentage
#thus with sett.minw=-0.03 we define the min weight of each leaf as 3% of the training data

#model_type: currently only "build_tree" and "boosted_tree" are supported
#we might add the implementation of bagging at a later stage

#niter is the number of trees

#mf is the moderation factor which is also know as shrinkage or learning rate

#subsampling_features_prop must be between 0 and 1 and defines the subsampling probability of the predictors for each tree 
#i.e. with subsampling_features_prop=0.5 each tree will only use half of the predictors

#boolCalculatePoissonError, is a boolean which we need to enable here in order for output to show the poisson error of the estimates
updateSettingsMod!(sett,minw=-0.03,model_type="boosted_tree",niter=40,mf=0.025,subsampling_features_prop=.7,boolCalculatePoissonError=true)

##############################
#Run single Boosting Model
##############################

#resM is the resulting Model
resultingFiles,resM=dtm(dtmtable,sett)
#consider the resulting Excel file for a descritipon of the model 
#it should be located here:
@show resultingFiles[2]

#if you prefer to analyze the data in Julia instead of excel, consider
sheetnames=map(x->x.name,resM.exceldata.sheets)
#the following is the ModelStatistics sheets as a Dataframe 
resM.exceldata.sheets[3].data
#this is the lift in the training data for each iteration of the boosting model
resM.exceldata.sheets[3].data[:x6]

############################################################
#Investigate the predictions (fitted values)
############################################################

#DecisionTrees Boosting Models provide a number of different estimates
#1. Fitted values based on the "raw relativities"
#this is simply the combination of all decision trees

#2. Smoothed fitted values per Score. 
#the user can specify the number of scores for the output (see sett.nscores which defaults to 1000)
#DecisionTrees runs a moving average over the observed data to smooth the noise

#3. Unsmoothed fitted values per Score. 
#As the score is merely an aggregation of the resulting 'relativities' (i.e. the spread of the fitted frequencies by the ensbemble)
#one can calculate the observed frequency for each score. Here, no smoothing is performed


#If you want, you can consider the function below to get the fitted value
# rawFittedValues,smoothedFittedValuesPerScore,unsmoothedFittedValuesPerScore=getFittedValues(reM)
#in the following lines we derive the fitted values from the resM struct 

#We will show the different fitted values in the following

#1. Fitted values based on the "raw relativities"
#for the fitted values, consider:
fitted=resM.meanobserved.*resM.rawrelativities
#resM.meanobserved is the mean observed value of the training data
#rawrelativities is the estimated relative deviation from the mean for each observation
#Note:as for the input data the fitted vector corresponds to all the data.
#dtmtable.validx determines which observations correspond to the validation data

#Let us calculate the Poisson Error for these estimates
errs=poissonError(dtmtable.numerator,dtmtable.weight,fitted);
trnAvgError=sum(errs[dtmtable.trnidx])/length(dtmtable.trnidx) #should be around 0.30939804739941096
valAvgError=sum(errs[dtmtable.validx])/length(dtmtable.validx) #should be around 0.31950570163159764

#2. Smoothed fitted values per Score. 
scores=resM.scores
fittedThroughScores=zeros(Float64,length(scores))
for i=1:length(scores)
    @inbounds sc=scores[i]
    @inbounds fittedThroughScores[i]=resM.ScoreToSmoothedEstimate[sc]
end

fitted=fittedThroughScores
#notably for these fitted values we have less unique values
@show unique(fitted)
#this should correspond to sett.nscores
errs=poissonError(dtmtable.numerator,dtmtable.weight,fitted);
trnAvgError=sum(errs[dtmtable.trnidx])/length(dtmtable.trnidx) #0.3087311765403276
valAvgError=sum(errs[dtmtable.validx])/length(dtmtable.validx) #0.31904324177240373

#3. (raw) Observed frequency per Score
#As the score is merely an aggregation of the resulting 'relativities' (i.e. the spread of the fitted frequencies by the ensbemble)
#one can calculate the observed frequency for each score. Here, no smoothing is performed
scores=resM.scores
rawFittedThroughScores=zeros(Float64,length(scores))
for i=1:length(scores)
    @inbounds sc=scores[i]
    @inbounds rawFittedThroughScores[i]=resM.rawObservedRatioPerScore[sc]
end

#Again, let us consider the Poisson Error
fitted=rawFittedThroughScores
errs=poissonError(dtmtable.numerator,dtmtable.weight,fitted);
trnAvgError=sum(errs[dtmtable.trnidx])/length(dtmtable.trnidx) #0.30466923428172205
valAvgError=sum(errs[dtmtable.validx])/length(dtmtable.validx) #0.3192841010476779


#Unseen data
#Note: For any out of sample data (which MUST have exactly the same structure as dtmtable.features!)
#you can also use the predict function
estimatedFrequencieForFirst20Obs=predict(resM,dtmtable.features[1:20,:])
#getting new/unseen data into the right format is not trivial at the moment. It is easiset to use prepare_dataframe_for_dtm! to the 
#whole data set and then discrard possible hold out observations  (or mark as validation through validx)



############################################################
#Additional options
############################################################

#You can modify the settings directly
sett.niter=4
#however it is preferred to use the update function to avoid creating an 'invalid' setting (such as a negative number of iterations)
updateSettingsMod!(sett,niter=-3) #e.g. this will throw an error (whereas sett.niter=-3 does not)


#SOME COMMENT ON MANY UNUSED AND EXPERMIANTAL options
#if you want to use the resutling tree in a different programming language you may want to
#consider these options which give you the core structure of the model in other languages
# updateSettingsMod!(sett,write_csharp_code="true",write_vba_code="true",write_sas_code="true")

#mention DOT graph

############################################################
#Prepare a grid search over the parameter space
############################################################

#set some default settings
updateSettingsMod!(sett,
niter=100,
model_type="boosted_tree",
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

#performanceMeasure will determine which metric is considered for the output summary 
#performanceMeasure must be in DecisionTrees.GLOBALperfMeasursesWhereMaxIsBest
#you may want to consider the Excel file of any boosting to get an idea of the available options for performanceMeasure 
#consider the first row of the sheet ModelStatistics of the Excel file which lists all measures
#We note that some missing (or constant) measures are either not fully implemented yet, or were not calculated (e.g. because a boolean in the settings was set to false)
@show sett.performanceMeasure
@assert in(sett.performanceMeasure,DecisionTrees.GLOBALperfMeasursesWhereMaxIsBest)

#create a vector of ModelSettings
#note, we can provide an arbitrary number of vectors over which we want to loop
#the function createGridSearchSettings will copy the 'input' settings X times into a vector of settings
#X is the size of the cartesian product of all parameter vectors over which we want to loop.

settV=createGridSearchSettings(sett,    
    minw=[-0.2,-0.1,-0.05,-0.025,-0.01]
    ,mf=[0.01,0.05,0.1],
    subsampling_features_prop=[1.0,.5],
    subsampling_prop=[1.0,.5],
    smoothEstimates=["0"]);

#consider the length(settV) which is the number of models that will be run
length(settV)

#depending on the size of settV, you may want to restart Julia with additional workers for parallelization, i.e. 
#shell> julia -p 8 
#to start julia with 8 workers. Depending on your CPU you may want to choose a different number

Sys.CPU_CORES #might be give you an indication of the number of workers() you could use

############################################################
#Run grid search
############################################################

@warn "This may take quiet some time depending on length(settV)=$(length(settV))"
@info "Starting grid search..."

tt0=time_ns()
#dtm(dtmtable,settV,fn="R:\\temp\\1\\dtm.CSV")
@show ela=(-tt0+time_ns())/1e9
@info ".....done"

warn("todo:check best tree with poisson error too!")
    #crit::SplittingCriterion # fn version of 4
