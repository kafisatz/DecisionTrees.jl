#Date: July 3, 2018
#Author: Bernhard König

############################################################
#Title: Decision Tree models used for the Whitepaper
############################################################

# If you do not yet have a suitable editor for Julia, you may want to consider VS Code.
# Consider this link
# https://github.com/JuliaEditorSupport/julia-vscode
# We note that there are various alternative editors for Julia.

##############################
#Load packages
##############################
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

#the data is described in the R package CASdatasets
#also the Noll, Salzmann, Wüthrich Paper has some descriptive graphs of the data.

#We will slightly modify the data in the following
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
    println("Variable number ",x,": ",names(fullData)[x])
    @show size(unique(fullData[:,x]))[1]
end

##############################
#Prepare the data
##############################

#dtmtable::DTMTable is a custom data format for DecisionTrees
#sett are the model settings
#dfprepped is an intermediary dataframe which is generally not needed
dtmtable,sett,dfprepped=prepare_dataframe_for_dtm!(fullData,keycol="IDpol",trnvalcol="trnTest",numcol="ClaimNb",denomcol="Exposure",weightcol="Exposure",independent_vars=selected_explanatory_vars);

#=
#NOTE:if the input data (i.e. the DataFrame fullData in this case) contains a row which is currently 'numerical' but should be treated as a categorical variable,
#then you can use the parameter treat_as_categorical_variable
#consider this code example, which is NOT MEANINGFUL in this context, but showcases the optional parameter
    dtmtableAlt,settAlt,dfpreppedAlt=prepare_dataframe_for_dtm!(fullData,keycol="IDpol",trnvalcol="trnTest",numcol="ClaimNb",denomcol="Exposure",weightcol="Exposure",independent_vars=selected_explanatory_vars,treat_as_categorical_variable=["DrivAge","VehAge"]);
    dtmtable.features[:DrivAge].pool #the data was interpreted as a numerical feature
    dtmtableAlt.features[:DrivAge].pool #the data was interpreted as a categorical feature (i.e. a String)
=#

#Consider the following regarding the data structure of dtmtable
fieldnames(typeof(dtmtable))
dtmtable.key #an identifier (String type), ideally it is a unique identifier for each row
dtmtable.numerator #a vector
dtmtable.denominator #a vector
dtmtable.weight #a vector
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
#We realize that it may be sutiable to define the splitting in a different manner (than uniformly spaced).
#This feature might be added at a later time. You can consider the function add_coded_numdata!(...) to see how splitting points are chosen

#We note that 'pooling' the data into X (up to sett.max_splitting_points_num) buckets has a few advantages:
#one can do sorting in linear time (countsort)
#the data takes up less space in RAM and caches which generally increases performance

#keep a copy of original train and test split
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

#boolCalculatePoissonError, is a boolean which we need to enable here in order for output to show the poisson error of the estimates
updateSettingsMod!(sett,minw=-0.03,model_type="build_tree",boolCalculatePoissonError=true)

##############################
#Run single tree model
##############################

#resM is the resulting Model
resultingFiles,resM=dtm(dtmtable,sett)
#time to build the tree (3% min exposure, difference measure): 0.8s

#consider the resulting Excel file for a description of the model 
#it should be located here:
@show resultingFiles[1]

############################################################
#Investigate the predictions (fitted values)
############################################################

#For a single decision tree, the mean observed value for each leaf is one choice for the predicted value.

fitted,leafnumbers=predict(resM,dtmtable.features)

#by definition there are only N distinct fitted values, where N = 'number of leaves'
unique(fitted)
unique(leafnumbers)

#todo: add the poisson error calc here in julia
#notably, we are using the difference as splitting criterion, which has (to our knowledge) no explicit relationship
#to the poisson distribution which is often assumed for count (and frequency) data.

#the leaves are numbered in such a manner that the lowest fitted value corresponds to leaf number 1
#conversly, the highest fitted value corresponds to the highest leaf number.

############################################################
#Additional options
############################################################

#Changing the impurity function
#this is detailed further below

#Tree representation
#If you have Graphviz installed, you can generate a plot (on a PDF) of the tree
updateSettingsMod!(sett,write_dot_graph="true",graphvizexecutable="c:\\Program Files (x86)\\Graphviz2.38\\bin\\dot.exe")
resultingFiles,resM=dtm(dtmtable,sett)
#You will also get a *.dot.txt file with a text representation of the graph

#If you do not have graphviz installed, you can do the following: 
    #you can create a tree representation like this
    dot_graph=graph(resM) #this is the content of the .dot.txt file

    #copy tree graph to Clipboard
    clipboard(dot_graph)
    #go to 'a graphviz online too':
    https://dreampuf.github.io/GraphvizOnline/ 
    #paste the tree 


#output location
#you can provide a file and folder name in the following manner:
resultingFiles,resM=dtm(dtmtable,sett,file="C:\\temp\\somename.txt")
#the extension is not relevant (but we suggest to provide it)
#you should ensure that the folder exists, otherwise the model will not run.

#you can set to disable much of the text in the julia window (i.e. the log in stdout which gets printed by the model)
sett.print_details=false
#notably thenames of the output files will still be shown.

#let us keep it at true for the following
sett.print_details=true

#WE NOTE THAT THERE ARE A NUMBER OF UNDOCUMENTED OPTIONS. MOST OF THESE ARE EXPERIMENTAL OR UNUSED. 

#if you want to use the resutling tree in a different programming language you may want to
#consider these options which give you the core structure of the model in other languages
# updateSettingsMod!(sett,write_csharp_code="true",write_vba_code="true",write_sas_code="true")
# some of these options may only be applicable to boosting OR single tree models

#Consider the Poisson splitting criterion
settPoisson=deepcopy(sett)
updateSettingsMod!(settPoisson,crit="poisson")
dtm(dtmtable,settPoisson) #time to build tree: 12s
#time to build the tree (3% min exposure, Poisson deviance): about 12s


############################################################
#Prepare a grid search over the parameter space
############################################################

#set some default settings
updateSettingsMod!(sett,
model_type="build_tree",
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
#performanceMeasure must be in DecisionTrees.global_statsperiter_header
#you may want to consider the Excel file of any boosting to get an idea of the available options for performanceMeasure 
#consider the first row of the sheet ModelStatistics of the Excel file which lists all measures
#We note that some missing (or constant) measures are either not fully implemented yet, or were not calculated (e.g. because a boolean in the settings was set to false)
@show sett.performanceMeasure
@assert in(sett.performanceMeasure,DecisionTrees.global_statsperiter_header)

#create a vector of ModelSettings
#note, we can provide an arbitrary number of vectors over which we want to loop
#the function createGridSearchSettings will copy the 'input' settings X times into a vector of settings
#X is the size of the cartesian product of all parameter vectors over which we want to loop.

#Let us consider different stopping criteria (minimum weight per leaf)
#minweight_list=-collect(range(0.0001,stop=0.1,length=25))
minweight_list=Float64[50, 100, 150, 200, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 6000, 7000, 8000, 9000, 10000, 12500, 15000, 17500, 20000, 22500, 25000]

settV=createGridSearchSettings(sett,    
    minw=minweight_list
    );

#consider the length(settV) which is the number of models that will be run
length(settV)

#notably all these models will be constructed independently of each other#
#if we are only varying the minimum weight this could of course be improved (by pruning a granular tree)
#however, the functionality is much more flexible (see the boosting example). Therefore make no assmptions on the
#relationship between settV[1] and settV[2], etc.

#depending on the size of settV (and the size of the data set!), you may want to restart Julia with additional workers for parallelization, i.e. 
#shell> julia -p 8 
#to start julia with 8 workers. Depending on your CPU you may want to choose a different number

Sys.CPU_CORES #might be give you an indication of the number of workers() you could use

############################################################
#Run grid search
############################################################

outputfolderandFile="R:\\temp\\1\\MTPLsingleTreeDifference.CSV"
outputfolder=splitdir(outputfolderandFile)[1]
@assert isdir(outputfolder)

tt0=time_ns()
gridResult=dtm(dtmtable,settV,file=outputfolderandFile)
@show ela=(-tt0+time_ns())/1e9


# The excelfile 'MTPLsingleTreeDifference_multistats.xlsx' will give you a summary of the statistics of each model run
# The validation error may not be evaluated for the most granular tree
# This is because the tree has found leaves with a fitted frequency of zero (on the training data)


#Perform grid search for Poisson
settPoissonV=createGridSearchSettings(settPoisson,    
    minw=minweight_list
    );

outputfolderandFile2="R:\\temp\\2\\MTPLsingleTreePoisson.CSV"
outputfolder2=splitdir(outputfolderandFile2)[1]
@assert isdir(outputfolder2)

@warn("This may take quite some time!")
tt0=time_ns()
gridResult=dtm(dtmtable,settPoissonV,file=outputfolderandFile2)
@show ela=(-tt0+time_ns())/1e9

    
