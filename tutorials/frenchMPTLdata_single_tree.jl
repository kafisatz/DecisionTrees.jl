############################################################
#Decision Tree Model for the French MTPL Data
############################################################

# If you do not yet have a suitable editor for Julia, you may want to consider VS Code.
# Consider this link
# https://github.com/JuliaEditorSupport/julia-vscode
# We note that there are various alternative editors for Julia.

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

#boolCalculatePoissonError, is a boolean which we need to enable here in order for output to show the poisson error of the estimates
updateSettingsMod!(sett,minw=-0.03,model_type="build_tree",boolCalculatePoissonError=true)

##############################
#Run single tree model
##############################

#resM is the resulting Model
resultingFiles,resM=dtm(dtmtable,sett)
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

#the leaves are numbered in such a manner that the lowest fitted value corresponds to leaf number 1
#conversly, the highest fitted value corresponds to the highest leaf number.

############################################################
#Additional options
############################################################

#If you have Graphviz installed, you can generate a plot (on a PDF) of the tree
updateSettingsMod!(sett,write_dot_graph="true",graphvizexecutable="c:\\Program Files (x86)\\Graphviz2.38\\bin\\dot.exe")
resultingFiles,resM=dtm(dtmtable,sett)

#You will also get a *.dot.txt file with a text representation of the graph

#WE NOTE THAT THERE ARE A NUMBER OF UNDOCUMENTED OPTIONS. MOST OF THESE ARE EXPERIMENTAL OR UNUSED. 

#if you want to use the resutling tree in a different programming language you may want to
#consider these options which give you the core structure of the model in other languages
# updateSettingsMod!(sett,write_csharp_code="true",write_vba_code="true",write_sas_code="true")
# some of these options may only be applicable to boosting OR single tree models

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

minweight_list=-collect(range(0.001,stop=0.5,length=20))

settV=createGridSearchSettings(sett,    
    minw=minweight_list
    );

#consider the length(settV) which is the number of models that will be run
length(settV)

#depending on the size of settV (and the size of the data set!), you may want to restart Julia with additional workers for parallelization, i.e. 
#shell> julia -p 8 
#to start julia with 8 workers. Depending on your CPU you may want to choose a different number

Sys.CPU_CORES #might be give you an indication of the number of workers() you could use

############################################################
#Run grid search
############################################################

tt0=time_ns()
dtm(dtmtable,settV,fn="R:\\temp\\1\\dtm.CSV")
@show ela=(-tt0+time_ns())/1e9
@info ".....done"

warn("todo:check best tree with poisson error too!")
    #crit::SplittingCriterion # fn version of 4
