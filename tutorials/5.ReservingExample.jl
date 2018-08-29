#Date: July 1, 2018
#Author: Bernhard König
#Title: An application of Regression Trees for Reserving 

#We have used the code in RCode\GenerateIndividualClaimsData.R 
#to generate the data Simulated.Cashflows.csv
#We refer to the work of Wuethrich and Gabrielli for details on the
#simulated claims data, see the following references:
    ### Source: https://people.math.ethz.ch/~wmario/simulation.html
    ### Reference: Individual Claims History Simulation Machine. ###
    ###            SSRN Manuscript ID 3130560.                   ###  
    ###  Authors: Andrea Gabrielli, Mario V. Wuthrich    ###

#
@assert splitdir(pwd())[2]=="DecisionTrees.jl" """Error. The working directory should be the DecisionTrees.jl folder.\r\nUse 'cd("C:\\Users\\bernhard.konig\\Documents\\ASync\\home\\Code\\Julia\\DecisionTrees.jl")' or similar to change the working dir in Julia"""
#if the working dir is not correct, the include command further down will fail 

@warn("You need to install DecisionTrees.jl. Consider the readme of this page: 'https://github.com/kafisatz/DecisionTrees.jl'")
@warn("You need to make sure that the required packages are are installed for this code to run.")

using Revise
#stdlib packages:
    using Statistics
    using Distributed
    using Random
    using DelimitedFiles 
    using Pkg

Distributed.@everywhere import CSV
Distributed.@everywhere import DataFrames
Distributed.@everywhere import DataFrames: DataFrame,groupby,combine,names!,aggregate

@time Distributed.@everywhere using DecisionTrees #first time usage (precompilation) may take some time here 

#folderForOutput="H:\\Privat\\"
folderForOutput="C:\\Users\\bernhard.konig\\Documents\\ASync\\publicnl\\Personal\\Bernhard\\Projects & Meetings\\2018 SAV Vortrag MV\\ReservingTreesData\\20180829\\"
@assert isdir(folderForOutput) "Directory does not exist: $(folderForOutput)"
#define functions
#include(joinpath(@__DIR__,"..","tutorials","5.ReservingExample_functions.jl"))
#include(joinpath(Pkg.dir("DecisionTrees"),"tutorials","5.ReservingExample_functions.jl")) #this may not use the most current file but rather the file in the folder with compiled julia code, e.g. C:\\Users\\bernhard.konig\\.julia\\packages\\DecisionTrees\\3AXAR\\.... or similar
fn_file=joinpath(pwd(),"tutorials","5.ReservingExample_functions.jl") 
@assert isfile(fn_file) 
include(fn_file)

##############################
#Read the data
##############################
#note: the underyling data was generated with R (consider the readme.txt in the folder data\\reservingAccident\\)
datafile=joinpath("data","reservingAccident","Simulated.Cashflows.csv")
@assert isfile(datafile) "File not found: $(datafile). Consider the readme.txt in the folder data\\reservingAccident\\"
elt=[String,String,String,Int,Int,Int,String,Int]
elt=vcat(elt,repeat([Float64],12))
@info("Reading data...")
@time fullData=CSV.read(datafile,types=elt,rows_for_type_detect=100000,allowmissing=:none,categorical=false);

#for development only: define random subset of data
    debug=false
    if debug 
        @warn("Subsetting data:")
        rnd=rand(size(fullData,1))
        srtTmp=sortperm(rnd)
        subsetIdx=srtTmp[1:200000]
        fullData=fullData[subsetIdx,:]
        folderForOutput="C:\\temp\\"
    end

#Define a random split into training and validation data 
#training shall be 70% of the data
Random.seed!(1239)
trnValCol=BitArray(rand() < 0.7 for x in 1:size(fullData,1));
fullData[:trnValCol]=trnValCol

@warn("We may want to redo the whole calculation with no validation data -> the CL factors will be based on all of the data (not only the training porportion).")

#The data should match this total for the column Paycum11
expectedSum=9061763108
if !debug
    @assert sum(fullData[:PayCum11])==expectedSum "Expected $(expectedSum), have: $(sum(fullData[:PayCum11]))"
end
#The data is described in the paper mentioned above

#We provide a short summary from the Readme of Wuethrich and Gabrielli
#The variable ClNr is the claim number, this is a unique numerical identifier for each individual claim. 
#LoB is the line of business which should be of factor type (categorical) and may take the labels 1; : : : ; 4. 
#cc is the claims code that should be of factor type (categorical) and may take the labels 1; : : : ; 53 (with possible gaps).
#AY is the year of claims occurrence (accident year) that should be of integer type and may take the values 1994; : : : ; 2005. 
#AQ is the quarter of claims occurrence (accident quarter) that should be of integer type and may take the values 1; : : : ; 4.
#age is the age of the injured that should be of integer type and may take the values 15; : : : ; 70.
#inj_part is the body part injured that should be of factor type (categorical) and may take the labels 1; : : : ; 99 (with gaps).

#set independent variables
selected_explanatory_vars=["AQ","LoB","age","cc","inj_part"]
categoricalVars=["LoB","inj_part","cc"]

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
    println(names(fullData)[x])
    @show size(unique(fullData[:,x]))
end

#=
    #Let us see whether all claims are 'increasing' in time (i.e. whether all incremental payments are non-negative)
    #NOTE: This may take a few minutes to compute!
    #->there are some rows/claims which have non increasing payments over time
    #99.7% of the claims are nondecreasing though

    isNonDecreasing=checkIfPaymentsAreIncreasing(fullData[1:30]) 
    sum(isNonDecreasing)/length(isNonDecreasing) # 99.7% of the claims are nondecreasing
    nonDecreasingRows=.!isNonDecreasing
    nonDecreasingRowsInteger=find(nonDecreasingRows)
    @show size(nonDecreasingRowsInteger)
    #claims which are decreasing over time (at some point):
    #fullData[nonDecreasingRowsInteger,:]
=#


##############################
#Prepare the data
##############################
#consider paid to date amount
paidToDateFullDataPerRow=getPaidToDatePerRow(fullData)
fullData[:paidToDate]=paidToDateFullDataPerRow

#consider the Accident Years
ayears=sort(unique(fullData[:AY]))
maxAY=maximum(ayears)
@assert maxAY==2005

#prepare the whole dataset in the same manner as the training data set
    #indicator index whether claim is known by year end 2005
    indexOfClaimIsReportedByYE2005=fullData[:AY]+fullData[:RepDel].<=2005
    dataKnownByYE2005=fullData[indexOfClaimIsReportedByYE2005,:]
    #paycumcolumns=[Symbol(string("PayCum",lpad(ii,2,0))) for ii=0:11]
    #map(x->sum(fullData[.!indexOfClaimIsReportedByYE2005,x]),paycumcolumns)
    #aggregate(df,:someColumn,[sum, mean])    
    dtmKnownByYE2005,sett,dfpreppedUnused=prepare_dataframe_for_dtm!(dataKnownByYE2005,keycol="ClNr",numcol="PayCum00",trnvalcol="trnValCol",independent_vars=selected_explanatory_vars,treat_as_categorical_variable=categoricalVars);
    #define general model settings
    updateSettingsMod!(sett,ignoreZeroDenominatorValues=true,model_type="build_tree",write_dot_graph=true,writeTree=false,graphvizexecutable="C:\\Program Files (x86)\\Graphviz2.38\\bin\\dot.exe")
    #a note on the data
    #dataKnownByYE2005 is a DataFrame (it contains more data such as all PayCum columns and the AY)
    #dtmKnownByYE2005 is a DTMTable object (it has a different structre than the DataFrame)    

#consider paidToDate and 'truth'
paidToDatePerRow=getPaidToDatePerRow(dataKnownByYE2005)
#trueUltimatePerRow (this is not meaningful because we have two sets of claims: fullData[:PayCum11] AND dtmKnownByYE2005[:PayCum11])
#trueReservePerRow=trueUltimatePerRow.-paidToDatePerRow #this is not meaningful either

#build 'triangle' (in this case, we will actually get a rectangle rather than a triangle)
@info("Fitting 'plain' CL to all data (including 'unknown' part of the triangle)")
triangleInclIBNyRClaims=buildTriangle(fullData)
trueUltimate=triangleInclIBNyRClaims[:,end]

#CL
#consider aggregate CL Model
@info("Fitting 'plain' CL to all data")
triangle=buildTriangle(dataKnownByYE2005)
clAllLOBs=chainLadder(triangle)

trueReserves=trueUltimate.-clAllLOBs[:paidToDate]

#Consider CL Ultimate per known claim
#LDF = Loss Development Factor (also chain ladder factor)
ayperRow=dataKnownByYE2005[:AY]
CLfactorsToUltimate=clAllLOBs[:factorsToUltimate]
CLLDFperRow=map(x->CLfactorsToUltimate[x-1993],ayperRow)
CLUltimatePerRow=CLLDFperRow.*paidToDatePerRow
cldf=DataFrame(ldf=CLLDFperRow,ultimate=CLUltimatePerRow,paidToDate=paidToDatePerRow)
cldf[:reserves]=cldf[:ultimate].-cldf[:paidToDate]
CLcomparisons=CL_est_per_variable(fullData,cldf)
summaryByVar=vcat(collect(values(CLcomparisons)))
CSV.write(string(folderForOutput,"resultsCL.csv"),summaryByVar)

if !debug 
    @assert isapprox(0,sum(CLUltimatePerRow)-sum(clAllLOBs[:ultimate]),atol=1e-4) #should be zero 
    for i in ayears
        idx=ayperRow.==i
        @assert isapprox(0,sum(CLUltimatePerRow[idx])-clAllLOBs[:ultimate][i-1993],atol=1e-4) #should be zero for each year
    end
end
CLTotalUltimate=sum(CLUltimatePerRow)

#consider CL per LoB
@info("Fitting 'plain' CL models to each LoB")
LOB="1";triangleLOB1=buildTriangle(dataKnownByYE2005[dataKnownByYE2005[:LoB].==LOB,:]);clLOB1=chainLadder(triangleLOB1);trueUltimateLOB1=buildTriangle(fullData[fullData[:LoB].==LOB,:])[:,end];clLOB1[:trueReserves]=trueUltimateLOB1.-clLOB1[:paidToDate];clLOB1[:error]=clLOB1[:trueReserves].-clLOB1[:reserves];
LOB="2";triangleLOB2=buildTriangle(dataKnownByYE2005[dataKnownByYE2005[:LoB].==LOB,:]);clLOB2=chainLadder(triangleLOB2);trueUltimateLOB2=buildTriangle(fullData[fullData[:LoB].==LOB,:])[:,end];clLOB2[:trueReserves]=trueUltimateLOB2.-clLOB2[:paidToDate];clLOB2[:error]=clLOB2[:trueReserves].-clLOB2[:reserves];
LOB="3";triangleLOB3=buildTriangle(dataKnownByYE2005[dataKnownByYE2005[:LoB].==LOB,:]);clLOB3=chainLadder(triangleLOB3);trueUltimateLOB3=buildTriangle(fullData[fullData[:LoB].==LOB,:])[:,end];clLOB3[:trueReserves]=trueUltimateLOB3.-clLOB3[:paidToDate];clLOB3[:error]=clLOB3[:trueReserves].-clLOB3[:reserves];
LOB="4";triangleLOB4=buildTriangle(dataKnownByYE2005[dataKnownByYE2005[:LoB].==LOB,:]);clLOB4=chainLadder(triangleLOB4);trueUltimateLOB4=buildTriangle(fullData[fullData[:LoB].==LOB,:])[:,end];clLOB4[:trueReserves]=trueUltimateLOB4.-clLOB4[:paidToDate];clLOB4[:error]=clLOB4[:trueReserves].-clLOB4[:reserves];

addTotal!(clLOB1)
addTotal!(clLOB2)
addTotal!(clLOB3)
addTotal!(clLOB4)

trueReservesPerLOB=Dict{String,Array{Float64}}()
trueReservesPerLOB["1"]=clLOB1[:trueReserves]
trueReservesPerLOB["2"]=clLOB2[:trueReserves]
trueReservesPerLOB["3"]=clLOB3[:trueReserves]
trueReservesPerLOB["4"]=clLOB4[:trueReserves]

clPerLOB=Dict{String,DataFrame}()
clPerLOB["1"]=clLOB1
clPerLOB["2"]=clLOB2
clPerLOB["3"]=clLOB3
clPerLOB["4"]=clLOB4

if !debug 
    @assert all(isapprox.(vcat(trueReserves,sum(trueReserves)).-clLOB1[:trueReserves].-clLOB2[:trueReserves].-clLOB3[:trueReserves].-clLOB4[:trueReserves],0))
    @assert all(isapprox.(trueUltimate.-trueUltimateLOB1.-trueUltimateLOB2.-trueUltimateLOB3.-trueUltimateLOB4,0))
end
##############################
#Consider CL decision trees
##############################

#Define minimum weight per leaf (depending on the development year)
#note:these weights are calibrated for a total training weight of 3477583 (Notably the number of claims will be smaller for higher LDFs!)
@show sum(dataKnownByYE2005[:trnValCol]) #expected value: 3477583
@assert debug||(sum(dataKnownByYE2005[:trnValCol])==3477583)

#each entry corresponds to an LDF Year
modelsWeightsPerLDF=Vector{Vector{Float64}}()
#the values were selected by expert judgement
#primarily we considered reversals (i.e when the observed ratio in the validation data is not monotonic in the number of leaves (as the leaf number is defined such that the training observed ratios are increasing))
push!(modelsWeightsPerLDF,[150000,200000,270000,220000,235000,210000,200000,240000,130000,130000,15000000])

#alternatives

#old weights for sse
#push!(modelsWeightsPerLDF,[130000,170000,240000,220000,235000,210000,200000,240000,130000,130000,15000000])

#You can provide additional sets of minimum weights to see the differences between two models
#push!(modelsWeightsPerLDF,[70000,170000,200000,220000,235000,210000,200000,240000,200000,130000,15000000])

#WEIGHTS FOR DIFFERENCE FN: push!(modelsWeightsPerLDF,[130000,170000,240000,220000,235000,210000,200000,240000,130000,130000,15000000])
#version where the weight is increasing
#push!(modelsWeightsPerLDF,[80000,170000,240000,220000,240000,240000,240000,240000,240000,240000,15000000])

for i=1:length(modelsWeightsPerLDF)
    @assert length(modelsWeightsPerLDF[i])==length(ayears)-1 # == number of LDFs
end

#goal: build a tree for each LDF (i.e. a tree for 12-24 months, a second tree for 24-36 months, etc.)

treeResults=Dict{Float64,DataFrame}()
treeResultsAgg=Dict{Float64,DataFrame}()

@info("Starting model fit...")
@show size(modelsWeightsPerLDF)
#folderForOutput should be defined further above
@assert isdir(folderForOutput) "Directory does not exist: $(folderForOutput)"
#actual model run (this may take a few minutes)
updateSettingsMod!(sett,crit="sse",model_type="build_tree")
@time this_ldfarr=runModels!(dataKnownByYE2005,dtmKnownByYE2005,modelsWeightsPerLDF,treeResults,treeResultsAgg,selected_explanatory_vars,categoricalVars,folderForOutput,clAllLOBs,paidToDatePerRow,sett);

@info("Finished building models.")

mdl_stats=write_results(treeResults,treeResultsAgg,folderForOutput)


#boosted model
@info("Starting Boosting models...")

treeResultsBoosting=Dict{Float64,DataFrame}()
treeResultsAggBoosting=Dict{Float64,DataFrame}()
updateSettingsMod!(sett,crit="sse",iterations=0*2+1*8,learningRate=.15,model_type="boosted_tree",subsampling_features_prop=0.6)
@time ldfArrBoosting=runModels!(dataKnownByYE2005,dtmKnownByYE2005,modelsWeightsPerLDF,treeResultsBoosting,treeResultsAggBoosting,selected_explanatory_vars,categoricalVars,folderForOutput,clAllLOBs,paidToDatePerRow,sett);

mdl_statsBoosting=write_results(treeResultsBoosting,treeResultsAggBoosting,folderForOutput)

@warn("If we were to update the 'tree-CL' factors such that they are based on all data (currently they are only using the training data), the model might further improve.")
@warn("Make a note that 5m claims is a large data set! In practice less data might be availabe but the concept can be applied nevertheless")

#=    
#try some alternative models
NN=20
this_ldfarr[1:NN,1]
ldfArr[1:NN,1]
    ldfYear=4
    selectedWeight=250000
    settC=deepcopy(sett)
    updateSettingsMod!(settC,crit="sse")
    resultingFiles,resM= runSingleModel(dataKnownByYE2005,dtmKnownByYE2005,selectedWeight,ldfYear,2005,               selected_explanatory_vars,categoricalVars,"C:\\temp\\331\\",settC)
    0


#tbd: compare mse versus sse
selectedWeight=150000
while (selectedWeight < 180000)
    ldfYear=2
    #selectedWeight=130000
    settC=deepcopy(sett)
    updateSettingsMod!(settC,crit="mse")
    resultingFiles,resM= runSingleModel(dataKnownByYE2005,dtmKnownByYE2005,selectedWeight,ldfYear,2005,               selected_explanatory_vars,categoricalVars,"C:\\temp\\331\\",settC)

    settC=deepcopy(sett)
    updateSettingsMod!(settC,crit="sse")
    resultingFiles,resM= runSingleModel(dataKnownByYE2005,dtmKnownByYE2005,selectedWeight,ldfYear,2005,               selected_explanatory_vars,categoricalVars,"C:\\temp\\331\\",settC)
    selectedWeight+=10000
end
   
ldfYear=1
selectedWeight=170000
settC=deepcopy(sett)
updateSettingsMod!(settC,crit="difference",iterations=8,learningRate=.15,model_type="boosted_tree",subsampling_features_prop=0.6)
resultingFiles,resM= runSingleModel(dataKnownByYE2005,dtmKnownByYE2005,selectedWeight,ldfYear,2005,               selected_explanatory_vars,categoricalVars,"C:\\temp\\331\\",settC)
dfpred=predict(resM,dtmKnownByYE2005.features)  
fittedValues=dfpred[:RawEstimate]
@show mean(fittedValues),median(fittedValues)
0

ldfYear=1
selectedWeight=130000
settC=deepcopy(sett)
updateSettingsMod!(settC,crit="gamma")
resultingFiles,resM= runSingleModel(dataKnownByYE2005,dtmKnownByYE2005,selectedWeight,ldfYear,2005,               selected_explanatory_vars,categoricalVars,"C:\\temp\\331\\",settC)
0

    =#


#=
#test algorithm on smaller data Set

if debug 
    @warn("Subsetting data:")
    this_ldfy=1
    this_selectedWeight=Int(round(1/size(fullData,1)*modelsWeightsPerLDF[1][this_ldfy]*trnsize,digits=0))
    this_selectedWeight=5500
    trnsize=40000
    cvsampler=CVOptions(0,0.0,true)
    test_on_subset(dataKnownByYE2005,this_ldfy,folderForOutput,sett,trnsize,this_selectedWeight,cvsampler)

    cvsampler=CVOptions(-2,0.0,true)
    test_on_subset(dataKnownByYE2005,this_ldfy,folderForOutput,sett,trnsize,this_selectedWeight,cvsampler)
0
end

=#