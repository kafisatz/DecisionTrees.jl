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

import Distributed
import Random  
import DelimitedFiles
Distributed.@everywhere import CSV
Distributed.@everywhere import DataFrames
Distributed.@everywhere import DataFrames: DataFrame,groupby,combine,names!,aggregate

@time Distributed.@everywhere using DecisionTrees #200 - 300 seconds in 0.7beta INCLUDING precompilation (124s w/o precompilation) times are for my notebook (precompilation allocates 20GB of memory?)

#folderForOutput="C:\\temp\\"
folderForOutput="H:\\Privat\\SAV\\Fachgruppe Data Science\\ReservingTrees\\20180702\\"
@assert isdir(folderForOutput) "Directory does not exist: $(folderForOutput)"
#define functions
#=
    include(joinpath(@__DIR__,"tutorials","5.ReservingExample_functions.jl"))
=#
#include(joinpath(@__DIR__,"..","tutorials","5.ReservingExample_functions.jl"))
include(joinpath(Pkg.dir("DecisionTrees"),"tutorials","5.ReservingExample_functions.jl"))

##############################
#Read the data
##############################
datafile=joinpath("data","reservingAccident","Simulated.Cashflows.csv")
@assert isfile(datafile) "File not found: $(datafile). Consider the readme.txt in the folder data\\reservingAccident\\"
elt=[String,String,String,Int,Int,Int,String,Int]
elt=vcat(elt,repeat([Float64],12))
@time fullData=CSV.read(datafile,types=elt,rows_for_type_detect=100000,allowmissing=:none,categorical=false);

#Define a random split into training and validation data 
#training shall be 70% of the data
Random.srand(1239)
trnValCol=BitArray(rand() < 0.7 for x in 1:size(fullData,1));
fullData[:trnValCol]=trnValCol

@warn("We may want to redo the whole calculation with no validation data -> the CL factors will be based on all of the data (not only the training porportion).")

#The data should match this total for the column Paycum11
expectedSum=9061763108 
@assert sum(fullData[:PayCum11])==expectedSum "Expected $(expectedSum), have: $(sum(fullData[:PayCum11]))"
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

#consider the Accient Years
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
    dtmtableKnownData,settUnused,dfpreppedUnused=prepare_dataframe_for_dtm!(dataKnownByYE2005,keycol="ClNr",numcol="PayCum00",trnvalcol="trnValCol",independent_vars=selected_explanatory_vars,treat_as_categorical_variable=categoricalVars);

#consider paidToDate and 'truth'
paidToDatePerRow=getPaidToDatePerRow(dataKnownByYE2005)
#trueUltimatePerRow (this is not meaningful because we have two sets of claims: fullData[:PayCum11] AND dtmtableKnownData[:PayCum11])
#trueReservePerRow=trueUltimatePerRow.-paidToDatePerRow #this is not meaningful either

#build 'triangle' (in this case, we will atualls get a retangle rather than a triangle)
triangleInclIBNyRClaims=buildTriangle(fullData)
trueUltimate=triangleInclIBNyRClaims[:,end]

#CL
#consider aggregate CL Model
triangle=buildTriangle(dataKnownByYE2005)
clAllLOBs=chainLadder(triangle)

trueReserves=trueUltimate.-clAllLOBs[:paidToDate]

#Consider CL Ultimate per known claim
#LDF = Loss Development Factor (also chain ladder factor)
ayperRow=dataKnownByYE2005[:AY]
CLfactorsToUltimate=clAllLOBs[:factorsToUltimate]
CLLDFperRow=map(x->CLfactorsToUltimate[x-1993],ayperRow)
CLUltimatePerRow=CLLDFperRow.*paidToDatePerRow
@assert isapprox(0,sum(CLUltimatePerRow)-sum(clAllLOBs[:ultimate]),atol=1e-4) #should be zero 
for i in ayears
    idx=ayperRow.==i
    @assert isapprox(0,sum(CLUltimatePerRow[idx])-clAllLOBs[:ultimate][i-1993],atol=1e-4) #should be zero for each year
end
CLTotalUltimate=sum(CLUltimatePerRow)

#consider CL per LoB
LOB="1";triangleLOB1=buildTriangle(dataKnownByYE2005[dataKnownByYE2005[:LoB].==LOB,:]);clLOB1=chainLadder(triangleLOB1);trueUltimateLOB1=buildTriangle(fullData[fullData[:LoB].==LOB,:])[:,end]
LOB="2";triangleLOB2=buildTriangle(dataKnownByYE2005[dataKnownByYE2005[:LoB].==LOB,:]);clLOB2=chainLadder(triangleLOB2);trueUltimateLOB2=buildTriangle(fullData[fullData[:LoB].==LOB,:])[:,end]
LOB="3";triangleLOB3=buildTriangle(dataKnownByYE2005[dataKnownByYE2005[:LoB].==LOB,:]);clLOB3=chainLadder(triangleLOB3);trueUltimateLOB3=buildTriangle(fullData[fullData[:LoB].==LOB,:])[:,end]
LOB="4";triangleLOB4=buildTriangle(dataKnownByYE2005[dataKnownByYE2005[:LoB].==LOB,:]);clLOB4=chainLadder(triangleLOB4);trueUltimateLOB4=buildTriangle(fullData[fullData[:LoB].==LOB,:])[:,end]

trueReservesLOB1=trueUltimateLOB1.-clLOB1[:paidToDate]
trueReservesLOB2=trueUltimateLOB2.-clLOB2[:paidToDate]
trueReservesLOB3=trueUltimateLOB3.-clLOB3[:paidToDate]
trueReservesLOB4=trueUltimateLOB4.-clLOB4[:paidToDate]

trueReservesPerLOB=Dict{String,Array{Float64}}()
trueReservesPerLOB["1"]=trueReservesLOB1
trueReservesPerLOB["2"]=trueReservesLOB2
trueReservesPerLOB["3"]=trueReservesLOB3
trueReservesPerLOB["4"]=trueReservesLOB4

clPerLOB=Dict{String,DataFrame}()
clPerLOB["1"]=clLOB1
clPerLOB["2"]=clLOB2
clPerLOB["3"]=clLOB3
clPerLOB["4"]=clLOB4

@assert all(isapprox.(trueReserves.-trueReservesLOB1.-trueReservesLOB2.-trueReservesLOB3.-trueReservesLOB4,0))
@assert all(isapprox.(trueUltimate.-trueUltimateLOB1.-trueUltimateLOB2.-trueUltimateLOB3.-trueUltimateLOB4,0))


#consider CL decision trees
#each entry corresponds to an LDF Year
modelsWeightsPerLDF=Vector{Vector{Float64}}()
push!(modelsWeightsPerLDF,[80000,130000,160000,200000,350000,500000,15000000,15000000,15000000,15000000,15000000])
push!(modelsWeightsPerLDF,[80000,130000,160000,175000,300000,450000,750000,15000000,15000000,15000000,15000000])

for i=1:length(modelsWeightsPerLDF)
    @assert length(modelsWeightsPerLDF[i])==length(ayears)-1 # == number of LDFs
end

#array of LDFs
LDFArray=zeros(size(dataKnownByYE2005,1),length(ayears))
LDFArray[:,end].=1.0 #tail factor 1.0

#goal: build a tree for each LDF (i.e. a tree for 12-24 months, a second tree for 24-36 months, etc.)

treeResults=Dict{Float64,DataFrame}()
treeResultsAgg=Dict{Float64,DataFrame}()
error("redo this loop such that we have one minw per LDF. then consider the fit in Excel/reversals, etc.")
error("then consdier a boosting model?")

#actual model run (this may take a few minutes)
@time runModels!(dataKnownByYE2005,modelsWeightsPerLDF,treeResults,treeResultsAgg,LDFArray,selected_explanatory_vars,categoricalVars,folderForOutput,dtmtableKnownData);

#=
    resultingFiles,resM= runSingleModel(dataKnownByYE2005,160000,4,2005,selected_explanatory_vars,categoricalVars,folderForOutput)
=#

#define  struct to hold resulting data
struct ResStats
    comparisonByAY::DataFrame
    comparisonByLOB::Dict
    absErrorsByAY::Array 
    comparisons::Dict
end

#compare results (CL versus Tree versus Truth)
minws=sort(collect(keys(treeResultsAgg)))
modelStatistics=Dict{Float64,ResStats}()
for thiskey in minws
    @show thiskey
    resultTuple=customSummary(treeResultsAgg[thiskey],treeResults[thiskey])
    obj=ResStats(resultTuple[1],resultTuple[2],resultTuple[3],resultTuple[4])
    deltaTotalPerLOB=map(x->1/1e6.*(obj.comparisonByLOB[x][end,2].-obj.comparisonByLOB[x][end,3]),string.(collect(1:4)))
    @show deltaTotalPerLOB
    modelStatistics[thiskey]=obj
end

for kk in keys(modelStatistics)    
    @show kk
    obj=modelStatistics[kk]
    trueTotal=map(x->1/1e6.*(obj.comparisonByLOB[x][end,2]),string.(collect(1:4)))
    treeTotal=map(x->1/1e6.*(obj.comparisonByLOB[x][end,3]),string.(collect(1:4)))
    clTotal=map(x->1/1e6.*(obj.comparisonByLOB[x][end,4]),string.(collect(1:4)))
    deltaTotalPerLOB=map(x->1/1e6.*(obj.comparisonByLOB[x][end,2].-obj.comparisonByLOB[x][end,3]),string.(collect(1:4)))
    @show deltaTotalPerLOB,sum(deltaTotalPerLOB)
    @show clTotal,sum(clTotal)
    @show treeTotal,sum(treeTotal)
    @show trueTotal,sum(trueTotal)
end

@warn("If we were to update the 'tree-CL' factors such that they are based on all data (currently they are only using the training data), the model might further improve.")
@warn("make a note that 5m claims is a large data set!")