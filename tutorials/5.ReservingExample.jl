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
Distributed.@everywhere import DataFrames: DataFrame,groupby,combine,names!

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
ayears=unique(fullData[:AY])
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
minwListForEachYear=[75000,125000,150000,200000,350000,500000,15000000,15000000,15000000,15000000,15000000,15000000]

minwList=[40000,50000,75000,100000,200000]
LDF_Cutoff=5 #for the years >=LDF_Cutoff, we use the CL factors

#array of LDFs
LDFArray=zeros(size(dataKnownByYE2005,1),length(ayears))
LDFArray[:,end].=1.0 #tail factor 1.0

#goal: build a tree for each LDF (i.e. a tree for 12-24 months, a second tree for 24-36 months, etc.)

#selectedWeight=-0.15
treeMSE=zeros(length(minwList))
treeResults=Dict{Float64,DataFrame}()
treeResultsAgg=Dict{Float64,DataFrame}()
kk=0
Random.srand(1240)
error("redo this loop such that we have one minw per LDF. then consider the fit in Excel/reversals, etc.")
error("then consdier a boosting model?")
for selectedWeight in minwList
    kk+=1    
    for ldfYear=1:length(ayears)-1        
        if ldfYear>=LDF_Cutoff #for the tail of the triangle the models do not validate well, thus we want 'no tree' but a single leaf => same as CL factors
            #actual_minimum_weight=-0.51
            #use Chain Ladder LDF factor
            LDFArray[:,ldfYear].= 1/reverse(clAllLOBs[:LDFS])[ldfYear]
        else 
            #subset the data
            thisdata=copy(dataKnownByYE2005)        
            #discard recent AYs (i.e ensure the data corresponds to a rectangle rather than a triangle)        
            thisdata=thisdata[thisdata[:AY].<=(maxAY-ldfYear),:]
            
            cumulativePaymentCol=Symbol(string("PayCum",lpad(ldfYear,2,0)))
            cumulativePaymentColPrev=Symbol(string("PayCum",lpad(ldfYear-1,2,0)))

            #discard all claims which are zero for 'both columns'
            discardIdx=(thisdata[cumulativePaymentCol].==0) .& (thisdata[cumulativePaymentColPrev] .==0)        
            if any(discardIdx)
                thisdata=thisdata[.!discardIdx,:]
                println("Discarded $(sum(discardIdx)) observations with PayCum zero for both development years $(ldfYear) and $(ldfYear-1)")
            end 
            #Note: in a previous version of the model, we were discarding all claims which have zero payment in the more recent year which determines the CL factor
            #It turns out that we can actaully keep these claims we can keep the claims 
            #However, the algorithm will fail when a segment with only such claims is identified (which should be unlikely if sett.minWeight is large enough) 
            #thisdata=thisdata[thisdata[cumulativePaymentCol].>0,:]
            dtmtable,sett,dfprepped=prepare_dataframe_for_dtm!(thisdata,trnvalcol="trnValCol",keycol="ClNr",numcol=string(cumulativePaymentColPrev),denomcol=string(cumulativePaymentCol),independent_vars=selected_explanatory_vars,treat_as_categorical_variable=categoricalVars);        

            #run a single tree model
            updateSettingsMod!(sett,ignoreZeroDenominatorValues=true,minWeight=selectedWeight,model_type="build_tree",write_dot_graph=true,writeTree=false,graphvizexecutable="C:\\Program Files (x86)\\Graphviz2.38\\bin\\dot.exe")
            resultingFiles,resM=dtm(dtmtable,sett,file=joinpath(folderForOutput,string("minw_",sett.minWeight,"_ldfYear_",ldfYear,".txt")))
            #apply tree to the known claims
            fittedValues,leafNrs=predict(resM,dtmtableKnownData.features)
            #save estimated ldf per observation 
            LDFArray[:,ldfYear].=copy(fittedValues)
        end    
    end 

    #consider per Row Model 
    estPerRow,estAgg=calculateEstimates(dataKnownByYE2005,LDFArray,clAllLOBs,paidToDatePerRow)
    treeResults[selectedWeight]=deepcopy(estPerRow)
    treeResultsAgg[selectedWeight]=deepcopy(estAgg)
    
    reservesPerAY=estAgg[:reserves]
    errPerAY=estAgg[:error]
    errTotal=sum(errPerAY)            
    totalUltimate=sum(estPerRow[:ultimate])        
end 

#create struct to hold resulting data
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