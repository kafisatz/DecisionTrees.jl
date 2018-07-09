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
@time triangleInclIBNyRClaims=buildTriangle(fullData)

#CL
#consider aggregate CL Model
triangle=buildTriangle(dataKnownByYE2005)
aggCL=chainLadder(triangle)
truthperAY=triangleInclIBNyRClaims[:,end]

#Consider CL MSE
ayperRow=dataKnownByYE2005[:AY]
CLfactorsToUltimate=aggCL[:factorsToUltimate]
CLLDFperRow=map(x->CLfactorsToUltimate[x-1993],ayperRow)
CLUltimatePerRow=CLLDFperRow.*paidToDatePerRow
@assert isapprox(0,sum(CLUltimatePerRow)-sum(aggCL[:ultimate]),atol=1e-4) #should be zero 
for i in ayears
    idx=ayperRow.==i
    @assert isapprox(0,sum(CLUltimatePerRow[idx])-aggCL[:ultimate][i-1993],atol=1e-4) #should be zero for each year
end
CLTotalUltimate=sum(CLUltimatePerRow)

minwList=[30000,40000,50000,75000,100000,200000]
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
for selectedWeight in minwList
    kk+=1    
    for ldfYear=1:length(ayears)-1        
        #subset the data
        thisdata=copy(dataKnownByYE2005)        
        #discard recent AYs (i.e ensure the data corresponds to a rectangle rather than a triangle)
        #filter!(row->row[:AY]<=maxAY-ldfYear, thisdata)
        thisdata=thisdata[thisdata[:AY].<=(maxAY-ldfYear),:]
        #Discard claims which had 0 cumulative payment by year 'ldfYear'
        #this is currently a requirement for the algorithms we apply (although one may be able to weaken this condition and the algorithms would still run)
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

        #run a single tree         
        if ldfYear>=LDF_Cutoff #for the tail of the triangle the models do not validate well, thus we want 'no tree' but a single leaf => same as CL factors
            actual_minimum_weight=-0.51
        else 
            actual_minimum_weight=copy(selectedWeight)
        end    
        updateSettingsMod!(sett,ignoreZeroDenominatorValues=true,minWeight=actual_minimum_weight,model_type="build_tree",writeTree=false)
        resultingFiles,resM=dtm(dtmtable,sett,file=joinpath(folderForOutput,string("minw_",sett.minWeight,"_ldfYear_",ldfYear,".txt")))

        #apply tree to the known claims
        fittedValues,leafNrs=predict(resM,dtmtableKnownData.features)

        #save estimated ldf per observation 
        LDFArray[:,ldfYear].=copy(fittedValues)
        if ldfYear>=LDF_Cutoff
            #@show fittedValues[1]
            #use Chain Ladder LDF factor
            LDFArray[:,ldfYear].= 1/reverse(aggCL[:LDFS])[ldfYear]
        end
    end 

    #consider per Row Model 
    estPerRow,estAgg=calculateEstimates(dataKnownByYE2005,LDFArray,aggCL,paidToDatePerRow)
    treeResults[selectedWeight]=deepcopy(estPerRow)
    treeResultsAgg[selectedWeight]=deepcopy(estAgg)
    
    reservesPerAY=estAgg[:reserves]
    errPerAY=estAgg[:error]
    errTotal=sum(errPerAY)        
    #errPerRow=estPerRow[:ultimate]-trueUltimatePerRow
    totalUltimate=sum(estPerRow[:ultimate])    
    #quantile(errPerRow,qtl_range)
    #median(errPerRow)
    #@show selectedWeight,mean(errPerRow.^2)
    #treeMSE[kk]=mean(errPerRow.^2)
    
    #what if we correct the ultimate total to match the CL total?
    #corrFactor=1/totalUltimate*CLTotalUltimate    
    #@assert isapprox(0,sum(corrFactor.*estPerRow[:ultimate])-CLTotalUltimate,atol=1e-6)    
end 


#Consider the quantiles of the error
qtls=map(x->quantile(treeResults[x][:ultimate].-trueUltimatePerRow,qtl_range),minwList)
expectedQtls35k=  [-1.4657e5, -28191.4, -882.096, 0.0, 0.0, 0.0, 0.0, 9.02446, 53.2051, 78.8788, 221.046, 513.365, 2300.84, 12311.4, 44644.0][:]
@assert all(isapprox.(qtls[minwList.==35000][1].-expectedQtls35k,0,atol=1)) 

#Consider the quantiles of the error when we correct the total ultimate such that it matches the CL Ultimate
#qtlsWithCorrectedTotalUltimate=map(x->getcorrectedQuantilesOfError(treeResults[x],trueUltimatePerRow,qtl_range,CLTotalUltimate),minwList)

function summaryByAY(fullData,trueUltimatePerRow,trueReservePerRow,paidToDatePerRow,ayears,treeEstimate::DataFrame,clEstimatedUltimate::Vector;by="LoB")
    @assert in(by,selected_explanatory_vars)
    byS=Symbol(by)
    uniqueValues=sort(unique(fullData[byS]))
    ayShift=minimum(ayears)-1
    
    ayColumn=fullData[:AY]
    treeReserve=treeEstimate[:reserves]    
    clReserve=clEstimatedUltimate.-paidToDatePerRow

    res=Dict{eltype(uniqueValues),Array{Float64,2}}()

    for thisValue in uniqueValues
        res0=zeros(length(ayears)+1,4)        
        res0[1:end-1,1]=ayears        
        #note: this is very inefficiently programmed! (a countsort would be preferable)
        for i=1:size(fullData,1)
            @inbounds rowMatchesValue=fullData[byS][i]==thisValue
            if rowMatchesValue                                    
                @inbounds row=ayColumn[i]-ayShift
                @inbounds res0[row,2]+=trueReservePerRow[i]
                @inbounds res0[row,3]+=treeReserve[i]
                @inbounds res0[row,4]+=clReserve[i]                    
            end            
        end        
        resTot=sum(res0,1)
        #consider total 
        res0[end,1]=9999
        res0[end,2:end].=resTot[2:end]
        res[thisValue]=deepcopy(res0)
    end 
    
    return res
end

selectedTree=treeResults[35000]
selectedTree[:ultimate]

trueReserves=sum(trueReservePerRow)
trueReserves/1000
CLReserves=sum(CLUltimatePerRow.-paidToDatePerRow)
treeReserves=sum(selectedTree[:reserves])
(treeReserves-trueReserves)/1000
(CLReserves-trueReserves)/1000

268338+208376+214768+424738

r1=summaryByAY(fullData,trueUltimatePerRow,trueReservePerRow,paidToDatePerRow,ayears,treeResults[35000],CLUltimatePerRow,by="LoB")


@warn("If we were to update the 'tree-CL' factors such that they are based on all data (currently they are only using the training data), the model might further improve.")
@warn("make a note that 5m claims is a large data set!")