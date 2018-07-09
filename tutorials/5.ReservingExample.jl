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

#Define a random split into training and validaiton data 
#training shall be 70% of the data
Random.srand(1239)
trnValCol=BitArray(rand() < 0.7 for x in 1:size(fullData,1));
fullData[:trnValCol]=trnValCol

@warn("We may want to redo the whole calculation with no validation data -> the CL factors will be based on all of the data (not only the training porportion).")

#The data should match this total for the column Paycum11
@assert sum(fullData[:PayCum11])==1477250848 "Expected 1464014219, have: $(sum(fullData[:PayCum11]))"
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
    println(names(fullData)[x])
    @show size(unique(fullData[:,x]))
end

#Let us see whether all claims are 'increasing' in time (i.e. whether all incremental payments are non-negative)
function checkIfPaymentsAreIncreasing(data)
    #consider payment columns
    paymCols=[Symbol(string("PayCum",lpad(ii,2,0))) for ii=0:11]
    isNonDecreasing=BitVector(undef,size(data,1))
    @time for ii=1:length(isNonDecreasing)
        @inbounds thisrow=vec(convert(Array,view(data,ii,paymCols)))
        @inbounds isNonDecreasing[ii]=issorted(thisrow)
    end
    return isNonDecreasing
end

isNonDecreasing=checkIfPaymentsAreIncreasing(fullData) #takes about 18 seconds
sum(isNonDecreasing)/length(isNonDecreasing) # 99.67% of the claims are nondecreasing
nonDecreasingRows=.!isNonDecreasing
nonDecreasingRowsInteger=find(nonDecreasingRows)

#claims which are decreasing over time (at some point)
fullData[nonDecreasingRowsInteger,:]

##############################
#Prepare the data
##############################
#consider the Accient Years
ayears=unique(fullData[:AY])
maxAY=maximum(ayears)
#For the first model we perform we consider the Cumulative Payment after 1 year and compare it to the Cumulative Payment after 2 years

#array of LDFs
LDFArray=zeros(size(fullData,1),length(ayears))
LDFArray[:,end].=1.0 #tail factor 1.0
#prepare the whole dataset in the same manner as the training data set
dtmtablefullData,settfullData,dfpreppedfullData=prepare_dataframe_for_dtm!(fullData,keycol="ClNr",numcol="PayCum00",trnvalcol="trnValCol",independent_vars=selected_explanatory_vars,treat_as_categorical_variable=["LoB","inj_part","cc"]);

#consider paidToDate and 'truth'
paidToDatePerRow=getPaidToDatePerRow(fullData)
truthPerRow=fullData[:PayCum11]

#CL
#consider aggregate CL Model
triangle=buildTriangle(fullData)
aggCL=chainLadder(triangle)
truthperAY=triangle[:,end]

#Consider CL MSE
ayperRow=fullData[:AY]
CLfactorsToUltimate=aggCL[:factorsToUltimate]
CLLDFperRow=map(x->CLfactorsToUltimate[x-1993],ayperRow)
CLUltimatePerRow=CLLDFperRow.*paidToDatePerRow
@assert isapprox(0,sum(CLUltimatePerRow)-sum(aggCL[:ultimate])) #should be zero 
for i in ayears
    idx=ayperRow.==i
    @assert isapprox(0,sum(CLUltimatePerRow[idx])-aggCL[:ultimate][i-1993],atol=1e-7) #should be zero for each year
end
CLerrPerRow=CLUltimatePerRow-truthPerRow
CLTotalUltimate=sum(CLUltimatePerRow)
qtl_range=[.0001,.001,.01,.05,.1,.2,.25,.5,.75,.8,.9,.95,.99,.999,.9999] 
CLqtls=quantile(CLerrPerRow,qtl_range)
median(CLerrPerRow)
CLmse=mean(CLerrPerRow.^2)
#CL LDFs 1.624586303	1.143772133	1.056608075	1.029755909	1.019553623	1.014561794	1.011186622	1.008085459	1.007157648	1.004935975	1.003334179

minwList=[35000]
LDF_Cutoff=4 #for the years >=LDF_Cutoff, we use the CL factors

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
        thisdata=copy(fullData)
        #discard recent AYs (i.e ensure the data corresponds to a rectangle rather than a triangle)
        thisdata=thisdata[thisdata[:AY].<=maxAY-ldfYear,:]
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
        #NOTE for now we are discarding all claims which have zero payment in the more recent year which dermines the CL factor
        #we can keep the claims with zero values for cumulativePaymentCol
        #however, the algorithm will fail when a segment with only such claims is identified (which should be unlikely if sett.minWeight is large enough) 
        #thisdata=thisdata[thisdata[cumulativePaymentCol].>0,:]
        dtmtable,sett,dfprepped=prepare_dataframe_for_dtm!(thisdata,trnvalcol="trnValCol",keycol="ClNr",numcol=string(cumulativePaymentColPrev),denomcol=string(cumulativePaymentCol),independent_vars=selected_explanatory_vars,treat_as_categorical_variable=["LoB","inj_part","cc"]);

        #run a single tree         
        if ldfYear>=LDF_Cutoff #for the tail of the triangle the models do not validate well, thus we want 'no tree' but a single leaf => same as CL factors
            actual_minimum_weight=-0.51
        else 
            actual_minimum_weight=copy(selectedWeight)
        end    
        updateSettingsMod!(sett,ignoreZeroDenominatorValues=true,minWeight=actual_minimum_weight,model_type="build_tree",writeTree=false)
        resultingFiles,resM=dtm(dtmtable,sett,file=joinpath(folderForOutput,string("minw_",sett.minWeight,"_ldfYear_",ldfYear,".txt")))

        #apply tree to the FULL data set
        fittedValues,leafNrs=predict(resM,dtmtablefullData.features)

        #save estimated ldf per observation 
        LDFArray[:,ldfYear].=copy(fittedValues)
        if ldfYear>=LDF_Cutoff
            #@show fittedValues[1]
            #use Chain Ladder LDF factor
            LDFArray[:,ldfYear].= 1/reverse(aggCL[:LDFS])[ldfYear]
        end
    end 

    #consider per Row Model 
    estPerRow,estAgg=calculateEstimates(fullData,LDFArray,triangle,aggCL,paidToDatePerRow)
    treeResults[selectedWeight]=deepcopy(estPerRow)
    treeResultsAgg[selectedWeight]=deepcopy(estAgg)
    
    reservesPerAY=estAgg[:reserves]
    errPerAY=estAgg[:error]
    errTotal=sum(errPerAY)        
    errPerRow=estPerRow[:ultimate]-truthPerRow
    totalUltimate=sum(estPerRow[:ultimate])    
    #quantile(errPerRow,qtl_range)
    #median(errPerRow)
    #@show selectedWeight,mean(errPerRow.^2)
    treeMSE[kk]=mean(errPerRow.^2)
    
    #what if we correct the ultimate total to match the CL total?
    corrFactor=1/totalUltimate*CLTotalUltimate    
    @assert isapprox(0,sum(corrFactor.*estPerRow[:ultimate])-CLTotalUltimate,atol=1e-6)
    errPerRowCorrected=corrFactor.*estPerRow[:ultimate]-truthPerRow    
    quantile(errPerRowCorrected,qtl_range)    
end 


#Consider the quantiles of the error
qtls=map(x->quantile(treeResults[x][:ultimate].-truthPerRow,qtl_range),minwList)
expectedQtls35k= [-3.00292e5, -42808.2, -966.453, 0.0, 0.0, 0.0, 0.0, 4.32545, 45.9819, 74.31, 247.741, 647.852, 3773.03, 23917.0, 90376.7][:]
@assert all(isapprox.(qtls[minwList.==35000][1].-expectedQtls35k,0,atol=1)) 

error("mi")

#Consider the quantiles of the error when we correct the total ultimate such that it matches the CL Ultimate
qtlsWithCorrectedTotalUltimate=map(x->getcorrectedQuantilesOfError(treeResults[x],truthPerRow,qtl_range,CLTotalUltimate),minwList)

@show CLqtls

qtls[1]./CLqtls
qtls[2]./CLqtls

treeMSE
CLmse
treeMSE./CLmse
#this is the same as treeMSE:
mean.(map(x->abs2.(treeResults[x][:ultimate]-truthPerRow),minwList))

#we can correct the aggregate total in order for it to match the CL ultimate
@warn("correct the ultimate estimate?!? with a certain factor")

#the chain ladder volume weighted LDF for Y1-Y2 is 1.62458536172099 
#its inverse is 0.615541678241308 
#this is the average observed value of the root node (if we consider the validation and training data together)
#the tree has identified several segments (i.e. leaves) where this ratio is either higher or lower than that value

#consider the fitted values and the leaf numbers

#fit model 
#=
updateSettingsMod!(sett,minWeight=-0.2,model_type="build_tree")
resultingFiles,resM=dtm(dtmtable,sett)
#predict the fittedValues and Leaf Numbers on the total data
fittedValues,leafNrs=predict(resM,dtmtablefullData.features)
reservesCombined=aggregateReservesOfIndividualCLModels(fullData,leafNrs)

DelimitedFiles.writedlm("C:\\temp\\1.csv",reservesCombined,',')
=#

@warn("show scores y1/y2 ldf versus y2/y3 ldf in an xy scatter plot")
