
#define  struct to hold resulting data
struct ModelStats
    comparisonByAY::DataFrame
    comparisonByLOB::Dict
    absErrorsByAY::Array 
    comparisons::Dict
end

function buildTriangle(data;rowIdentifier=:AY,columns=[Symbol(string("PayCum",lpad(string(i),2,"0"))) for i=0:11])
    res=DataFrames.aggregate(data[vcat(columns,rowIdentifier)],rowIdentifier,sum);
    res2=convert(Array,res)
    return res2[:,2:end]    
end

function buildTriangleOLD(data;rowIdentifier=:AY,columns=[Symbol(string("PayCum",lpad(string(i),2,"0"))) for i=0:11])
    minRow,MaxRow=extrema(data[rowIdentifier])
    nRows=MaxRow-minRow+1
    nColumns=length(columns)
    res=zeros(Int,nRows,nColumns)
    have=names(data)
    @assert issubset(columns,have)
    for i=1:size(data,1)
        for j=1:length(columns)     
            @inbounds thisAY = data[rowIdentifier][i]
            @inbounds res[thisAY-minRow+1,j] += data[columns[j]][i]
        end  
    end
    return res
end

function chainLadder(x;tail=1.0)
    #x is the cumulative triangle
    k=1
    height=size(x,1)
    width=size(x,2) 
    LDFs=zeros(width)
    factorsToUltimate=ones(width)
    ultimate=zeros(height)
    paidToDate=zeros(height)
    for j=1:size(x,2)-1
        LDFs[j]=sum(x[1:height-j,j+1])/sum(x[1:height-j,j])        
    end
    LDFs[end]=tail
    factorsToUltimate[end]=LDFs[end]
    for j=size(x,2):-1:2
        factorsToUltimate[j-1] *= (factorsToUltimate[j] * LDFs[j-1])
        ultimate[height-j+2] = x[height-j+2,j-1] * factorsToUltimate[j-1]
        paidToDate[height-j+2] = x[height-j+2,j-1]
    end
    j=1+size(x,2)
    ultimate[height-j+2] = x[height-j+2,j-1] * factorsToUltimate[j-1]
    paidToDate[height-j+2] = x[height-j+2,j-1]
    reserves = ultimate .- paidToDate
    #construct output DataFrame
    reverse!(LDFs)
    reverse!(factorsToUltimate)
    df=DataFrame(LDFS=LDFs,factorsToUltimate=factorsToUltimate,paidToDate=paidToDate,reserves=reserves,ultimate=ultimate)
    return df #LDFs,factorsToUltimate,ultimate,reserves,paidToDate
end

function aggregateReservesOfIndividualCLModels(fullData,leafNrs)
    reservesCombined=zeros(length(unique(fullData[:AY])))
    uqLeafNrs=sort(unique(leafNrs))
    for leaf in uqLeafNrs
        partialData=fullData[leafNrs.==leaf,:]
        triangle=buildTriangle(partialData)
        resDF=chainLadder(triangle)
        reservesCombined .+= resDF[:reserves]
    end
    return reservesCombined
end

function calculateEstimates(fullData,LDFArray,aggCL,paidToDatePerRow)
    ayears=unique(fullData[:AY])
    #apply Chainladder to each row
    ultimatePerRow=zeros(size(fullData,1))
    #paidToDatePerRow=zeros(size(fullData,1))
    cumulativefactors=zeros(size(fullData,1))
    ultimatePerAY=zeros(length(ayears))
    truthPerAY=zeros(length(ayears))
    paycumcols=[Symbol(string("PayCum",lpad(string(i),2,"0"))) for i=0:11]
    minAY,maxAY=extrema(ayears)
    for iRow=1:size(fullData,1)
        @inbounds row=fullData[:AY][iRow]-minAY+1
        @inbounds rowInversed = length(ayears) - row + 1
        @inbounds cumulativefactor=prod(1 ./ LDFArray[iRow,rowInversed:end])
        @inbounds cumulativefactors[iRow] = cumulativefactor
        #@inbounds paidToDatePerRow[iRow] = fullData[paycumcols[rowInversed]][iRow]
        @inbounds ultimatePerRow[iRow] = cumulativefactor*paidToDatePerRow[iRow]        
        @inbounds ultimatePerAY[row] += ultimatePerRow[iRow]
        @inbounds truthPerAY[row] += fullData[:PayCum11][iRow]
    end    
    reservesPerAY=ultimatePerAY.-aggCL[:paidToDate]
    errPerAY=ultimatePerAY.-truthPerAY #aggCL[:ultimate]
    dfPerRow=DataFrame(paidToDate=paidToDatePerRow,reserves=ultimatePerRow.-paidToDatePerRow,ultimate=ultimatePerRow,cumulativefactors=cumulativefactors)
    df=DataFrame(paidToDate=deepcopy(aggCL[:paidToDate]),reserves=reservesPerAY,ultimate=ultimatePerAY,truth=truthPerAY,error=errPerAY)
    return dfPerRow,df 
end

function getPaidToDatePerRow(fullData)
    ayears=unique(fullData[:AY])
    paidToDatePerRow=zeros(size(fullData,1))    
    minAY,maxAY=extrema(ayears)
    paycumcols=[Symbol(string("PayCum",lpad(string(i),2,"0"))) for i=0:11]
    for iRow=1:size(fullData,1)
        @inbounds row=fullData[:AY][iRow]-minAY+1
        @inbounds rowInversed = length(ayears) - row + 1
        @inbounds paidToDatePerRow[iRow] = fullData[paycumcols[rowInversed]][iRow]
    end
    return paidToDatePerRow
end

function getcorrectedQuantilesOfError(estPerRow,truthPerRow,qtl_range,CLTotalUltimate)
    errPerRow=estPerRow[:ultimate]-truthPerRow
    totalUltimate=sum(estPerRow[:ultimate])    
    #what if we correct the ultimate total to match the CL total?
    corrFactor=1/totalUltimate*CLTotalUltimate
    #@show sum(corrFactor.*estPerRow[:ultimate])-CLTotalUltimate
    @assert isapprox(0,sum(corrFactor.*estPerRow[:ultimate])-CLTotalUltimate,atol=1e-6)
    errPerRowCorrected=corrFactor.*estPerRow[:ultimate]-truthPerRow    
    return Statistics.quantile(errPerRowCorrected,qtl_range)    
end

function checkIfPaymentsAreIncreasing(data)
    #consider payment columns
    paymCols=[Symbol(string("PayCum",lpad(ii,2,"0"))) for ii=0:11]
    isNonDecreasing=BitVector(undef,size(data,1))
    data2=data[paymCols]
    dataAsArray=transpose(convert(Array,data2))    
    for ii=1:size(data2,1)        
        @inbounds isNonDecreasing[ii]=issorted(dataAsArray[:,ii])
        #@inbounds isNonDecreasing[ii]=issorted(view(dataAsArray,:,ii))        
    end
    return isNonDecreasing
end



function customSummary(treeEstimateAgg::DataFrame,treeEstimatePerRow::DataFrame;writeResultToTemp=true)
    #NOTE: this function is using several global variables!

    #overall Result
    comparisonByAY=hcat(ayears,trueReserves,treeEstimateAgg[:reserves],clAllLOBs[:reserves])
    comparisonByAY=vcat(comparisonByAY,sum(comparisonByAY,dims=1))
    comparisonByAY[end,1]=9999
    comparisonByAY2=DataFrame(comparisonByAY)
    names!(comparisonByAY2,[:AY,Symbol("True Outstanding Amount"),Symbol("Estimated Reserves"),Symbol("CL Estimated Reserves")])
    
    absErrorsByAY=ceil.(Int,abs.(comparisonByAY[:,2].-comparisonByAY[:,3])/1000)
    @show absErrorsByAY
    @show absErrorsByAY[end]

    #result by LOB versus four CL  models (one per LOB)
    treeEstimatePerRow[:AY]=dataKnownByYE2005[:AY]
    comparisonByLOB=Dict{String,DataFrame}()
    for lob in sort(unique(dataKnownByYE2005[:LoB]))
        sumRes=DataFrames.aggregate(treeEstimatePerRow[dataKnownByYE2005[:LoB].==lob,[:reserves,:AY]],:AY,sum)
        #names!(sumRes,[:AY,:reserves])
        sumRes2=hcat(ayears,trueReservesPerLOB[lob][1:end-1],sumRes[:reserves_sum],clPerLOB[lob][:reserves][1:end-1])
        #sumRes[:CLReserves]=clPerLOB[lob][:reserves]
        #sumRes[:trueReserves]=trueReservesPerLOB[lob]
        sumRes2=vcat(sumRes2,sum(sumRes2,dims=1))
        sumRes2[end,1]=9999
        sumRes3=DataFrame(sumRes2)
        names!(sumRes3,[:AY,Symbol("True Outstanding Amount"),Symbol("Estimated Reserves"),Symbol("CL Estimated Reserves")])
        comparisonByLOB[lob]=deepcopy(sumRes3)
    end
    

    comparisons=Dict{String,DataFrame}()
    byVars=["cc","inj_part","age","LoB"]
    for vv in byVars
        #attach explanatory data to the tree estimates
        treeEstimatePerRow[Symbol(vv)]=dataKnownByYE2005[Symbol(vv)]
        #create comparison table
        comparison=DataFrames.aggregate(fullData[vcat(Symbol(vv),:PayCum11, :paidToDate)],Symbol(vv),sum)
        #rename table
        names!(comparison,vcat(Symbol(vv),:PayCum11, :paidToDate))
        #calculate reserves
        comparison[Symbol("True Outstanding Amount")]=comparison[:PayCum11].-comparison[:paidToDate]
        sort!(comparison,Symbol(vv))
        #calculate estimated reserves (tree model)
        est=DataFrames.aggregate(treeEstimatePerRow[vcat(Symbol(vv),:ultimate, :reserves,:paidToDate)],Symbol(vv),sum)
        sort!(est,Symbol(vv))
        #rename table
        names!(est,vcat(Symbol(vv),:ultimate,:reserves, :paidToDate))
        #check
        @assert all(isapprox.(est[:paidToDate].-comparison[:paidToDate],0))
        #attach estimated data (this should be a join (but as long as all values occurr in the data, we are fine))
        comparison[Symbol("Estimated Reserves")]=est[:reserves]
        #reorder columns
        comparison=comparison[[1,4,5,3,2]]
        #sort the table
        if !isa(comparison[1,Symbol(vv)],Number)
            comparison[:sort]=Meta.parse.(comparison[Symbol(vv)])
        else 
            comparison[:sort]=deepcopy(comparison[Symbol(vv)])
        end 
        sort!(comparison,:sort)
        #add variable name column
        comparison[:variable]=vv        
        comparison[:value]=string.(comparison[Symbol(vv)])
        delete!(comparison, Symbol(vv))
        #reorder columns
        comparison = comparison[[7,1,2,3,4,5,6]]
        
        comparisons[vv]=deepcopy(comparison)
        if writeResultToTemp
            if (isdir("C:\\temp\\"))
                CSV.write(string("C:\\temp\\",vv,".csv"),comparison)
            end
        end
        println(vv," : ",ceil(Int,sum(abs.(comparison[Symbol("True Outstanding Amount")].-comparison[Symbol("Estimated Reserves")]))/1000))
    end
    
    return comparisonByAY2,comparisonByLOB,absErrorsByAY,comparisons
    end
    
    
function runModels!(dataKnownByYE2005,dtmKnownByYE2005,modelsWeightsPerLDF,treeResults,treeResultsAgg,selected_explanatory_vars,categoricalVars,folderForOutput,clAllLOBs,paidToDatePerRow,settOrig::ModelSettings)

    kk=0
    Random.srand(1240)
    ayears=sort(unique(dataKnownByYE2005[:AY]))
    maxAY=maximum(ayears)
    
    #array of LDFs
    LDFArray=zeros(size(dataKnownByYE2005,1),length(ayears))
    LDFArray[:,end].=1.0 #tail factor 1.0
    LDFArraySmoohted=deepcopy(LDFArray)
    LDFArrayUnSmoohted=deepcopy(LDFArray)
        
    for minWeightPerLDF in modelsWeightsPerLDF
        fill!(LDFArray,0.0)
        LDFArray[:,end].=1.0 #tail factor 1.0
        kk+=1
        for ldfYear=1:length(ayears)-1 
            selectedWeight=minWeightPerLDF[ldfYear]  
            @show ldfYear
            @show selectedWeight
        
            resultingFiles,resM=runSingleModel(dataKnownByYE2005,dtmKnownByYE2005,selectedWeight,ldfYear,maxAY,selected_explanatory_vars,categoricalVars,folderForOutput,settOrig)
            #apply tree to the known claims
            if sett.model_type=="build_tree"
                fittedValues,leafNrs=predict(resM,dtmKnownByYE2005.features)  
            else
                dfpred=predict(resM,dtmKnownByYE2005.features)  
                fittedValues=dfpred[:RawEstimate]                                
                LDFArraySmoohted[:,ldfYear].=dfpred[:SmoothedEstimate]
                LDFArrayUnSmoohted[:,ldfYear].=dfpred[:UnsmoothedEstimate]
            end           
            #save estimated ldf per observation 
            LDFArray[:,ldfYear].=copy(fittedValues)
        end   

    #consider estimated ultimates per claim 
    estPerRow,estAgg=calculateEstimates(dataKnownByYE2005,LDFArray,clAllLOBs,paidToDatePerRow)    
    treeResults[kk]=deepcopy(estPerRow)
    treeResultsAgg[kk]=deepcopy(estAgg)
    #for boosting also store the other estimates
    if !(sett.model_type=="build_tree")
        estPerRow,estAgg=calculateEstimates(dataKnownByYE2005,LDFArraySmoohted,clAllLOBs,paidToDatePerRow)    
        treeResults[kk+0.2]=deepcopy(estPerRow)
        treeResultsAgg[kk+0.2]=deepcopy(estAgg)
        estPerRow,estAgg=calculateEstimates(dataKnownByYE2005,LDFArrayUnSmoohted,clAllLOBs,paidToDatePerRow)    
        treeResults[kk+0.4]=deepcopy(estPerRow)
        treeResultsAgg[kk+0.4]=deepcopy(estAgg)        
    end

    end 

return LDFArray
end 

function runSingleModel(dataKnownByYE2005,dtmKnownByYE2005,selectedWeight,ldfYear,maxAY,selected_explanatory_vars,categoricalVars,folderForOutput,settOrig::ModelSettings)    
    #a note on the data
    #dataKnownByYE2005 is a DataFrame (it contains more data such as all PayCum columns and the AY)
    #dtmKnownByYE2005 is a DTMTable object (it has a different structre than the DataFrame)    
    cumulativePaymentCol=Symbol(string("PayCum",lpad(ldfYear,2,"0")))
    cumulativePaymentColPrev=Symbol(string("PayCum",lpad(ldfYear-1,2,"0")))

    #subset data, conditions are: 
    #1) AY < 2005-ldfYear
    #2) positive denominator
    #3) cumPay is NOT zero for both ldfYear and ldfYear-1 (these rows do not add any value and can be discarded)    
    subsetIdx=(dataKnownByYE2005[:AY].<=(maxAY-ldfYear)) .& (dataKnownByYE2005[cumulativePaymentCol] .> 0) .& (.! ((dataKnownByYE2005[cumulativePaymentCol].==0) .& (dataKnownByYE2005[cumulativePaymentColPrev].==0)))    
    #subsetIdx is a boolean index, we need to convert it to an integer index for the next step
    subsetIdxInteger=findall(subsetIdx)
    #this step also creates a copy of the data (which is intended)
    dtmSubset=dtmKnownByYE2005[subsetIdxInteger]
    #redefine numerator and denominator (we are considering different PayCum columns of the data for each development year)
    dtmSubset.numerator=dataKnownByYE2005[cumulativePaymentColPrev][subsetIdxInteger]
    dtmSubset.denominator=dataKnownByYE2005[cumulativePaymentCol][subsetIdxInteger]    
    
    sett=deepcopy(settOrig)
    critSTR=string(sett.crit)[15:22]
    updateSettingsMod!(sett,minWeight=selectedWeight)
    resultingFiles,resM=dtm(dtmSubset,sett,file=joinpath(folderForOutput,string(sett.model_type[1:7],"_LDF_Year_",ldfYear,"crt_",critSTR,"_minw_",sett.minWeight,".txt")))    
    
    return resultingFiles,resM
end


function addTotal!(d::DataFrame)
    sm=[sum(d[i]) for i in 1 : size(d,2)]
    totRow=d[1,:]    
    nm=names(totRow)
    for j=1:size(d,2)
        totRow[j]=0.0
        if in(nm,[:LDFS,:factorsToUltimate])
        else 
            totRow[j]=sm[j]
        end
    end
    totRow[:LDFS]=1.0
    totRow[:factorsToUltimate]=totRow[:ultimate][end]/totRow[:paidToDate][end]
    append!(d,totRow)
    return nothing
end
