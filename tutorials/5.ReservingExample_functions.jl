
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
        @inbounds cumulativefactor=prod(1./LDFArray[iRow,rowInversed:end])
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
    return quantile(errPerRowCorrected,qtl_range)    
end

function checkIfPaymentsAreIncreasing(data)
    #consider payment columns
    paymCols=[Symbol(string("PayCum",lpad(ii,2,0))) for ii=0:11]
    isNonDecreasing=BitVector(undef,size(data,1))
    data2=data[paymCols]
    dataAsArray=transpose(convert(Array,data2))    
    for ii=1:size(data2,1)        
        @inbounds isNonDecreasing[ii]=issorted(dataAsArray[:,ii])
        #@inbounds isNonDecreasing[ii]=issorted(view(dataAsArray,:,ii))        
    end
    return isNonDecreasing
end



function customSummary(treeEstimateAgg::DataFrame,treeEstimatePerRow::DataFrame)
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
        sumRes2=hcat(ayears,trueReservesPerLOB[lob],sumRes[:reserves_sum],clPerLOB[lob][:reserves])
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
            sort!(comparison,:sort)
        end
        comparisons[vv]=deepcopy(comparison)
        #CSV.write(string("C:\\temp\\",vv,".csv"),comparison)
        println(vv," : ",ceil(Int,sum(abs.(comparison[Symbol("True Outstanding Amount")].-comparison[Symbol("Estimated Reserves")]))/1000))
    end
    
    return comparisonByAY2,comparisonByLOB,absErrorsByAY,comparisons
    end
    
    
function runModels!(dataKnownByYE2005,modelsWeightsPerLDF,treeResults,treeResultsAgg,LDFArray,selected_explanatory_vars,categoricalVars,folderForOutput,dtmtableKnownData)
    kk=0
    Random.srand(1240)
    ayears=sort(unique(dataKnownByYE2005[:AY]))
    maxAY=maximum(ayears)
    
for minWeightPerLDF in modelsWeightsPerLDF
    fill!(LDFArray,0.0)
    LDFArray[:,end].=1.0 #tail factor 1.0
    @show kk+=1    
    for ldfYear=1:length(ayears)-1 
        selectedWeight=minWeightPerLDF[ldfYear]  
        @show ldfYear
        @show selectedWeight
       
        resultingFiles,resM=runSingleModel(dataKnownByYE2005,selectedWeight,ldfYear,maxAY,selected_explanatory_vars,categoricalVars,folderForOutput)
        
        #=
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
        =#
            #apply tree to the known claims
            fittedValues,leafNrs=predict(resM,dtmtableKnownData.features)            
            #save estimated ldf per observation 
            LDFArray[:,ldfYear].=copy(fittedValues)
        end    
    end 

    #consider per Row Model 
    estPerRow,estAgg=calculateEstimates(dataKnownByYE2005,LDFArray,clAllLOBs,paidToDatePerRow)
    treeResults[kk]=deepcopy(estPerRow)
    treeResultsAgg[kk]=deepcopy(estAgg)
    
    #reservesPerAY=estAgg[:reserves]
    #errPerAY=estAgg[:error]
    #errTotal=sum(errPerAY)            
    #totalUltimate=sum(estPerRow[:ultimate])        
    return nothing
end 

function runSingleModel(dataKnownByYE2005,selectedWeight,ldfYear,maxAY,selected_explanatory_vars,categoricalVars,folderForOutput)
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
    return resultingFiles,resM
end