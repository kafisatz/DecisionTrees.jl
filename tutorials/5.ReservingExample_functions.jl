
function buildTriangle(data;rowIdentifier=:AY,columns=[Symbol(string("PayCum",lpad(i,2,0))) for i=0:11])
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

function calculateEstimates(fullData,LDFArray,triangle,aggCL)
    ayears=unique(fullData[:AY])
    #apply Chainladder to each row
    ultimatePerRow=zeros(size(fullData,1))
    paidToDatePerRow=zeros(size(fullData,1))
    cumulativefactors=zeros(size(fullData,1))
    ultimatePerAY=zeros(length(ayears))
    truthPerAY=zeros(length(ayears))
    paycumcols=[Symbol(string("PayCum",lpad(i,2,0))) for i=0:11]
    minAY,maxAY=extrema(ayears)
    for iRow=1:size(fullData,1)
        @inbounds row=fullData[:AY][iRow]-minAY+1
        @inbounds rowInversed = length(ayears) - row + 1
        @inbounds cumulativefactor=prod(1./LDFArray[iRow,rowInversed:end])
        @inbounds cumulativefactors[iRow] = cumulativefactor
        @inbounds paidToDatePerRow[iRow] = fullData[paycumcols[rowInversed]][iRow]
        @inbounds ultimatePerRow[iRow] = cumulativefactor*paidToDatePerRow[iRow]        
        @inbounds ultimatePerAY[rowInversed] += ultimatePerRow[iRow]
        @inbounds truthPerAY[row] += fullData[:PayCum11][iRow]
    end
    #reserves=fullData[:PayCum11]
    ultimate
    reservesPerAY=ultimatePerAY.-aggCL[:paidToDate]
    errPerAY=ultimatePerAY.-truthPerAY #aggCL[:ultimate]
    dfPerRow=DataFrame(paidToDate=paidToDatePerRow,reserves=ultimatePerRow.-paidToDatePerRow,ultimate=ultimatePerRow,cumulativefactors=cumulativefactors)
    df=DataFrame(paidToDate=deepcopy(aggCL[:paidToDate]),reserves=reservesPerAY,ultimate=ultimatePerAY,truth=truthPerAY,error=errPerAY)
    return dfPerRow,df 
    end