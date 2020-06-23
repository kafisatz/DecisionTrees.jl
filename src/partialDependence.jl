export partialDependence 

function partialDependence(mdl::Tree, dtmtable;thisIdx::Vector{Int}=dtmtable.trnidx,listOfFeatures::Vector{Symbol}=propertynames(dtmtable.features))
    @assert issubset(listOfFeatures, propertynames(dtmtable.features)) "DTM: At least one of the feature names you provided were not found in the names of the feature dataframe\r\nYou provided: $(listOfFeatures)\rpropertynames(dtmtable.features)=$(propertynames(dtmtable.features))"

    copyOfFeatures = deepcopy(dtmtable.features[thisIdx,:])

    nobs = size(copyOfFeatures, 1)
    fitReused = zeros(Float64, nobs)
    leafnrsReused = zeros(Int, nobs)
    idx = collect(1:nobs)
    # listOfFeatures=[:Density,:BonusMalus]

    res = Dict{Symbol,DataFrames.DataFrame}()
    for thisFeature in listOfFeatures
        forig = deepcopy(copyOfFeatures[thisFeature])
        pool = copyOfFeatures[thisFeature].pool
        predictedValues = zeros(Float64, length(pool))
        for i = 1:length(pool)
            copyOfFeatures[thisFeature] .= pool[i]        
            DecisionTrees.apply_tree_by_leaf_iteration!(idx, mdl.rootnode, copyOfFeatures, fitReused, leafnrsReused) 
            predictedValues[i] = Statistics.mean(fitReused)
        end
        # save result in dict
        df1 = DataFrame(Values=deepcopy(pool), AveragePredictions=copy(predictedValues))
        res[thisFeature] = deepcopy(df1)
        # reset values to original values
        copyOfFeatures[thisFeature] = deepcopy(forig)  
    end
    return res
end

# this function is currently very slow....
function partialDependence(mdl::BoostedTree, dtmtable;fitSelection=:RawEstimate,thisIdx::Vector{Int}=dtmtable.trnidx,listOfFeatures::Vector{Symbol}=propertynames(dtmtable.features))
    @assert issubset(listOfFeatures, propertynames(dtmtable.features)) "DTM: At least one of the feature names you provided were not found in the names of the feature dataframe\r\nYou provided: $(listOfFeatures)\rpropertynames(dtmtable.features)=$(propertynames(dtmtable.features))"

    copyOfFeatures = deepcopy(dtmtable.features[thisIdx,:])
   
    nobs = size(copyOfFeatures, 1)
    # fitReused=zeros(Float64,nobs)
    # leafnrsReused=zeros(Int,nobs)
    idx = collect(1:nobs)
    # listOfFeatures=[:Density,:BonusMalus]

    res = Dict{Symbol,DataFrames.DataFrame}()
    for thisFeature in listOfFeatures
        forig = deepcopy(copyOfFeatures[thisFeature])
        pool = copyOfFeatures[thisFeature].pool
        predictedValues = zeros(Float64, length(pool))
        for i = 1:length(pool)
            copyOfFeatures[thisFeature] .= pool[i]        
            result = predict(mdl, copyOfFeatures)            
            predictedValues[i] = Statistics.mean(result[fitSelection])  
        end
        # save result in dict
        df1 = DataFrame(Values=deepcopy(pool), AveragePredictions=copy(predictedValues))
        res[thisFeature] = deepcopy(df1)
        # reset values to original values        
        copyOfFeatures[thisFeature] = deepcopy(forig)  
    end
    return res
end

# partialDependence(resM,dtmtable,listOfFeatures=[:Density,:VehGas])