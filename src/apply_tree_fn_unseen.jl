export predict

#similar functions as in apply_tree_fn.jl
#however, for the functions in this file we do not require the featurepools to match
#i.e. DecisionTrees.assert_consistent_features could evaluate to false

function predict(x::Tree,f::DataFrame)
    asrt=DecisionTrees.assert_consistent_features(x.featurepools,f)
    if asrt
        return predictConsistentFeatures(x,f)
    else         
        @warn("DTM: Applying model even though model features do not match data features. Results might be inconsistent, you may want to review the predictions.")
        fp=x.featurepools
        res=apply_tree_by_leaf_inconsistent_features(x.rootnode,f,fp)
        return res 
    end
    #return apply_tree_by_leaf(x.rootnode,f)
end
    
function predict(x::BoostedTree,f::DataFrame;prroduceEstAndLeafMatrices::Bool=false)
    asrt=DecisionTrees.assert_consistent_features(x.featurepools,f)
    if asrt
        return predictConsistentFeatures(x::BoostedTree,f::DataFrame;prroduceEstAndLeafMatrices=prroduceEstAndLeafMatrices)
    else
        @show DecisionTrees.assert_consistent_features(x.featurepools,f)
        @error("in the works....this is not yet done for boosting")
    end
end

function apply_tree_by_leaf_inconsistent_features(t::Union{Leaf,Node{UInt8},Node{UInt16}},features::DataFrame,fp)
  nobs=size(features,1)
  fit=zeros(Float64,nobs)
  leafnrs=zeros(Int,nobs)
  idx=collect(1:nobs)
  apply_tree_by_leaf_inconsistent_features_iteration!(idx,t,features,fit,leafnrs,fp)  
  return fit,leafnrs
end


function apply_tree_by_leaf_inconsistent_features_iteration!(idx::Vector{Int},t::Node{T},features::DataFrame,fit::Vector{Float64},leaf::Vector{Int},fp) where T<:Unsigned
    if length(idx)==0
      return nothing
    end
   subset=t.subset 
   subsetWhichIsStringOrFloatVector=fp[t.featid_new_positive][subset]
   idxlInteger,idxrInteger=lrIndicesInconsitentFeatures(idx,features[t.featid_new_positive],subsetWhichIsStringOrFloatVector)  
   apply_tree_by_leaf_inconsistent_features_iteration!(idxlInteger,t.left,features,fit,leaf,fp)
   apply_tree_by_leaf_inconsistent_features_iteration!(idxrInteger,t.right,features,fit,leaf,fp)
   return nothing 
end

function apply_tree_by_leaf_inconsistent_features_iteration!(idx,t::Leaf,features::DataFrame,fit,leaf,fp)
     @inbounds for i in idx
       leaf[i]=t.id
       fit[i]=t.fitted
   end
   return nothing
end	 


"""
function lrIndicesInconsitentFeatures(trnidx::Vector{Int},f,subset::Array)
modifies trnidx to 'become' the left child index, returns r
l and r are defined by subset. l corresponds to the observations 'in' subset
this function is used for new/unseen features which have a different pool definition than the model data
"""
function lrIndicesInconsitentFeatures(trnidx::Vector{Int},f,subset::Array)
if eltype(f.pool)<:Number
    return lrIndicesForNumericalVariableInconsitentFeatures(trnidx,f,subset)
else
    return lrIndicesDefaultInconsitentFeatures(trnidx,f,subset)
end

end

"""
function lrIndicesDefaultInconsitentFeatures(trnidx::Vector{Int},f,subset::Array)
modifies trnidx to 'become' the left child index, returns r
l and r are defined by subset. l corresponds to the observations 'in' subset
this function is used for new/unseen features which have a different pool definition than the model data
"""
function lrIndicesDefaultInconsitentFeatures(trnidx::Vector{Int},f,subset::Array)
	l=Vector{Int}(undef,0)
	r=Vector{Int}(undef,0)
	sizehint!(l,length(trnidx))
	sizehint!(r,length(trnidx))
	for i in trnidx
		@inbounds thisref=f[i]
		if in(thisref,subset)
			push!(l,i)
		else
			push!(r,i)
		end
	end
	return l,r
end


"""
function lrIndicesForNumericalVariableInconsitentFeatures(trnidx::Vector{Int},f,subset::Array)
modifies trnidx to 'become' the left child index, returns r
l and r are defined by subset. l corresponds to the observations 'in' subset
this function is used for new/unseen features which have a different pool definition than the model data
"""
function lrIndicesForNumericalVariableInconsitentFeatures(trnidx::Vector{Int},f,subset::Array)
	l=Vector{Int}(undef,0)
	r=Vector{Int}(undef,0)
	sizehint!(l,length(trnidx))
    sizehint!(r,length(trnidx))
    subsetEnd=subset[end]
	for i in trnidx
		@inbounds thisref=f[i]                
        if thisref<=subsetEnd
            push!(l,i)
        else
            push!(r,i)
        end
	end
	return l,r
end
