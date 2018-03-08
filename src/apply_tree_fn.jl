export predict

#Bagging
function apply_tree_by_leaf(idx::Vector{Int},leaves_of_tree::Array{Leaf,1},t::Tree,features::DataFrame)	
	return apply_tree_by_leaf(idx,leaves_of_tree,t.rootnode,features)
end

################################################################################################################
# these functions are used internally (e.g. by the boosting algorithm)  
################################################################################################################

"""
This is the core apply_tree_by_leaf! function which is used for boosting (and other models)
apply_tree_by_leaf!(fit::Vector{Float64},leaf::Vector{Int},idx::Vector{Int},t::Union{Leaf,Node{UInt8},Node{UInt16}},features::DataFrame)  
"""
function apply_tree_by_leaf!(fit::Vector{Float64},leaf::Vector{Int},idx::Vector{Int},t::Union{Leaf,Node{UInt8},Node{UInt16}},features::DataFrame)  
  nobs=size(features,1) 
  @assert nobs==length(fit)==length(leaf)
  #rpvector=[x.rule_path for x in leaves_of_tree] #rule_path of each leaf is used to identify a leaf in the following (this could probably be done more efficiently) todo/tbd
  apply_tree_by_leaf_iteration!(idx,t,features,fit,leaf)  
return nothing
end
  
function apply_tree_by_leaf(idx::Vector{Int},t::Union{Leaf,Node{UInt8},Node{UInt16}},features::DataFrame)  
  #as we always supply an index here, the function will always return two arrays of length size(features,2)
  #thus the "training" observations will be empty if this function is called on a validation index only!!
  nobs=size(features,1) 
  #rpvector=[x.rule_path for x in leaves_of_tree] #rule_path of each leaf is used to identify a leaf in the following (this could probably be done more efficiently) todo/tbd
  fit=zeros(Float64,nobs)
  leaf=zeros(Int,nobs)
  apply_tree_by_leaf_iteration!(idx,t,features,fit,leaf)  
  return fit,leaf
end

function apply_tree_by_leaf_iteration!(idx::Vector{Int},t::Node{T},features::DataFrame,fit::Vector{Float64},leaf::Vector{Int}) where T<:Unsigned
	 if length(idx)==0
       #info("DTM: No data (an empty index) was provided to apply tree function")
       return nothing
	 end
	 featid=t.featid
   subset=t.subset   
   idxlInteger,idxrInteger=lrIndices(idx,features[t.featid_new_positive],subset)   
       #here each iteration will write on certain elements of the two vectors fit and leaf
     #@code_warntype apply_tree_by_leaf_iteration!(idxlInteger,rpvector,t.left,features,fit,leaf)
      apply_tree_by_leaf_iteration!(idxlInteger,t.left,features,fit,leaf)
      apply_tree_by_leaf_iteration!(idxrInteger,t.right,features,fit,leaf)
    return nothing 
 end
 
function apply_tree_by_leaf_iteration!(idx,t::Leaf,features::DataFrame,fit,leaf)
	nobs=length(idx)
	#numm=firstmatch(rpvector,t.rule_path) #todo/tbd: can this be done in a better way? firstmatch is sort of a hack since I could not get findin working on the first try.... the match should be unique by definition
    @inbounds for i in idx
        leaf[i]=t.id
        fit[i]=t.fitted
    end
    return nothing
end	 

################################################################################################################
# functions for predicting estimates on unseen data 
################################################################################################################

"""
Apply tree by leaf function for 'new' data
"""
function apply_tree_by_leaf(t::Union{Leaf,Node{UInt8},Node{UInt16}},features::DataFrame)
  #fit::Vector{Float64},leaf::Vector{Int},
  nobs=size(features,1)
  fit=zeros(Float64,nobs)
  leafnrs=zeros(Int64,nobs)
  idx=collect(1:nobs)
  apply_tree_by_leaf_iteration!(idx,t,features,fit,leafnrs)  
  return fit,leafnrs
end


function predict(x::Tree,f::DataFrame)
  #warn("need to ensure, data is consistent with the data used to derive the tree!")
  @assert assert_consistent_features(x.featurepools,f)
  return apply_tree_by_leaf(x.rootnode,f)
end

function predict(x::BoostedTree,f::DataFrame,boolProduceEstAndLeafMatrices::Bool=false)
  #warn("need to ensure, data is consistent with the data used to derive the tree!")
  @assert assert_consistent_features(x.featurepools,f)
  
  obs=size(f,1)
  sett=x.settings
  trn_meanobservedvalue=x.meanobserved
  niter=sett.niter
  
  estimatedRatio=ones(obs).*trn_meanobservedvalue  
  currentRelativity=ones(obs)
  indicatedRelativityForApplyTree_reused=zeros(obs)
  reused_fitted_leafnr_vector=zeros(Int64,obs)
  est_SmoothedEstFromScores=zeros(obs)
  est_UnSmoothedEstFromScores=zeros(obs)
  
  thisidx=collect(1:obs)
  empty_validx=ones(Int64,0)

  if boolProduceEstAndLeafMatrices
		est_matrix=Array{Float64}(obs,iterations+1)
		est_matrixFromScores=deepcopy(est_matrix)
		MatrixOfLeafNumbers=Array{Int}(obs,iterations+1)
		MatrixOfLeafNumbers[:,1]=0
		est_matrix[:,1]=copy(transpose(estimatedRatio))
		est_matrixFromScores[:,1]=copy(transpose(estimatedRatio))		
	else
		est_matrix=Array{Float64}(0,0)
		est_matrixFromScores=deepcopy(est_matrix)
		MatrixOfLeafNumbers=Array{Int}(0,0)
	end
	
  for iter=1:niter    
    current_mdf=x.moderationvector[iter]
    
    #estimatedNumerator=estimatedRatio.*denominator #these are for instance the estimated losses for a LR model
    #vectorOfLeafArrays[iter+1]=create_leaves_array(res[iter])

    apply_tree_by_leaf!(indicatedRelativityForApplyTree_reused,reused_fitted_leafnr_vector,thisidx,x.trees[iter],f)    
    #Moderate the estimate
    _moderate!(estimatedRatio,indicatedRelativityForApplyTree_reused,current_mdf)
    update_current_rels!(currentRelativity,estimatedRatio,trn_meanobservedvalue)

    if boolProduceEstAndLeafMatrices
      write_column!(est_matrix,iter+1,estimatedRatio)
      MatrixOfLeafNumbers[:,iter+1]=reused_fitted_leafnr_vector	
    end
  end

  #derive scores from raw rels
  scores=map(z->searchsortedfirst(x.maxRawRelativityPerScoreSorted,z),currentRelativity)
  #derive estimates
  update_mapped_estFromScores!(est_UnSmoothedEstFromScores,x.rawObservedRatioPerScore,scores)
  update_mapped_estFromScores!(est_SmoothedEstFromScores,x.ScoreToSmoothedEstimate,scores)

  #NOTE: we have three different estimates
  est_RawRelativities=currentRelativity.*trn_meanobservedvalue
  est_UnSmoothedEstFromScores
  est_SmoothedEstFromScores

  result=DataFrame()
  result[:scores]=scores
  result[:SmoothedEstimate]=est_SmoothedEstFromScores
  result[:UnsmoothedEstimate]=est_UnSmoothedEstFromScores
  result[:RawEstimate]=est_RawRelativities
  
  return result
end

function assert_consistent_features(fp,f::DataFrame)
  @assert length(fp)==size(f,2)
  #the above condition could also be weakend (as could the loop below) we only need to check the variables that were actually used by the model. 
  for j=1:length(fp)
    if !(fp[j]==f[j].pool)
      @show !(fp[j]==f[j].pool)
      @show issubset(fp[j],f[j].pool)
      error("Features do not match model features. Note: this condition could be weakend: the features of the data which is provided only needs to be a subset of the pools used during modelling")      
      error("Abort.")
    end
  end
  return true
end