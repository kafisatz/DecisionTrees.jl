function apply_tree_by_row(leaves_of_tree::Vector{Leaf},tree::Tree,numfeatures::Array{pdaMod,1},charfeatures::Array{pdaMod,1})
  t=tree.rootnode
  length(numfeatures)>0 ? nobs=length(numfeatures[1]) : nobs=length(charfeatures[1])
  @assert (size(numfeatures,1)==0 || length(numfeatures[1])==nobs)
  @assert (size(charfeatures,1)==0 || length(charfeatures[1])==nobs)
  #Julia’s pmap() is designed for the case where each function call does a large amount of work.
  #In contrast, @parallel for can handle situations where each iteration is tiny, perhaps merely summing two numbers.
  #Only worker processes are used by both pmap() and @parallel for for the parallel computation.
  #In case of @parallel for, the final reduction is done on the calling process
  res=Array{Float64}(nobs)
  leafNumberVal=Array{Int64}(nobs)
  rpvector=[x.rule_path for x in leaves_of_tree] #rule_path of each leaf is used to identify a leaf in the following (this could probably be done more efficiently) todo/tbd
  for row=1:nobs
    res[row],leafNumberVal[row]=apply_tree_by_row(rpvector,t,numfeatures,charfeatures,row)
  end
  return res,leafNumberVal
end

function apply_tree_by_row(rpvector::Array{Array{Rulepath,1},1},t::Node,numfeatures::Array{pdaMod,1},charfeatures::Array{pdaMod,1},row::Int64)
#function for a single observation only
  if t.featid<0
    #character split
    if charfeatures[-t.featid].pda.refs[row] in t.subset
      return apply_tree_by_row(rpvector,t.left,numfeatures,charfeatures,row)
    else
      return apply_tree_by_row(rpvector,t.right,numfeatures,charfeatures,row)
    end
  else
    #numeric split
    if numfeatures[t.featid].pda.refs[row] <= t.subset[end] #if t.subset is not of the form [1:n] then this means that during the trainig, there was "no data left which fell between two splitting points" . Thus the splitting candidates collapsed
      return apply_tree_by_row(rpvector,t.left,numfeatures,charfeatures,row)
    else
      return apply_tree_by_row(rpvector,t.right,numfeatures,charfeatures,row)
    end
  end
end

function apply_tree_by_row(rpvector::Array{Array{Rulepath,1},1},t::Leaf,numfeatures::Array{pdaMod,1},charfeatures::Array{pdaMod,1},row::Int64)
	return t.fitted,firstmatch(rpvector,t.rule_path) #todo/tbd: can this be done in a better way? fristmatch is sort of a hack since I could not get findin working on the first try.... the match should be unique by definition
end

#for rankoptnew rankslost minimization
function apply_tree_by_row(numeratorval::Array{Float64,1},leaves_of_tree::Vector{Leaf},t::Node,numfeatures::Array{pdaMod,1},charfeatures::Array{pdaMod,1},leafIndexWhichIsAdjusted::Int64,premStep::Float64)
#only the fitted values for the leaf with index leafIndexWhichIsAdjusted are adjusted
length(numfeatures)>0 ? nobs=length(numfeatures[1]) : nobs=length(charfeatures[1])
  @assert (size(numfeatures,1)==0 || length(numfeatures[1])==nobs)
  @assert (size(charfeatures,1)==0 || length(charfeatures[1])==nobs)
  res=Array{Float64}(nobs)
  leafNumberVal=Array{Int64}(nobs)
  #note this can be done much more efficiently by first identifying the leaf where adjustmetns are made, and then only adjust the estimater for this leaf (the rest will be left as is)
  rpvector=[x.rule_path for x in leaves_of_tree] #rule_path of each leaf is used to identify a leaf in the following (this could probably be done more efficiently) todo/tbd
  for row=1:nobs
    res[row],leafNumberVal[row]=apply_tree_by_row(numeratorval,rpvector,t,numfeatures,charfeatures,row,leafIndexWhichIsAdjusted,premStep)
  end  
  return res,leafNumberVal
end

function apply_tree_by_row(numeratorval::Array{Float64,1},rpvector::Array{Array{Rulepath,1},1},t::Node,numfeatures::Array{pdaMod,1},charfeatures::Array{pdaMod,1},row::Int64,leafIndexWhichIsAdjusted::Int64,premStep::Float64)
#function for a single observation only
  if t.featid<0
    #character split
    if charfeatures[-t.featid].pda.refs[row] in t.subset
      return apply_tree_by_row(numeratorval,rpvector,t.left,numfeatures,charfeatures,row,leafIndexWhichIsAdjusted,premStep)
    else
      return apply_tree_by_row(numeratorval,rpvector,t.right,numfeatures,charfeatures,row,leafIndexWhichIsAdjusted,premStep)
    end
  else
    #numeric split
    if numfeatures[t.featid].pda.refs[row] <= t.subset[end] #if t.subset is not of the form [1:n] then this means that during the trainig, there was "no data left which fell between two splitting points" . Thus the splitting candidates collapsed
      return apply_tree_by_row(numeratorval,rpvector,t.left,numfeatures,charfeatures,row,leafIndexWhichIsAdjusted,premStep)
    else
      return apply_tree_by_row(numeratorval,rpvector,t.right,numfeatures,charfeatures,row,leafIndexWhichIsAdjusted,premStep)
    end
  end
end

function apply_tree_by_row(numeratorval::Array{Float64,1},rpvector::Array{Array{Rulepath,1},1},t::Leaf,numfeatures::Array{pdaMod,1},charfeatures::Array{pdaMod,1},row::Int64,leafIndexWhichIsAdjusted::Int64,premStep::Float64)
	leafNumberOfThisObservation=firstmatch(rpvector,t.rule_path) #todo/tbd: can this be done in a better way? fristmatch is sort of a hack since I could not get findin working on the first try.... the match should be unique by definition
	if leafNumberOfThisObservation==leafIndexWhichIsAdjusted 
		estimatedValue= premStep<0 ? numeratorval[row]*(1.0-premStep) : numeratorval[row]+premStep  #todo/tbd this can be done much more efficiently, we only need to assess once whether premStep is a percentage (i.e. negative) or a EURO (or other currency) amount
	else
		estimatedValue=numeratorval[row]
	end
	return estimatedValue,leafNumberOfThisObservation
end

function apply_tree_by_row(tree::BoostedTree,numfeatures::Array{pdaMod,1},charfeatures::Array{pdaMod,1})
	ntrees=size(tree.trees,1)
    length(numfeatures)>0 ? nobs=length(numfeatures[1]) : nobs=length(charfeatures[1])
    @assert (size(numfeatures,1)==0 || length(numfeatures[1])==nobs)
    @assert (size(charfeatures,1)==0 || length(charfeatures[1])==nobs)
    relativities=ones(Float64,nobs)
    scores=ones(Int64,nobs)
    est_matrix=Array{Float64}(nobs,ntrees+1)
    tree.BoolStartAtMean ? (basis=tree.meanobserved) :  (basis=1.0)
    est_matrix[:,1]=basis
    estimates=Array{Float64}(nobs)
    fill!(estimates,basis)
    p = Progress(ntrees, 2, "Applying model to validation data ...") # minimum update interval: 2 second
    for i=1:ntrees
		leaves_of_tree=create_leaves_array(tree.trees[i])
        indicated,unusedLeafNumber=apply_tree_by_row(leaves_of_tree,tree.trees[i],numfeatures,charfeatures)		
        estimates=estimates.*_moderate(indicated,tree.moderationvector[i])
        est_matrix[:,i+1]=estimates
        next!(p)
    end
    relativities=est_matrix[:,ntrees+1]./est_matrix[:,1]  
    scores=derive_scores(tree.maxRawRelativityPerScoreSorted,relativities)  
    #TODO: TBD: there is something wrong here if we have very small trees. If we test this function with the training data scores will not be equal to scores
    return hcat(scores,relativities.*basis), est_matrix,scores #we also return the scores as an integer array here
end

#for rankoptnew rankslost minimization
function apply_tree_by_row(numeratorval::Array{Float64,1},leafIndexWhichIsAdjustedArray::Array{Int64,1},leavesPerTree::Array{Array{Leaf,1},1},tree::BoostedTree,numfeatures::Array{pdaMod,1},charfeatures::Array{pdaMod,1},premStep::Float64)
    ntrees=size(tree.trees,1)
    length(numfeatures)>0 ? nobs=length(numfeatures[1]) : nobs=length(charfeatures[1])
    @assert (size(numfeatures,1)==0 || length(numfeatures[1])==nobs)
    @assert (size(charfeatures,1)==0 || length(charfeatures[1])==nobs)
    relativities=ones(Float64,nobs)
    scores=ones(Int64,nobs)
    est_matrix=Array{Float64}(nobs,ntrees+1)        
    currentEstimatedNumerator=copy(numeratorval)
	est_matrix[:,1]=currentEstimatedNumerator	
    p = Progress(ntrees, 2, "Applying model to validation data ...") # minimum update interval: 2 second    
	for i=1:ntrees
		leaves_of_tree=leavesPerTree[i]
		leafIndexWhichIsAdjusted=leafIndexWhichIsAdjustedArray[i]
		indicated,leafNumberVal=apply_tree_by_row(numeratorval,leaves_of_tree,tree.trees[i],numfeatures,charfeatures,leafIndexWhichIsAdjusted,premStep)		
		currentEstimatedNumerator=currentEstimatedNumerator.+(tree.moderationvector[i].*(indicated.-currentEstimatedNumerator)) 				
        est_matrix[:,i+1]=currentEstimatedNumerator 
        next!(p)
    end	
    relativities=est_matrix[:,ntrees+1]./est_matrix[:,1]  
    scores=derive_scores(tree.maxRawRelativityPerScoreSorted,relativities)	  
    #TODO: TBD: there is something wrong here if we have very small trees. If we test this function with the training data scores will not be equal to scores
    return hcat(scores,est_matrix[:,end]), est_matrix,scores #we also return the scores as an integer array here
end

function apply_tree_by_leaf(tree::BoostedTree,numfeatures::Array{pdaMod,1},charfeatures::Array{pdaMod,1})    
	ntrees=size(tree.trees,1)
    length(numfeatures)>0 ? nobs=length(numfeatures[1]) : nobs=length(charfeatures[1])
    @assert (size(numfeatures,1)==0 || length(numfeatures[1])==nobs)
    @assert (size(charfeatures,1)==0 || length(charfeatures[1])==nobs)
    relativities=ones(Float64,nobs)
    scores=ones(Int64,nobs)
    est_matrix=Array{Float64}(nobs,ntrees+1)
    tree.BoolStartAtMean ? (basis=tree.meanobserved) :  (basis=1.0)
    est_matrix[:,1]=basis
    estimates=Array{Float64}(nobs)
    fill!(estimates,basis)
    p = Progress(ntrees, 2, "Applying model to validation data ...") # minimum update interval: 2 second
    for i=1:ntrees
		leaves_of_tree=create_leaves_array(tree.trees[i])
        indicated,unusedLeafNumber=apply_tree_by_leaf(leaves_of_tree,tree.trees[i],numfeatures,charfeatures)		
        estimates=estimates.*_moderate(indicated,tree.moderationvector[i])
        est_matrix[:,i+1]=estimates
        next!(p)
    end
    relativities=est_matrix[:,ntrees+1]./est_matrix[:,1]  
    scores=derive_scores(tree.maxRawRelativityPerScoreSorted,relativities)      
    return hcat(scores,relativities*basis), est_matrix,scores #we also return the scores as an integer array here
end

function apply_tree_by_leaf_iteration(rpvector::Array{Array{Rulepath,1},1},t::Node,numfeatures::Array{pdaMod,1},charfeatures::Array{pdaMod,1},numeratorval::Array{Float64,1},leafIndexWhichIsAdjusted::Int64,premStep::Float64)	 
	 length(numfeatures)>0 ? nobs=length(numfeatures[1]) : nobs=length(charfeatures[1])
	 if nobs==0
		return Array{Float64}(0),Array{Int64}(0)
	 end
	 featid=t.featid
	 subset=t.subset
	 idxl,idxr,idxlInteger,idxrInteger=lrIndices(charfeatures,numfeatures,subset,featid)
	 charfeaturesl,charfeaturesr=split_pda_left_right(idxlInteger,idxrInteger,charfeatures)
	 numfeaturesl,numfeaturesr=split_pda_left_right(idxlInteger,idxrInteger,numfeatures)	 
	 numeratorvall=numeratorval[idxlInteger]
	 numeratorvalr=numeratorval[idxrInteger]
	 fitr,leafr=apply_tree_by_leaf_iteration(rpvector,t.right,numfeaturesr,charfeaturesr,numeratorvalr,leafIndexWhichIsAdjusted,premStep)
	 fitl,leafl=apply_tree_by_leaf_iteration(rpvector,t.left,numfeaturesl,charfeaturesl,numeratorvall,leafIndexWhichIsAdjusted,premStep)
	 
	fit=Array{Float64}(nobs)
	leaf=Array{Int64}(nobs)
	loopl=loopr=1	
	for jj=1:nobs
	if idxl[jj]
		fit[jj]=fitl[loopl]
		leaf[jj]=leafl[loopl]
		loopl+=1
	else
		fit[jj]=fitr[loopr]
		leaf[jj]=leafr[loopr]
		loopr+=1
		end
	end		
 return fit,leaf 
 end
 
function apply_tree_by_leaf_iteration(rpvector::Array{Array{Rulepath,1},1},t::Leaf,numfeatures::Array{pdaMod,1},charfeatures::Array{pdaMod,1},numeratorval::Array{Float64,1},leafIndexWhichIsAdjusted::Int64,premStep::Float64)
	length(numfeatures)>0 ? nobs=length(numfeatures[1]) : nobs=length(charfeatures[1])
	leafNumberOfThisObservation=firstmatch(rpvector,t.rule_path) #todo/tbd: can this be done in a better way? fristmatch is sort of a hack since I could not get findin working on the first try.... the match should be unique by definition
	if leafNumberOfThisObservation==leafIndexWhichIsAdjusted 
		fit= premStep<0 ? numeratorval.*(1.0-premStep) : numeratorval.+premStep  #todo/tbd this can be done much more efficiently, we only need to assess once whether premStep is a percentage (i.e. negative) or a EURO (or other currency) amount
	else
		fit=numeratorval
	end	
	nr=Array{Int64}(nobs)
	fill!(nr,leafNumberOfThisObservation)
	return fit,nr	
end	 

#Bagging
function apply_tree_by_leaf(rpvector::Array{Array{Rulepath,1},1},t::Node,numfeatures::Array{pdaMod,1},charfeatures::Array{pdaMod,1})
	return apply_tree_by_leaf_iteration(rpvector,t,numfeatures,charfeatures)
end

function apply_tree_by_leaf(idx::Vector{Int},leaves_of_tree::Array{Leaf,1},t::Tree,features::DataFrame)	
	return apply_tree_by_leaf(idx,leaves_of_tree,t.rootnode,features)
end

function apply_tree_by_leaf!(fit::Vector{Float64},leaf::Vector{Int},idx::Vector{Int},t::Union{Leaf,Node},features::DataFrame)  
  nobs=size(features,1) 
  @assert nobs==length(fit)==length(leaf)
  #rpvector=[x.rule_path for x in leaves_of_tree] #rule_path of each leaf is used to identify a leaf in the following (this could probably be done more efficiently) todo/tbd
  apply_tree_by_leaf_iteration!(idx,t,features,fit,leaf)  
return nothing
end
  
function apply_tree_by_leaf(idx::Vector{Int},t::Union{Leaf,Node},features::DataFrame)  
  #as we always supply an index here, the function will always return two arrays of length size(features,2)
  #thus the "training" observations will be empty if this function is called on a validation index only!!
  nobs=size(features,1) 
  #rpvector=[x.rule_path for x in leaves_of_tree] #rule_path of each leaf is used to identify a leaf in the following (this could probably be done more efficiently) todo/tbd
  fit=zeros(Float64,nobs)
  leaf=zeros(Int,nobs)
  apply_tree_by_leaf_iteration!(idx,t,features,fit,leaf)  
  return fit,leaf
end

function apply_tree_by_leaf_iteration!(idx::Vector{Int},t::Node,features::DataFrame,fit::Vector{Float64},leaf::Vector{Int})
	 if length(idx)==0
       info("DTM: No data (an empty index) was provided to apply tree function")  #@show "mix"
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
 
function apply_tree_by_row(tree::BoostedTree,numfeatures::Array{pdaMod,1},charfeatures::Array{pdaMod,1},labels_orig_val::Array{Float64,1})
	ntrees=size(tree.trees,1)
    length(numfeatures)>0 ? nobs=length(numfeatures[1]) : nobs=length(charfeatures[1])
    @assert (size(numfeatures,1)==0 || length(numfeatures[1])==nobs)
    @assert (size(charfeatures,1)==0 || length(charfeatures[1])==nobs)
    relativities=ones(Float64,nobs)
    scores=ones(Int64,nobs)
    est_matrix=Array{Float64}(nobs,ntrees+1)
    #tree.BoolStartAtMean ? (basis=tree.meanobserved) :  (basis=1.0)
    estimates=copy(labels_orig_val)
	  est_matrix[:,1]=copy(estimates)
    p = Progress(ntrees, 2, "Applying model to validation data ...") # minimum update interval: 2 second
    for i=1:ntrees
		indicated=apply_tree_by_row(tree.trees[i],numfeatures,charfeatures)
        estimates=estimates.+indicated
        est_matrix[:,i+1]=estimates
        next!(p)
    end
	relativities=est_matrix[:,ntrees+1]./est_matrix[:,1]
    scores=derive_scores(tree.maxRawRelativityPerScoreSortedSorted,relativities)
    #TODO: TBD: there is something wrong here if we have very small trees. If we test this function with the training data scores will not be equal to scores
    return hcat(scores,relativities), est_matrix,scores #we also return the scores as an integer array here
end

function apply_tree_additiv(tree::BoostedTree,numfeatures::Array{pdaMod,1},charfeatures::Array{pdaMod,1},labels_orig::Array{Float64,1})
    ntrees=size(tree.trees,1)
    length(numfeatures)>0 ? nobs=length(numfeatures[1]) : nobs=length(charfeatures[1])
    relativities=zeros(Float64,nobs)
    scores=ones(Int64,nobs)
    est_matrix=Array{Float64}(nobs,ntrees+1)
    est_matrix[:,1]=copy(labels_orig)
    #fill!(estimates,labels_orig)
    p = Progress(ntrees, 2, "Applying model to validation data ...") # minimum update interval: 2 second
    for i=1:ntrees
        indicated=apply_tree_by_row(tree.trees[i],numfeatures,charfeatures)
        relativities=relativities+indicated.*tree.moderationvector[i]
        est_matrix[:,i+1]=relativities+labels_orig
        next!(p)
    end
    scores=derive_scores(tree.maxRawRelativityPerScoreSortedSorted,relativities)
    return hcat(relativities+labels_orig,relativities,scores), est_matrix
end
