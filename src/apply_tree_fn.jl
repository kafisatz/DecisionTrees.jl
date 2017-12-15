#Bagging
function apply_tree_by_leaf(idx::Vector{Int},leaves_of_tree::Array{Leaf,1},t::Tree,features::DataFrame)	
	return apply_tree_by_leaf(idx,leaves_of_tree,t.rootnode,features)
end

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