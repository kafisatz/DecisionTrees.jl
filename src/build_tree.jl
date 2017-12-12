#trnidx,validx,candMatWOMaxValues,mappings,deepcopy(sett),actualNumerator,estimatedNumerator,weight,features,sampleSizeCanBeNEGATIVE,abssampleSize,sampleVector)																																																	  
function sample_data_and_build_tree!(trnidx::Vector{Int},validx::Vector{Int},fitted_values_all_data_this_vector_is_modified_by_build_tree::Vector{Float64},candMatWOMaxValues::Array{Array{Float64,1},1},mappings::Array{Array{String,1},1},settings::ModelSettings,numerator::Array{Float64,1},denominator::Array{Float64,1},weight::Array{Float64,1},features,sampleSizeCanBeNEGATIVE,abssampleSize,sampleVector)
#todo: automatically ignore varialbes with no splitting points!
#NOTE: this function will call build_tree!()
	#fitted_values_all_data_this_vector_is_modified_by_build_tree=zeros(numerator)
	if abs(settings.subsampling_prop)>=1.0
		n=build_tree!(trnidx,validx,candMatWOMaxValues,mappings,settings,numerator,denominator,weight,features,fitted_values_all_data_this_vector_is_modified_by_build_tree)
		return n 
	else
        reused_fitted_leafnr_vector=zeros(Int,length(weight))
		#do subsampling
		#sampleSizeCanBeNEGATIVE=convert(Int,round(settings.subsampling_prop*length(denominator)))
		#note: this can probably be done more efficiently (see also bagging)
		#abssampleSize,sampleVector,num,denom,w,numf,charf,ooBagnum,ooBagdenom,ooBagw,ooBagnumf,ooBagcharf,ooBagsize=initBoostrapSample(sampleSizeCanBeNEGATIVE,numerator,denominator,weight,charfeatures,numfeatures)
		unusedSamplePart=sampleData!(trnidx,sampleSizeCanBeNEGATIVE,sampleVector)
		#build tree	on subsample of data
			thistree=build_tree!(sampleVector,validx,candMatWOMaxValues,mappings,deepcopy(settings),numerator,denominator,weight,features,fitted_values_all_data_this_vector_is_modified_by_build_tree)
	#apply tree to ooBagSample
			leaves_of_tree=create_leaves_array(thistree.rootnode)
			#ooBagEstimates,leafNrooBag=eval(parse(settings.chosen_apply_tree_fn))(unusedSamplePart,leaves_of_tree,thistree.rootnode,features)
            apply_tree_by_leaf!(fitted_values_all_data_this_vector_is_modified_by_build_tree,reused_fitted_leafnr_vector,unusedSamplePart,thistree.rootnode,features)
		#Determine Goodness of Fit
			#I think it is best to consider the goodness of fit per leaf here
			#out of bag performance should be measured here, however this is not coded yet (see also bagging!)
				#ooBagEstNumerator=ooBagEstimates.*ooBagdenom
				#fittedPerLeaf=Float64[x.fitted for x in leaves_of_tree]
				#ooBagcntperLeaf,ooBagsumnumeratorEST,ooBagsumnumerator,ooBagsumdenominator,ooBagsumweight=aggregateByLeafNr(leafNrooBag,ooBagEstNumerator,ooBagnum,ooBagdenom,ooBagw)
				#todo tbd check this!
				#thisError=calcErrorStats(fittedPerLeaf,ooBagsumnumerator./ooBagsumdenominator,ooBagsumweight)
				#thisError=goodnessOfFit(ooBagEstimates,ooBagnum,ooBagdenom,ooBagw,meanobservedvalue,ngroupsInput)
		thisError=ErrorStats() #Empty/Default Error Stats
		return thistree
	end
end

function build_tree!(trnidx::Vector{Int},validx::Vector{Int},candMatWOMaxValues::Array{Array{Float64,1},1},mappings::Array{Array{String,1},1},settings::ModelSettings,numerator::Array{Float64,1},denominator::Array{Float64,1},weight::Array{Float64,1},features,fitted_values::Vector{Float64})
	intVarsUsed,inds,minweightcalculated,nDepthToStartParallelization=some_tree_settings(trnidx,validx,settings.fixedinds,candMatWOMaxValues,mappings,settings.minw,weight,settings.subsampling_features_prop,size(features,2))
	settings.minw=minweightcalculated #update minw
	@assert length(inds)>0 "Error: no features were selected length(num_inds)=$(length(num_inds)), length(char_inds)=$(length(char_inds))"
	settings.nDepthToStartParallelization=nDepthToStartParallelization #update nDepthToStartParallelization
	empty_xl_data=ExcelData(Array{ExcelSheet}(0),Array{Chart}(0))
	resultingTree=Tree(deepcopy(emptyNode),intVarsUsed,candMatWOMaxValues,mappings,inds,settings,deepcopy(empty_xl_data))
	resultingTree.rootnode=build_tree_iteration!(trnidx,validx,settings,resultingTree,numerator,denominator,weight,features,0,settings.randomw,Array{Rulepath}(0),settings.parallel_tree_construction,myid(),fitted_values)
	#set Leaf Numbers
	set_leaf_numbers!(resultingTree)
	return resultingTree
end

function build_tree_iteration!(trnidx::Vector{Int},validx::Vector{Int},settings::ModelSettings,thisTree::Tree,numerator::Array{Float64,1},denominator::Array{Float64,1},weight::Array{Float64,1},features::DataFrame,
									depth::Int64,randomweight::Float64,parent_rp::Array{Rulepath,1},parallel_tree_construction::Bool,parentid::Int64,fitted_values::Vector{Float64})
	#!!!! the current concept forsees that features is always the FULL DataFrame
	boolRandomizeOnlySplitAtTopNode=settings.boolRandomizeOnlySplitAtTopNode
	local inds
	#if fixedinds are provided we never change the selection of features for the whole tree! (bagged boosting model)
	if boolRandomizeOnlySplitAtTopNode||(length(settings.fixedinds)>0)
		#these settings are randomly chosen by some_tree_settings in this case
		inds=thisTree.inds_considered
	else
		#char_inds and num_inds need to be randomly picked for each iteration
		inds=randomFeatureSelection(size(features,2),settings.subsampling_features_prop)
	end
	intVarsUsed=thisTree.intVarsUsed

	minweight=settings.minw
	crit=settings.crit
	parallel_level_threshold=settings.parallel_level_threshold
	parallel_weight_threshold=settings.parallel_weight_threshold
	nDepthToStartParallelization=settings.nDepthToStartParallelization
	spawnsmaller=settings.spawnsmaller
	recursivespawning=settings.recursivespawning
	pminweightfactor=settings.pminweightfactor
	pminrelpctsize=settings.pminrelpctsize
	pflipspawnsmalllargedepth=settings.pflipspawnsmalllargedepth
	catSortByThreshold=settings.catSortByThreshold
	catSortBy=settings.catSortBy

	nobs=size(numerator,1)
	boolMDFPerLeaf=false #this is disabled here
    fnames=names(features)
	#@code_warntype _split(settings.number_of_num_features,trnidx,validx,numerator,denominator,weight,fnames,features, minweight, depth,randomweight,crit,parallel_level_threshold,parallel_weight_threshold,inds,catSortByThreshold,catSortBy)
	best_split = _split(settings.number_of_num_features,trnidx,validx,numerator,denominator,weight,fnames,features, minweight, depth,randomweight,crit,parallel_level_threshold,parallel_weight_threshold,inds,catSortByThreshold,catSortBy)
	id=best_split.featid
	subset=best_split.subset
	fname=best_split.featurename
	id2=best_split.featid_new_positive
	
	tmpsn=sum(view(numerator,trnidx))
    tmpsd=sum(view(denominator,trnidx))
    #todo/tbd check if this "sort!(subset)" is necessary and find out why.... is it intended?
	sort!(subset)

	#check if no split was found
  #this can happen (EVEN AT THE TOP NODE) if subsampling (of data and or features) is enabled: it may be that there is in fact no split possible (e.g. if all variables are constant)
  if id == 0
	if (depth==0)
		printover("No split was found at the top node. This may cause errors later on in the code!")
	end
		trnsumn=sum(view(numerator,trnidx))
		trnsumd=sum(view(denominator,trnidx))
        mean_observed=trnsumn/trnsumd;this_leaf_size=sum(view(weight,trnidx))
		fitted_increment=mean_observed
      return Leaf(length(trnidx),mean_observed,fitted_increment, this_leaf_size, depth, parent_rp,trnsumn,trnsumd,-1)::Leaf
  end;
#Rulepath
		this_left_rp=deepcopy(parent_rp);
		this_right_rp=deepcopy(parent_rp);
		push!(this_left_rp,Rulepath(id,subset,true));
		push!(this_right_rp,Rulepath(id,subset,false));
		if id<0
			#split by character variable       
            for u in subset
                intVarsUsed[-id+settings.number_of_num_features][u]+=1
            end
        else
            #split by numeric variable
            #todo check performance of this and improve
            intVarsUsed[id][subset[end]]+=1 #set this value to true to indicate, that is is used by the tree
        end
        column=features[id2]
        #indexOfMatchedUint8s=findin(column.ids,subset)
    	matched_strings=column.pool[subset]
        l,r=lrIndices(trnidx,column,subset)

    countl=size(l,1)
    countr=size(r,1)
	sumwl=sum(view(weight,l))
	sumwr=sum(view(weight,r))

	leftchildwillbefurthersplit=sumwl<2*minweight
	rightchildwillbefurthersplit=sumwr<2*minweight
	(sumwr>sumwl) ? boolSpawnLeft=false : boolSpawnLeft=true
    if spawnsmaller;boolSpawnLeft=!boolSpawnLeft;end;
    #flip spawning left or rigth childer at a certain depth
    if (depth>=pflipspawnsmalllargedepth);boolSpawnLeft=!boolSpawnLeft;end;

    boolRChildAllowedToFurtherSpawnProcesses=(recursivespawning&&(sumwr>(minweight))&&(pminrelpctsize<(sumwr/(sumwr+sumwl))))
    boolLChildAllowedToFurtherSpawnProcesses=(recursivespawning&&(sumwl>(minweight))&&(pminrelpctsize<(sumwl/(sumwr+sumwl))))
  	boolRandomizeOnlySplitAtTopNode ? newrandomweight=copy(randomweight) : newrandomweight=0.0

    #todo/tbd check the different if/then here, I think a few things are not quite accurate
	#Also, we probably do not need any parallelization for the construction of a single tree-> we should disable this functionality
	trnsumnl=sum(view(numerator,l))
	trnsumdl=sum(view(denominator,l))
	trnsumnr=sum(view(numerator,r))
	trnsumdr=sum(view(denominator,r))
	mean_observedl=trnsumnl/trnsumdl	
	mean_observedr=trnsumnr/trnsumdr
	fitted_labelsl=mean_observedl;
	fitted_labelsr=mean_observedr;
	if !(parallel_tree_construction)||(max(sumwr,sumwl)<minweight)
					if leftchildwillbefurthersplit			
						fill_some_elements!(fitted_values,l,mean_observedl)
							remote_ref_build_tree_leftchild = Leaf(length(l),mean_observedl,fitted_labelsl,sumwl, depth+1, this_left_rp,trnsumnl,trnsumdl,-1)
          else
            remote_ref_build_tree_leftchild = build_tree_iteration!(l,validx,settings,thisTree,numerator,denominator,weight,features,depth+1,newrandomweight,this_left_rp,boolLChildAllowedToFurtherSpawnProcesses,myid(),fitted_values)
		      end
					if rightchildwillbefurthersplit			    
						fill_some_elements!(fitted_values,r,mean_observedr)
            remote_ref_build_tree_rightchild = Leaf(length(r),mean_observedr,fitted_labelsr,sumwr, depth+1, this_right_rp,trnsumnr,trnsumdr,-1)
          else
            remote_ref_build_tree_rightchild = build_tree_iteration!(r,validx,settings,thisTree,numerator,denominator,weight,features,depth+1,newrandomweight,this_right_rp,boolRChildAllowedToFurtherSpawnProcesses,myid(),fitted_values)
          end
      else
       #here we spawn a process for the smaller "child"
        if (!boolSpawnLeft)
        #here it is important that the @spawn is executed before the build_tree_iteration! call!
						if rightchildwillbefurthersplit
							fill_some_elements!(fitted_values,r,mean_observedr)
             	  remote_ref_build_tree_rightchild = Leaf(length(r),mean_observedr,fitted_labelsr,sumwr, depth+1, this_right_rp,trnsumnr,trnsumdr,-1)
            else
              remote_ref_build_tree_rightchild = @spawn build_tree_iteration!(r,validx,settings,thisTree,numerator,denominator,weight,features,depth+1,newrandomweight,this_right_rp,boolRChildAllowedToFurtherSpawnProcesses,myid(),fitted_values)
            end
						if leftchildwillbefurthersplit
							fill_some_elements!(fitted_values,l,mean_observedl)
						  remote_ref_build_tree_leftchild = Leaf(length(l),mean_observedl,fitted_labelsl,sumwl, depth+1, this_left_rp,trnsumnl,trnsumdl,-1)
            else
             remote_ref_build_tree_leftchild = build_tree_iteration!(l,validx,settings,thisTree,numerator,denominator,weight,features,depth+1,newrandomweight,this_left_rp,boolLChildAllowedToFurtherSpawnProcesses,myid(),fitted_values)
            end
        else
						if leftchildwillbefurthersplit
							fill_some_elements!(fitted_values,l,mean_observedl)
						  remote_ref_build_tree_leftchild = Leaf(length(l),mean_observedl,fitted_labelsl,sumwl, depth+1, this_left_rp,trnsumnl,trnsumdl,-1)
            else
              remote_ref_build_tree_leftchild = @spawn build_tree_iteration!(l,validx,settings,thisTree,numerator,denominator,weight,features,depth+1,newrandomweight,this_left_rp,boolLChildAllowedToFurtherSpawnProcesses,myid(),fitted_values)
            end
						if rightchildwillbefurthersplit
							fill_some_elements!(fitted_values,r,mean_observedr)
               remote_ref_build_tree_rightchild = Leaf(length(r),mean_observedr,fitted_labelsr,sumwr, depth+1, this_right_rp,trnsumnr,trnsumdr,-1)
            else
              remote_ref_build_tree_rightchild = build_tree_iteration!(r,validx,settings,thisTree,numerator,denominator,weight,features,depth+1,newrandomweight,this_right_rp,boolRChildAllowedToFurtherSpawnProcesses,myid(),fitted_values)
            end
        end
     end

  fetched_left=fetch(remote_ref_build_tree_leftchild)
  fetched_right=fetch(remote_ref_build_tree_rightchild)
  #fitted_values=mergeLeftAndRightFittedValues(fitted_valuesl,fitted_valuesr,left)
	@assert (id<0)||(maximum(subset)==subset[end]) #for numeric variables subset should generally be of the form [1:n], if not this means that the data was split in a way such that several candidates collapsed as there was no more data between to (or more) candidates.
    id2 = id < 0 ? abs(id)+settings.number_of_num_features : id
return Node(id,id2,subset,fetched_left,fetched_right,parent_rp)::Node
end

function _split(number_of_num_features::Int,trnidx::Vector{Int},validx::Vector{Int},numerator::Array{Float64,1},denominator::Array{Float64,1},weight::Array{Float64,1},fnames::Vector{Symbol},features, minweight::Float64, depth::Int64,randomweight::Float64,crit::SplittingCriterion,parallel_level_threshold::Int64=9999999999,parallel_weight_threshold::Int64=9999999999,inds::Array{Int64,1}=Array{Int64}(0),catSortByThreshold::Int64=8,catSortBy::SortBy=SORTBYMEAN)
  #This function selects the maximal possible split defined by crit (thus depending on the impurity function, we need to put a minus sign in front of it)
		tmpsz::Int=0
		if sum(view(weight,trnidx))<2*minweight;
			return const_default_splitdef::Splitdef;
		end
		tmp_splitlist=Array{Splitdef}(0)
		for i in inds
			#ATTENTION: for char variables we pass the variable i with a negative sing!!
			#this allows us to distinguish whether we are working on a char or num variable later on
			modified_i = eltype(features[i])<:AbstractString ? number_of_num_features - i  :  i 
			tmplist::Vector{Splitdef}=_split_feature(number_of_num_features,trnidx,validx,numerator,denominator,weight,fnames[i],features[i],minweight,crit,modified_i,randomweight,catSortByThreshold,catSortBy)::Vector{Splitdef}
			append!(tmp_splitlist,tmplist)
		end

		#pick best "valid" (minweight feasible) split
		if !(size(tmp_splitlist,1)>0)
			return const_default_splitdef
		else 
				tmp_splitlist=subset_splitlist(tmp_splitlist,minweight)
				if !(size(tmp_splitlist,1)>0)
					return const_default_splitdef
				else
            splitlist_sorted::Vector{Splitdef}=sort_splitlist(tmp_splitlist)::Vector{Splitdef}
						tmpsz=size(splitlist_sorted,1)::Int
						if tmpsz>0 #can it be that splitlist_sorted is empty even though tmp_splitlist was not? I think so!
							#randomize choice
							spl::Splitdef=splitlist_sorted[1]
							if randomweight>0 #((depth==0) && (randomweight>0))
								rnd=Int64(max(1,min(tmpsz,ceil(rand()*randomweight*tmpsz)))) #I am not sure if it can happen that rnd becomes 0 if we do not impose the max(1,...) condition. But it seems safe to enforce it in case an incredibly small random number is generated
								spl=splitlist_sorted[rnd]::Splitdef
							else #randomweight<=0
								#deterministic choice (greedy)
								spl=splitlist_sorted[1]::Splitdef
							end #"if condition randomweight>0"

							if isfinite(spl.splitvalue) && (spl.featid<0)				
								#For Character variables: The split can be defined by the subset or its complement, We choose to define it via the subset which defines the smaller child node such that when new values arrive (in out of sample testing) they will "go with the larger child node"
								if spl.weightl>spl.weightr
										lvls::Vector{UInt8}=unique(view(features[spl.featid_new_positive].refs,trnidx))::Vector{UInt8}
										spl=Splitdef(spl.featid,spl.featid_new_positive,spl.featurename,setdiff(lvls,spl.subset),spl.splitvalue,spl.weightl,spl.weightr)::Splitdef
								end
							end #isfinite(best_value_split) && (best[1]<0)
						end #size(splitlist_sorted,1)>0
        end #size(tmp_splitlist,1)>0
		end #size(tmp_splitlist,1)>0
		
return spl::Splitdef
end


#new approach; first summarize by label, then iterate over the gray code such that only one category needs to switch classes and "online" update the metrics
function _split_feature(number_of_num_features::Int,trnidx::Vector{Int},validx::Vector{Int},numerator::Array{Float64,1},denominator::Array{Float64,1},weight::Array{Float64,1},fname::Symbol,features,minweight::Float64,crit::SplittingCriterion,feature_column_id::Int64,randomweight::Float64,catSortByThreshold::Int64,catSortBy::SortBy)
crit_type=typeof(crit)
#This function is now for numeric and character variables!
#feature_column_id is negative in case of character variables
best_value=-Inf
best_thresh=best_wl=best_wr=NaN
best_subset=Array{UInt8}(0)
trnfeatures=view(features,trnidx)
elt=eltype(trnfeatures.parent.refs)
#if randomweight>0;this_splitlist=Array{Splitdef}(0);end;
  labellist_sorted=collect(one(elt):convert(elt,length(trnfeatures.parent.pool))) #this used to be levels(features) #this also contains the val feature levels here! It is considerably faster than levels(view) \factor 100 or so
	# THIS IS CRITICAL
	# THE WAY A PooledArray IS CONSTRUCTED, WE WILL ALWAYS HAVE
	# pda.pool[2] == "some string" -> any pda.refs[x].==0x02 will be "some string" 
	#todo/tbd countsort here might be obsolete: we should check if levels is always sorted by construction
  #also we will later sort the labels in a different order anyway (then again the list probably needs to be sorted in the natural manner such that build_listOfMeanResponse is working properly)
  if size(labellist_sorted,1) <= 1
    return collect(Array{Splitdef}(0))
  else
	countsort!(labellist_sorted)
    #here I am (somewhat) misusing multiple dispatch since I was too lazy to parametrize this at the moment (todo/tbd in the future)
	if (crit_type==DifferenceSplit||crit_type==PoissonDevianceSplit||crit_type==MaxValueSplit||crit_type==MaxMinusValueSplit)
		labellist,sumnumerator,sumdenominator,sumweight,countlistfloat=build_listOfMeanResponse(crit,trnidx,validx,numerator,denominator,weight,trnfeatures,labellist_sorted,minweight)
	elseif (crit_type==NormalDevianceSplit)
		info("in the works...")
		labellist,sumnumerator,sumdenominator,sumweight,countlistfloat,moments_per_pdaclass=build_listOfMeanResponse(crit,numerator,denominator,weight,trnfeatures,labellist_sorted,minweight)
	else
		error("Invalid Splitting criterion $(crit)")
	end
  #todo/tbd
  #here we can introduce the possiblity to sort the labellist (e.g. by meanobserved or median (to be calculated)).
  #then we could only loop through the "increasing list" of sorted labels (instead of doing the 2^ncategories exhaustive search (bitflip_graycode_subsets))
  if feature_column_id>0 #id>0 -> we are working on a numeric column
		#only consider to split at the candidate split points
		subs=increasing_subsets(labellist)
  else #id<0 -> we are working on a character column
		#distinguish between exhaustive and "increasing" search for split point
	if size(labellist_sorted,1)>catSortByThreshold
		if (crit_type==DifferenceSplit||crit_type==MaxValueSplit||crit_type==MaxMinusValueSplit)
			sortlists!(catSortBy,labellist,sumnumerator,sumdenominator,sumweight,countlistfloat) #catSortBy::SortBy=SORTBYMEAN
		elseif (crit_type==NormalDevianceSplit)
			sortlists!(catSortBy,labellist,sumnumerator,sumdenominator,sumweight,countlistfloat,moments_per_pdaclass)
		end
		subs=increasing_subsets(labellist)
	else
	#perform exhaustive search
		subs=bitflip_graycode_subsetsHALF(labellist)
	end
	end

	if randomweight==0.0
	if (crit_type==DifferenceSplit||crit_type==MaxValueSplit||crit_type==MaxMinusValueSplit)
		tmp_result=calculateSplitValue(crit,fname,number_of_num_features,labellist,sumnumerator,sumdenominator,sumweight,countlistfloat,minweight,subs)
	elseif (crit_type==NormalDevianceSplit)
		tmp_result=calculateSplitValue(crit,fname,number_of_num_features,labellist,sumnumerator,sumdenominator,sumweight,countlistfloat,minweight,subs,moments_per_pdaclass)
	elseif (crit_type==PoissonDevianceSplit)
		tmp_result=calculateSplitValue(crit,fname,number_of_num_features,labellist,sumnumerator,sumdenominator,sumweight,countlistfloat,minweight,subs,numerator,denominator,weight,trnfeatures)
	end
    #error("the next step may be incorrect")
    if isfinite(tmp_result[1])
        feature_column_id2 = feature_column_id < 0 ? abs(feature_column_id) + number_of_num_features : feature_column_id
      return [Splitdef(feature_column_id,feature_column_id2,fname,tmp_result[2],tmp_result[1],tmp_result[3],tmp_result[4])]
    else
      return collect(Array{Splitdef}(0))
    end
  else
  #randomweight>0
	if (crit_type==DifferenceSplit||crit_type==MaxValueSplit||crit_type==MaxMinusValueSplit)
		tmpres=calculateSplitValue(crit,fname,number_of_num_features,labellist,sumnumerator,sumdenominator,sumweight,countlistfloat,minweight,subs,feature_column_id)
	elseif (crit_type==NormalDevianceSplit)
		tmpres=calculateSplitValue(crit,fname,number_of_num_features,labellist,sumnumerator,sumdenominator,sumweight,countlistfloat,minweight,subs,feature_column_id,moments_per_pdaclass)
	elseif (crit_type==PoissonDevianceSplit)
		tmp_result=calculateSplitValue(crit,fname,number_of_num_features,labellist,sumnumerator,sumdenominator,sumweight,countlistfloat,minweight,subs,numerator,denominator,weight,trnfeatures,feature_column_id)
		#error("PoissonDevianceSplit is not yet implemented for randomw>0")
	end
	return tmpres
  end
end
end
