warn("BK: this is in the works...")

#=
	include("analyze_pw_only_data_init_only.jl")


	info("use this function, specifically xlogy")
	info("""
	dres[i] = 2 * (xlogy(yi, yi / μi) - (yi - μi))

	for z in chlid_left
		# yi = 'number of claims for observation i'
		#\mu_i = mean number of claims in left child #<- determined in the first pass over the data
		deviance_left += (xlogy(yi, yi / μi) - (yi - μi))
	end
	deviance_left * = 2

	#when considering different splits -> consider the grey code 'shuffeling'

	# let us say we have 255 possible splits
	# this yields a vector mean_left of length 255 and a vector mean_right of length 255
	# for each mean left and mean right loop over the data and calculate deviance_left and deviance_right (which will also be both vectors of length 255!)
	# then we will have PAIRs deviance_left[43] and deviance_left[43] for the potential split number '43' (whatever that might be)
	# one of these pairs will correspond to the optimal split!

	#it might be that aggregate_data_pois_deviance does most (if not all) of the work by calculating deviance_left and deviance_right for each possible split!

	""")

 # import DTM: DTSubsets, PoissonDevianceSplit, countsort!,build_listOfMeanResponse,sortlists!,increasing_subsets,bitflip_graycode_subsetsHALF

		numerator=numeratortrn
		denominator=denominatortrn
		weight=weighttrn
		feature_column_id=-5
		features=trn_charfeatures_PDA[-feature_column_id]
		minweight=sett.minw
		randomweight=sett.randomw
		catSortByThreshold=sett.catSortByThreshold
		catSortBy=sett.catSortBy
		labellist_sorted=levels(features)

		crit=sett.crit
		#numerator::Array{Float64,1},denominator::Array{Float64,1},weight::Array{Float64,1},features::pdaMod,minweight::Float64,crit::SplittingCriterion,feature_column_id::Int64,randomweight::Float64,catSortByThreshold::Int64,catSortBy::SortBy)
		crit_type=typeof(crit)

		crit_type=typeof(crit)
		#This function is now for numeric and character variables!
		#feature_column_id is negative in case of character variables
		best_value=-Inf
		best_thresh=best_wl=best_wr=NaN
		best_subset=Array{UInt8}(0)
		#if randomweight>0;this_splitlist=Array{Splitdef}(0);end;
		labellist_sorted=levels(features)


		  countsort!(labellist_sorted)
			#here I am (somewhat) misusing multiple dispatch since I was too lazy to parametrize this at the moment (todo/tbd in the future)

				labellist,sumnumerator,sumdenominator,sumweight,countlistfloat=build_listOfMeanResponse(crit,numerator,denominator,weight,features,labellist_sorted,minweight)
			#@show feature_column_id
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
				tmp_result=calculateSplitValue(crit,labellist,sumnumerator,sumdenominator,sumweight,countlistfloat,minweight,subs)
			elseif (crit_type==NormalDevianceSplit)

=#


#function calculateSplitValue(a::DifferxxxxenceSplit,fname::Symbol,number_of_num_features::Int,labellist::Vector{T},sumnumerator::Array{Float64,1},sumdenominator::Array{Float64,1},sumweight::Array{Float64,1},countlistfloat::Array{Float64,1},minweight::Float64,subs::DTSubsets) where T<:Unsigned
function calculateSplitValue(a::PoissonDevianceSplit,fname::Symbol,number_of_num_features::Int,labellist::Vector{T},sumnumerator::Array{Float64,1},sumdenominator::Array{Float64,1},sumweight::Array{Float64,1},countlistfloat::Array{Float64,1},minweight::Float64,subs::DTSubsets,numerator::Array{Float64},denominator::Array{Float64},weight::Array{Float64},features) where T<:Unsigned
#  warn("xlogy will fail terribly if one of the arguments is zero! we need to amend the code to handle that case properly!")
  #here randomweight==0
  #for subsets, exhaustive search with flipping members (gray code) or "increasing" subset search ({1}, {1,2}, {1,2,3}, .... {1,2,3, ....., n-1,2})
  #all input lists (labellist,sumnumerator,sumdenominator,sumweight,countlistfloat) need to be sorted in the same manner

  #convention, in the beginning everything is on the right side
  elementsInLeftChildBV=BitVector(length(labellist));
  fill!(elementsInLeftChildBV,false) #indicates which classes belong to right child
  chosen_subset_bitarray=BitVector(0)
  numtot=sum(sumnumerator)
  denomtot=sum(sumdenominator)
  weighttot=sum(sumweight)
  weighttot_minw=weighttot-minweight

  #step 1: find out how 'big' subs is
  subs_size=length(subs)
    #@assert subs_size<255 #todo tbd: we can remove this later on

  #pass on over the data:
  #get mean of left chlid and right child for each possible split
  #warn("find out if we need/want freq or claimcount here!")
  meansl=zeros(subs_size) #mean (frequency or possibly claim count) of the left node in case we go for split 'number' i
  meansr=zeros(subs_size) #same for right child
  #same for weights
  weightsl=zeros(subs_size)
  weightsr=zeros(subs_size)

  sumnl=sumwl=sumdl=0.0

  counter=1
  for i in subs
    #@show "vv=$(i)" #should be the element which switches
    #i is an index, it indicates which element of labellist will flip sides for the next calculation
    if elementsInLeftChildBV[i]
    #update elementsInLeftChildBV, i.e. toggle the value of the ith component, $=XOR
      elementsInLeftChildBV[i]=xor(elementsInLeftChildBV[i],true) #updating needs to occur here before we copy the array (in case it is a better split than what we have seen so far)
      #move class from left to right side
      sumnl-=sumnumerator[i]
      sumdl-=sumdenominator[i]
      sumwl-=sumweight[i]
    else
    #update elementsInLeftChildBV, i.e. toggle the value of the ith component, $=XOR
      elementsInLeftChildBV[i]=xor(elementsInLeftChildBV[i],true)
        #move class from right to left side
      sumnl+=sumnumerator[i]
      sumdl+=sumdenominator[i]
      sumwl+=sumweight[i]
    end
    #save mean values (or should we save the frequency?)
    weightsl[counter]=sumwl
    weightsr[counter]=weighttot-sumwl
    meansl[counter]=sumnl/sumdl
    meansr[counter]=(numtot-sumnl)/(denomtot-sumdl)
    counter+=1
  end
  #reset bitvector
  fill!(elementsInLeftChildBV,false) #indicates which classes belong to right child

  #deviance vectors
  #the lenght of each vector is equal to the possible number of splits we consider!
  deviancel=zeros(subs_size) #deviance (frequency or possibly claim count) of the left node in case we go for split 'number' i
  deviancer=zeros(subs_size)

  lo=one(eltype(features.parent.refs)) #UInt8(1)
  ooo=one(lo)-lo

  val=-Inf
  chosen_sumwl=NaN

  #second pass over data
  #calculate deviances for each possible split
  counter=1
  for i in subs
      #update elementsInLeftChildBV, i.e. toggle the value of the ith component, $=XOR
      elementsInLeftChildBV[i]=xor(elementsInLeftChildBV[i],true)
      #calc deviances for the given split (which is defined by elementsInLeftChildBV)
      sumwl=weightsl[counter]
      #we can skip the calculation of the deviance, if we know that the leaves will be "too small"
      if (sumwl>minweight)&&(weighttot_minw>sumwl)        
        deviancel,deviancer=get_poisson_deviances(meansl[counter],meansr[counter],lo,ooo,features,numerator,denominator,weight,elementsInLeftChildBV)
      #end
      #if (sumwl>minweight)&&(weighttot_minw>sumwl) #do we have enough exposure? is the split valid?        
        valnew = -(deviancel+deviancer) #abs(sumnl/sumdl-(numtot-sumnl)/(denomtot-sumdl))
        if valnew>val
          val=valnew
          chosen_subset_bitarray=copy(elementsInLeftChildBV)
          chosen_sumwl=sumwl
        end
      end

    counter+=1
  end

    if isfinite(val)
      #warning: if the labellist is not sorted 1:n we need to be careful here!
      chosen_subset=labellist[chosen_subset_bitarray]
      else
      chosen_subset=Array{UInt}(0)
    end
    return val,chosen_subset,chosen_sumwl,weighttot-chosen_sumwl
end


#function calculateSplitValue(a::DifferenxxxxceSplit,fname::Symbol,number_of_num_features::Int,labellist::Vector{T},sumnumerator::Array{Float64,1},sumdenominator::Array{Float64,1},sumweight::Array{Float64,1},countlistfloat::Array{Float64,1},minweight::Float64,subs::DTSubsets,feature_column_id::Int64) where T<:Unsigned
function calculateSplitValue(a::PoissonDevianceSplit,fname::Symbol,number_of_num_features::Int,labellist::Vector{T},sumnumerator::Array{Float64,1},sumdenominator::Array{Float64,1},sumweight::Array{Float64,1},countlistfloat::Array{Float64,1},minweight::Float64,subs::DTSubsets,numerator::Array{Float64},denominator::Array{Float64},weight::Array{Float64},features,feature_column_id::Int64)  where T<:Unsigned
  error("need to add fname and the other unused argument")
  #here randomweight>0
#warn("xlogy will fail terribly if one of the arguments is zero! we need to amend the code to handle that case properly!")
  #here randomweight==0
  #for subsets, exhaustive search with flipping members (gray code) or "increasing" subset search ({1}, {1,2}, {1,2,3}, .... {1,2,3, ....., n-1,2})
  #all input lists (labellist,sumnumerator,sumdenominator,sumweight,countlistfloat) need to be sorted in the same manner

  #convention, in the beginning everything is on the right side
  elementsInLeftChildBV=BitVector(length(labellist));
  fill!(elementsInLeftChildBV,false) #indicates which classes belong to right child
  #chosen_subset_bitarray=BitVector(0)
  numtot=sum(sumnumerator)
  denomtot=sum(sumdenominator)
  weighttot=sum(sumweight)
  weighttot_minw=weighttot-minweight

  #step 1: find out how 'big' subs is
  subs_size=length(subs)
  #@assert subs_size<255 #todo tbd: we can remove this later on

  #pass on over the data:
  #get mean of left chlid and right child for each possible split
  #warn("find out if we need/want freq or claimcount here!")
  meansl=zeros(subs_size) #mean (frequency or possibly claim count) of the left node in case we go for split 'number' i
  meansr=zeros(subs_size) #same for right child
  #same for weights
  weightsl=zeros(subs_size)
  weightsr=zeros(subs_size)

  sumnl=sumwl=sumdl=0.0
  this_splitlist=Array{Splitdef}(0)

  counter=1
  for i in subs
    #@show "vv=$(i)" #should be the element which switches
    #i is an index, it indicates which element of labellist will flip sides for the next calculation
    if elementsInLeftChildBV[i]
    #update elementsInLeftChildBV, i.e. toggle the value of the ith component, $=XOR
      elementsInLeftChildBV[i]=xor(elementsInLeftChildBV[i],true) #updating needs to occur here before we copy the array (in case it is a better split than what we have seen so far)
      #move class from left to right side
      sumnl-=sumnumerator[i]
      sumdl-=sumdenominator[i]
      sumwl-=sumweight[i]
    else
    #update elementsInLeftChildBV, i.e. toggle the value of the ith component, $=XOR
      elementsInLeftChildBV[i]=xor(elementsInLeftChildBV[i],true)
        #move class from right to left side
      sumnl+=sumnumerator[i]
      sumdl+=sumdenominator[i]
      sumwl+=sumweight[i]
    end
    #save mean values (or should we save the frequency?)
    weightsl[counter]=sumwl # this is the exposure in an claim frequency model
    weightsr[counter]=weighttot-sumwl
    meansl[counter]=sumnl/sumdl #in a claim frequency model denominator is the exposure (i.e. euqual to weight)
    meansr[counter]=(numtot-sumnl)/(denomtot-sumdl)
    counter+=1
  end
  #reset bitvector
  fill!(elementsInLeftChildBV,false) #indicates which classes belong to right child

  #deviance vectors
  #the lenght of each vector is equal to the possible number of splits we consider!
  deviancel=zeros(subs_size) #deviance (frequency or possibly claim count) of the left node in case we go for split 'number' i
  deviancer=zeros(subs_size)

  lo=one(eltype(features.parent.refs)) #UInt8(1)
  ooo=one(lo)-lo

  #second pass over data
  #calculate deviances for each possible split
  counter=1
  for i in subs
      #update elementsInLeftChildBV, i.e. toggle the value of the ith component, $=XOR
      elementsInLeftChildBV[i]=xor(elementsInLeftChildBV[i],true)
      #calc deviances for the given split (which is defined by elementsInLeftChildBV)
      sumwl=weightsl[counter]
      #we can skip the calculation of the deviance, if we know that the leaves will be "too small"
      if (sumwl>minweight)&&(weighttot_minw>sumwl)        
        deviancel,deviancer=get_poisson_deviances(meansl[counter],meansr[counter],lo,ooo,features,numerator,denominator,weight,elementsInLeftChildBV)
      #end
      #if (sumwl>minweight)&&(weighttot_minw>sumwl) #do we have enough exposure? is the split valid?        				
        valnew = -(deviancel+deviancer) #abs(sumnl/sumdl-(numtot-sumnl)/(denomtot-sumdl))
		feature_column_id2 = feature_column_id < 0 ? abs(feature_column_id)+number_of_num_features : feature_column_id
        push!(this_splitlist,Splitdef(feature_column_id,feature_column_id2,fname,labellist[elementsInLeftChildBV],valnew,sumwl,weighttot-sumwl))
      end

    counter+=1
  end

  return this_splitlist
end

function get_poisson_deviances(current_meanl::Float64,current_meanr::Float64,lo,ooo,f,numerator::Array{Float64,1},denominator::Array{Float64,1},weight::Array{Float64,1},elementsInLeftChildBV)
	#consider the reacfin paper or the master thesis of 'axa'
	#this is the core function of the modelling process
	#besides copying of the data, the vast majority of time is spent in here!
	#most of the time is spent here, if we can improve the for loop below, that would improve performance greatly!
	#one possibility would be to introduce parallelization here (which is not straightforward, I think....)
	
	#dr and dl are 'reused' by each iteration
	dr=0.0
	dl=0.0
	#note: inbounds increases efficiency here (about a factor of 2), however if the bounds are violated something nasty might happen (quote: If the subscripts are ever out of bounds, you may suffer crashes or silent corruption.)
	#warn("need to add inbounds here!")
	#warn("need to treat ni == 0.0 specially....?")
	#warn("negative ni are impossible for poisson! -> DomainError")
	for count in 1:length(f)	
		@inbounds idx=f.parent.refs[count] + ooo		
		@inbounds ni = numerator[count]
		@inbounds wi = weight[count]
		if elementsInLeftChildBV[idx]
			dl += (xlogy(ni, ni / (wi*current_meanl)) - (ni - wi*current_meanl))
		else
			dr += (xlogy(ni, ni / (wi*current_meanr)) - (ni - wi*current_meanr))
		end
	end
	return dl,dr
end
