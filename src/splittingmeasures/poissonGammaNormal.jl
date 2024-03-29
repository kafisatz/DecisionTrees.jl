function calculateSplitValue(a::PGN, fname::Symbol, number_of_num_features::Int, labellist::Vector{T}, sumnumerator::Array{Float64,1}, sumdenominator::Array{Float64,1}, sumweight::Array{Float64,1}, countlistfloat::Array{Float64,1}, minweight::Float64, subs::DTSubsets, numerator::Array{Float64}, denominator::Array{Float64}, pointwiseRatio::Array{Float64,1}, weight::Array{Float64}, features) where {T <: Unsigned,PGN <: PoissonOrGamma}
    # here randomweight==0
    # for subsets, exhaustive search with flipping members (gray code) or "increasing" subset search ({1}, {1,2}, {1,2,3}, .... {1,2,3, ....., n-1,2})
    # all input lists (labellist,sumnumerator,sumdenominator,sumweight,countlistfloat) need to be sorted in the same manner
  
    # convention, in the beginning everything is on the right side
      elementsInLeftChildBV = BitVector(undef, length(labellist));
      fill!(elementsInLeftChildBV, false) # indicates which classes belong to right child
      chosen_subset_bitarray = BitVector(undef, 0)
      numtot = sum(sumnumerator)
      denomtot = sum(sumdenominator)
      weighttot = sum(sumweight)
      weighttot_minw = weighttot - minweight
  
    # step 1: find out how 'big' subs is
      subs_size = length(subs)
    # pass on over the data:
    # get mean of left child and right child for each possible split
      meansl = zeros(subs_size) # mean (frequency or possibly claim count) of the left node in case we go for split 'number' i
      meansr = zeros(subs_size) # same for right child
    # same for weights
      weightsl = zeros(subs_size)
      weightsr = zeros(subs_size)
  
      sumnl = sumwl = sumdl = 0.0
  
      counter = 1
      for i in subs
      # @show "vv=$(i)" #should be the element which switches
      # i is an index, it indicates which element of labellist will flip sides for the next calculation
          @inbounds if elementsInLeftChildBV[i]
      # update elementsInLeftChildBV, i.e. toggle the value of the ith component, $=XOR
              @inbounds elementsInLeftChildBV[i] = xor(elementsInLeftChildBV[i], true) # updating needs to occur here before we copy the array (in case it is a better split than what we have seen so far)
        # move class from left to right side
              @inbounds sumnl -= sumnumerator[i]
              @inbounds sumdl -= sumdenominator[i]
              @inbounds sumwl -= sumweight[i]
          else
      # update elementsInLeftChildBV, i.e. toggle the value of the ith component, $=XOR
              @inbounds elementsInLeftChildBV[i] = xor(elementsInLeftChildBV[i], true)
          # move class from right to left side
              @inbounds sumnl += sumnumerator[i]
              @inbounds sumdl += sumdenominator[i]
              @inbounds sumwl += sumweight[i]
          end
      # save mean values (or should we save the frequency?)
          @inbounds weightsl[counter] = sumwl # this is the exposure in an claim frequency model
          @inbounds weightsr[counter] = weighttot - sumwl
          @inbounds meansl[counter] = sumnl / sumdl # in a claim frequency model denominator is the exposure (i.e. equal to weight)
          @inbounds meansr[counter] = (numtot - sumnl) / (denomtot - sumdl)
          counter += 1
      end
    # reset bitvector
      fill!(elementsInLeftChildBV, false) # indicates which classes belong to right child
  
    # deviance vectors
    # the lenght of each vector is equal to the possible number of splits we consider!
      deviancel = zeros(subs_size) # deviance (frequency or possibly claim count) of the left node in case we go for split 'number' i
      deviancer = zeros(subs_size)
  
      lo = one(eltype(features.parent.refs)) # UInt8(1)
      ooo = one(lo) - lo
  
      val = -Inf
      chosen_sumwl = NaN
  
    # second pass over data
    # calculate deviances for each possible split
      counter = 1  
      for i in subs
        # update elementsInLeftChildBV, i.e. toggle the value of the ith component, $=XOR
          @inbounds elementsInLeftChildBV[i] = xor(elementsInLeftChildBV[i], true)
        # calc deviances for the given split (which is defined by elementsInLeftChildBV)
          @inbounds sumwl = weightsl[counter]      
        # we can skip the calculation of the deviance, if we know that the leaves will be "too small"
          if (sumwl > minweight) && (weighttot_minw > sumwl)        
              @inbounds deviancel, deviancer = get_deviances(a, meansl[counter], meansr[counter], lo, ooo, features, numerator, denominator, pointwiseRatio, weight, elementsInLeftChildBV)
              valnew = -(deviancel + deviancer) # abs(sumnl/sumdl-(numtot-sumnl)/(denomtot-sumdl))        
              if valnew > val
                  val = valnew
                  chosen_subset_bitarray = copy(elementsInLeftChildBV)
                  chosen_sumwl = sumwl
              end
          end
  
          counter += 1
      end
  
      if isfinite(val)
        # warning: if the labellist is not sorted 1:n we need to be careful here!
          chosen_subset = labellist[chosen_subset_bitarray]
      else
          chosen_subset = Array{T}(undef, 0)
      end
      return val, chosen_subset, chosen_sumwl, weighttot - chosen_sumwl
  end
  
  function calculateSplitValue(a::PGN, fname::Symbol, number_of_num_features::Int, labellist::Vector{T}, sumnumerator::Array{Float64,1}, sumdenominator::Array{Float64,1}, sumweight::Array{Float64,1}, countlistfloat::Array{Float64,1}, minweight::Float64, subs::DTSubsets, numerator::Array{Float64}, denominator::Array{Float64}, pointwiseRatio::Array{Float64,1}, weight::Array{Float64}, features, feature_column_id::Int) where {T <: Unsigned,PGN <: PoissonOrGamma}    
  # here randomweight>0

  # for subsets, exhaustive search with flipping members (gray code) or "increasing" subset search ({1}, {1,2}, {1,2,3}, .... {1,2,3, ....., n-1,2})
  # all input lists (labellist,sumnumerator,sumdenominator,sumweight,countlistfloat) need to be sorted in the same manner
  # convention, in the beginning everything is on the right side
    elementsInLeftChildBV = BitVector(undef, length(labellist));
    fill!(elementsInLeftChildBV, false) # indicates which classes belong to right child
  # chosen_subset_bitarray=BitVector(undef,0)
    numtot = sum(sumnumerator)
    denomtot = sum(sumdenominator)
    weighttot = sum(sumweight)
    weighttot_minw = weighttot - minweight

  # step 1: find out how 'big' subs is
    subs_size = length(subs)
  # pass on over the data:
  # get mean of left child and right child for each possible split
    meansl = zeros(subs_size) # mean (frequency or possibly claim count) of the left node in case we go for split 'number' i
    meansr = zeros(subs_size) # same for right child
  # same for weights
    weightsl = zeros(subs_size)
    weightsr = zeros(subs_size)

    sumnl = sumwl = sumdl = 0.0
    this_splitlist = Array{Splitdef{T}}(undef, 0)

    counter = 1
    for i in subs
    # @show "vv=$(i)" #should be the element which switches
    # i is an index, it indicates which element of labellist will flip sides for the next calculation
        @inbounds if elementsInLeftChildBV[i]
    # update elementsInLeftChildBV, i.e. toggle the value of the ith component, $=XOR
            @inbounds elementsInLeftChildBV[i] = xor(elementsInLeftChildBV[i], true) # updating needs to occur here before we copy the array (in case it is a better split than what we have seen so far)
      # move class from left to right side
            @inbounds sumnl -= sumnumerator[i]
            @inbounds sumdl -= sumdenominator[i]
            @inbounds sumwl -= sumweight[i]
        else
    # update elementsInLeftChildBV, i.e. toggle the value of the ith component, $=XOR
            @inbounds elementsInLeftChildBV[i] = xor(elementsInLeftChildBV[i], true)
        # move class from right to left side
            @inbounds sumnl += sumnumerator[i]
            @inbounds sumdl += sumdenominator[i]
            @inbounds sumwl += sumweight[i]
        end
    # save mean values (or should we save the frequency?)
        @inbounds weightsl[counter] = sumwl # this is the exposure in an claim frequency model
        @inbounds weightsr[counter] = weighttot - sumwl
        @inbounds meansl[counter] = sumnl / sumdl # in a claim frequency model denominator is the exposure (i.e. equal to weight)
        @inbounds meansr[counter] = (numtot - sumnl) / (denomtot - sumdl)
        counter += 1
    end
  # reset bitvector
    fill!(elementsInLeftChildBV, false) # indicates which classes belong to right child

  # deviance vectors
  # the lenght of each vector is equal to the possible number of splits we consider!
    deviancel = zeros(subs_size) # deviance (frequency or possibly claim count) of the left node in case we go for split 'number' i
    deviancer = zeros(subs_size)

    lo = one(eltype(features.parent.refs)) # UInt8(1)
    ooo = one(lo) - lo

  # second pass over data
  # calculate deviances for each possible split
    counter = 1  
    for i in subs
      # update elementsInLeftChildBV, i.e. toggle the value of the ith component, $=XOR
        @inbounds elementsInLeftChildBV[i] = xor(elementsInLeftChildBV[i], true)
      # calc deviances for the given split (which is defined by elementsInLeftChildBV)
        @inbounds sumwl = weightsl[counter]      
      # we can skip the calculation of the deviance, if we know that the leaves will be "too small"
        if (sumwl > minweight) && (weighttot_minw > sumwl)        
            @inbounds deviancel, deviancer = get_deviances(a, meansl[counter], meansr[counter], lo, ooo, features, numerator, denominator, pointwiseRatio, weight, elementsInLeftChildBV)
     		valnew = -(deviancel + deviancer) # abs(sumnl/sumdl-(numtot-sumnl)/(denomtot-sumdl))        
			feature_column_id2 = feature_column_id < 0 ? abs(feature_column_id) + number_of_num_features : feature_column_id
            push!(this_splitlist, Splitdef{T}(feature_column_id, feature_column_id2, fname, labellist[elementsInLeftChildBV], valnew, sumwl, weighttot - sumwl))
        end

        counter += 1
    end

    return this_splitlist
end

function get_deviances(a::PoissonDevianceSplit, current_meanl::Float64, current_meanr::Float64, lo, ooo, f, numerator::Array{Float64,1}, denominator::Array{Float64,1}, pointwiseRatio::Array{Float64,1}, weight::Array{Float64,1}, elementsInLeftChildBV)
    # for the poisson deviance consider the derivative of the poisson loss (or google, the reacfin paper or the 'axa' master thesis)
	# this is the core function of the modelling process
	# besides copying of the data, the vast majority of time is spent in here!
	# most of the time is spent here, if we can improve the for loop below, that would improve performance greatly!
	# one possibility would be to introduce parallelization here (which is not straightforward, though
		
	dr = 0.0
	dl = 0.0
	# note: inbounds increases efficiency here (about a factor of 2), however if the bounds are violated something nasty might happen (quote: If the subscripts are ever out of bounds, you may suffer crashes or silent corruption.)
	for count in 1:length(f)	
		@inbounds idx = f.parent.refs[count] + ooo		
		@inbounds ni = numerator[count]
        @inbounds di = denominator[count]
		@inbounds eli = elementsInLeftChildBV[idx]       
        if eli
            if iszero(ni)
                dl += -(ni - di * current_meanl)
            else 
                dl += -(ni - di * current_meanl) - ni * log(di * current_meanl / ni)
            end
        else            
            if iszero(ni)
                dr += -(ni - di * current_meanr)
            else 
                dr += -(ni - di * current_meanr) - ni * log(di * current_meanr / ni)
            end        
        end
	end
	return dl, dr
end


function get_deviances(a::GammaDevianceSplit, current_meanl::Float64, current_meanr::Float64, lo, ooo, f, numerator::Array{Float64,1}, denominator::Array{Float64,1}, pointwiseRatio::Array{Float64,1}, weight::Array{Float64,1}, elementsInLeftChildBV)    
	dr = 0.0
	dl = 0.0
	# note: inbounds increases efficiency here (about a factor of 2), however if the bounds are violated something nasty might happen (quote: If the subscripts are ever out of bounds, you may suffer crashes or silent corruption.)
	for count in 1:length(f)	
		@inbounds idx = f.parent.refs[count] + ooo		
		@inbounds ni = numerator[count]
        @inbounds di = denominator[count]
        @inbounds eli = elementsInLeftChildBV[idx]

        if eli
            if iszero(ni)
                dl += - (ni - current_meanl * di) / (di * current_meanl)
            else 
                dl += -(ni - di * current_meanl) / (di * current_meanl) + log(di * current_meanl / ni)
            end
        else            
            if iszero(ni)
                dr += - (ni - current_meanr * di) / (di * current_meanr)
            else 
                dr += -(ni - di * current_meanr) / (di * current_meanr) + log(di * current_meanr / ni)
            end        
        end        
        #= 
        #I think this is slightly slower
        diTimesMean=ifelse(eli,di*current_meanl,di*current_meanr)        
        if iszero(ni)
          addition = - (ni - diTimesMean)/diTimesMean 
        else
          addition = - (ni - diTimesMean)/diTimesMean + log(diTimesMean / ni)
        end        
        if eli            
            dl += addition
        else            
            dr += addition
        end =#
	end
	return dl, dr
end

#= 
#old functions
function get_deviances(a::mseSplit,current_meanl::Float64,current_meanr::Float64,lo,ooo,f,numerator::Array{Float64,1},denominator::Array{Float64,1},pointwiseRatio::Array{Float64,1},weight::Array{Float64,1},elementsInLeftChildBV)    
	dr=0.0
	dl=0.0
	#note: inbounds increases efficiency here (about a factor of 2), however if the bounds are violated something nasty might happen (quote: If the subscripts are ever out of bounds, you may suffer crashes or silent corruption.)
	for count in 1:length(f)	
		@inbounds idx = f.parent.refs[count] + ooo		
		#@inbounds ni = numerator[count]
        #@inbounds di = denominator[count]
		@inbounds ratioi = pointwiseRatio[count]
        @inbounds eli=elementsInLeftChildBV[idx]        
        if eli          
            dl += (ratioi - current_meanl)^2
        else            
            dr += (ratioi - current_meanr)^2
      end                
	end
	return dl,dr
end


function get_deviances(a::NormalDevianceDifferenceToMeanFitWEIGHTEDSplit,current_meanl::Float64,current_meanr::Float64,lo,ooo,f,numerator::Array{Float64,1},denominator::Array{Float64,1},weight::Array{Float64,1},elementsInLeftChildBV)    
	dr=0.0
	dl=0.0
	#note: inbounds increases efficiency here (about a factor of 2), however if the bounds are violated something nasty might happen (quote: If the subscripts are ever out of bounds, you may suffer crashes or silent corruption.)
	for count in 1:length(f)	
		@inbounds idx = f.parent.refs[count] + ooo		
		#@inbounds ni = numerator[count]
        @inbounds di = denominator[count]
		#@inbounds wi = weight[count]
		@inbounds ratioi = pointwiseRatio[count]
        @inbounds eli = elementsInLeftChildBV[idx]     
        if eli          
            dl += di*(ratioi - current_meanl)^2
        else            
            dr += di*(ratioi - current_meanr)^2
      end                
	end
	return dl,dr
end =#