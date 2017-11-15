function calculateSplitValue(a::DifferenceSplit,fname::Symbol,number_of_num_features::Int,labellist::Array{UInt8,1},sumnumerator::Array{Float64,1},sumdenominator::Array{Float64,1},sumweight::Array{Float64,1},countlistfloat::Array{Float64,1},minweight::Float64,subs::DTSubsets)
#here randomweight==0
#for subsets, exhaustive search with flipping members (gray code) or "increasing" subset search ({1}, {1,2}, {1,2,3}, .... {1,2,3, ....., n-1,2})
#all input lists (labellist,sumnumerator,sumdenominator,sumweight,countlistfloat) need to be sorted in the same manner
#convention, in the beginning everything is on the right side
elementsInLeftChildBV=BitVector(length(labellist));fill!(elementsInLeftChildBV,false) #indicates which classes belong to right child
chosen_subset_bitarray=similar(elementsInLeftChildBV) #BitVector(0)
numtot=sum(sumnumerator)
denomtot=sum(sumdenominator)
weighttot=sum(sumweight)
weighttot_minw=weighttot-minweight

sumnl=sumdl=sumwl=0.0
val=-Inf
chosen_sumwl=NaN
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
    if (sumwl>minweight)&&(weighttot_minw>sumwl) #do we have enough exposure? is the split valid?
	  #vold=valnew
	  valnew=abs(sumnl/sumdl-(numtot-sumnl)/(denomtot-sumdl))	  
      #if vold==valnew;@show i,valnew,sumnl,sumdl;end;
	  if valnew>val
        val=valnew
        #chosen_subset_bitarray=copy(elementsInLeftChildBV)        
        chosen_subset_bitarray.=elementsInLeftChildBV
        chosen_sumwl=sumwl
      end
    end
end
if isfinite(val)
  #warning: if the labellist is not sorted 1:n we need to be careful here!
  chosen_subset=labellist[chosen_subset_bitarray]
else
  chosen_subset=Array{UInt}(0)
end
return val,chosen_subset,chosen_sumwl,weighttot-chosen_sumwl
end


function calculateSplitValue(a::DifferenceSplit,fname::Symbol,number_of_num_features::Int,labellist::Array{UInt8,1},sumnumerator::Array{Float64,1},sumdenominator::Array{Float64,1},sumweight::Array{Float64,1},countlistfloat::Array{Float64,1},minweight::Float64,subs::DTSubsets,feature_column_id::Int64)
#here randomweight>0
#for subsets, exhaustive search with flipping members (gray code) or "increasing" subset search ({1}, {1,2}, {1,2,3}, .... {1,2,3, ....., n-1,2})
#all input lists (labellist,sumnumerator,sumdenominator,sumweight,countlistfloat) need to be sorted in the same manner
#convention, in the beginning everything is on the right side
elementsInLeftChildBV=BitVector(length(labellist));fill!(elementsInLeftChildBV,false) #indicates which classes belong to right child
numtot=sum(sumnumerator)
denomtot=sum(sumdenominator)
weighttot=sum(sumweight)
weighttot_minw=weighttot-minweight
sumnl=sumdl=sumwl=0.0
this_splitlist=Array{Splitdef}(0)

for i in subs
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
    if (sumwl>minweight)&&(weighttot_minw>sumwl) #do we have enough exposure? is the split valid?
      valnew=abs(sumnl/sumdl-(numtot-sumnl)/(denomtot-sumdl))	  
      feature_column_id2 = feature_column_id < 0 ? abs(feature_column_id)+number_of_num_features : feature_column_id
	  push!(this_splitlist,Splitdef(feature_column_id,feature_column_id2,fname,labellist[elementsInLeftChildBV],valnew,sumwl,weighttot-sumwl))
    end
end
return this_splitlist
end

