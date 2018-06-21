function calculateSplitValue(a::ROptMinRLostPctSplit,labellist::Array{UInt8,1},intSumrklost::Array{Int,1},sumweight::Array{Float64,1},countlistfloat::Array{Float64,1},minweight::Float64,subs::DTSubsets)
#here randomweight==0
#for subsets, exhaustive search with flipping members (gray code) or "increasing" subset search ({1}, {1,2}, {1,2,3}, .... {1,2,3, ....., n-1,2})
#all input lists (labellist,sumnumerator,sumdenominator,sumweight,countlistfloat) need to be sorted in the same manner
#convention, in the beginning everything is on the right side
elementsInLeftChildBV=BitVector(length(labellist));fill!(elementsInLeftChildBV,false) #indicates which classes belong to right child
chosen_subset_bitarray=BitVector(0)
rklosttot=float(sum(intSumrklost)) 
weighttot=sum(sumweight)
weighttot_minw=weighttot-minweight

sumrklostl=0
sumwl=0.0
val=-Inf
chosen_sumwl=NaN
  for i in subs
  #@show "vv=$(i)" #should be the element which switches
  #i is an index, it indicates which element of labellist will flip sides for the next calculation
  if elementsInLeftChildBV[i]
  #update elementsInLeftChildBV, i.e. toggle the value of the ith component, $=XOR
  elementsInLeftChildBV[i]=elementsInLeftChildBV[i]$true #updating needs to occur here before we copy the array (in case it is a better split than what we have seen so far)
    #move class from left to right side    
	sumrklostl-=intSumrklost[i]
	sumwl-=sumweight[i]    
  else
  #update elementsInLeftChildBV, i.e. toggle the value of the ith component, $=XOR
  elementsInLeftChildBV[i]=elementsInLeftChildBV[i]$true
      #move class from right to left side    
    sumrklostl+=intSumrklost[i]
	sumwl+=sumweight[i]    	
  end
    if (sumwl>minweight)&&(weighttot_minw>sumwl) #do we have enough exposure? is the split valid?
	  #vold=valnew
	  sumrklostlFLOAT=float(sumrklostl)
	  valnew=min(sumrklostlFLOAT/sumwl,(rklosttot-sumrklostlFLOAT)/(weighttot-sumwl)) #abs(sumnl/sumdl-(numtot-sumnl)/(denomtot-sumdl))	  
      #if vold==valnew;@show i,valnew,sumnl,sumdl;end;
	  if valnew>val
        val=valnew
        chosen_subset_bitarray=copy(elementsInLeftChildBV)
        chosen_sumwl=sumwl
      end
    end

end
if isfinite(val)
  #warning: if the labellist is not sorted 1:n we need to be careful here!
  chosen_subset=labellist[chosen_subset_bitarray]
else
  chosen_subset=Array(UInt,0)
end
return val,chosen_subset,chosen_sumwl,weighttot-chosen_sumwl
end


function calculateSplitValue(a::ROptMinRLostPctSplit,labellist::Array{UInt8,1},intSumrklost::Array{Int,1},sumweight::Array{Float64,1},countlistfloat::Array{Float64,1},minweight::Float64,subs::DTSubsets,feature_column_id::Int)
#here randomweight>0
#for subsets, exhaustive search with flipping members (gray code) or "increasing" subset search ({1}, {1,2}, {1,2,3}, .... {1,2,3, ....., n-1,2})
#all input lists (labellist,sumnumerator,sumdenominator,sumweight,countlistfloat) need to be sorted in the same manner
#convention, in the beginning everything is on the right side
elementsInLeftChildBV=BitVector(length(labellist));fill!(elementsInLeftChildBV,false) #indicates which classes belong to right child
rklosttot=float(sum(intSumrklost)) 
weighttot=sum(sumweight)
weighttot_minw=weighttot-minweight
sumrklostl=0
sumwl=0.0
val=-Inf
weighttot=sum(sumweight)
weighttot_minw=weighttot-minweight
this_splitlist=Array(Splitdef,0)

for i in subs
    #i is an index, it indicates which element of labellist will flip sides for the next calculation
 if elementsInLeftChildBV[i]
  #update elementsInLeftChildBV, i.e. toggle the value of the ith component, $=XOR
  elementsInLeftChildBV[i]=elementsInLeftChildBV[i]$true #updating needs to occur here before we copy the array (in case it is a better split than what we have seen so far)
    #move class from left to right side
    sumrklostl-=intSumrklost[i]
	sumwl-=sumweight[i]  
  else
  #update elementsInLeftChildBV, i.e. toggle the value of the ith component, $=XOR
  elementsInLeftChildBV[i]=elementsInLeftChildBV[i]$true
      #move class from right to left side    
    sumrklostl+=intSumrklost[i]	
	sumwl+=sumweight[i]
  end
    if (sumwl>minweight)&&(weighttot_minw>sumwl) #do we have enough exposure? is the split valid?
	  sumrklostlFLOAT=float(sumrklostl)
	  wr=weighttot-sumwl
	  valnew=min(sumrklostlFLOAT/sumwl,(rklosttot-sumrklostlFLOAT)/wr) #abs(sumnl/sumdl-(numtot-sumnl)/(denomtot-sumdl))	  	  
	  push!(this_splitlist,Splitdef(feature_column_id,copy(labellist[elementsInLeftChildBV]),valnew,sumwl,wr))
    end
end
return this_splitlist
end
