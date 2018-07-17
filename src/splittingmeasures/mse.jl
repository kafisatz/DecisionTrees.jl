#=
    this splitting criterion considers the  ratio = sum(numerator)/sum(denominator) (for each split§)
    we minimize the (ratio.-meanLeft)^2 of both children
=#


function calculateSplitValue(a::sseSplit,fname::Symbol,number_of_char_features::Int,labellist::Vector{T},sumnumerator::Array{Float64,1},sumdenominator::Array{Float64,1},sumweight::Array{Float64,1},countlistfloat::Array{Float64,1},minweight::Float64,subs::DTSubsets,moments_per_pdaclass) where T<:Unsigned
    #here randomweight==0
    #for subsets, exhaustive search with flipping members (gray code) or "increasing" subset search ({1}, {1,2}, {1,2,3}, .... {1,2,3, ....., n-1,2})
    #all input lists (labellist,sumnumerator,sumdenominator,sumweight,countlistfloat) need to be sorted in the same manner
    #convention: in the beginning everything is on the right side
    elementsInLeftChildBV=BitVector(undef,length(labellist));fill!(elementsInLeftChildBV,false) #indicates which classes belong to right child
    chosen_subset_bitarray=BitVector(undef,0)
    
    totalMeanAndVariance=copy(emptyCustomVariance)
    for idx=1:length(moments_per_pdaclass)
        @inbounds merge!(totalMeanAndVariance,moments_per_pdaclass[idx])
    end
    
    weighttot=sum(sumweight)
    weighttot_minw=weighttot-minweight
    
    momentsl=copy(emptyCustomVariance) 
    #momentsr needs to be the 'total' at the start of this.
    momentsr=totalMeanAndVariance 
    
    sumwl=0.0
    varl=varr=val=-Inf
    chosen_sumwl=NaN
      for i in subs
      #@show "vv=$(i)" #should be the element which switches
      #i is an index, it indicates which element of labellist will flip sides for the next calculation
      @inbounds switching_class=moments_per_pdaclass[i]
      @inbounds if elementsInLeftChildBV[i]
        #update elementsInLeftChildBV, i.e. toggle the value of the ith component, $=XOR
        @inbounds  elementsInLeftChildBV[i]=xor(elementsInLeftChildBV[i],true) #updating needs to occur here before we copy the array (in case it is a better split than what we have seen so far)
        #move class from left to right side
        unmerge!(momentsl,switching_class) #switching_class is 'removed' from the left variance/mean statistcs momentsl
        merge!(momentsr,switching_class) #switching_class is 'added' to the right variance/mean statistcs momentsl
        #sumnl-=sumnumerator[i]
        #sumdl-=sumdenominator[i]
        @inbounds  sumwl-=sumweight[i]
      else
          #update elementsInLeftChildBV, i.e. toggle the value of the ith component, $=XOR
          @inbounds  elementsInLeftChildBV[i]=xor(elementsInLeftChildBV[i],true)
          #move class from right to left side
          merge!(momentsl,switching_class)
          unmerge!(momentsr,switching_class)
          #sumnl+=sumnumerator[i]
            #sumdl+=sumdenominator[i]
            @inbounds  sumwl+=sumweight[i]
      end
        if (sumwl>minweight)&&(weighttot_minw>sumwl) #do we have enough exposure? is the split valid?
          #vold=valnew
          varl=momentsl.m2
          varr=momentsr.m2
          #@show fname,sumwl,weighttot-sumwl
          valnew = -(varl+varr)
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
      chosen_subset=Array{UInt}(undef,0)
    end
    return val,chosen_subset,chosen_sumwl,weighttot-chosen_sumwl
    end

    
function calculateSplitValue(a::mseSplit,fname::Symbol,number_of_char_features::Int,labellist::Vector{T},sumnumerator::Array{Float64,1},sumdenominator::Array{Float64,1},sumweight::Array{Float64,1},countlistfloat::Array{Float64,1},minweight::Float64,subs::DTSubsets,moments_per_pdaclass) where T<:Unsigned
#here randomweight==0
#for subsets, exhaustive search with flipping members (gray code) or "increasing" subset search ({1}, {1,2}, {1,2,3}, .... {1,2,3, ....., n-1,2})
#all input lists (labellist,sumnumerator,sumdenominator,sumweight,countlistfloat) need to be sorted in the same manner
#convention: in the beginning everything is on the right side
elementsInLeftChildBV=BitVector(undef,length(labellist));fill!(elementsInLeftChildBV,false) #indicates which classes belong to right child
chosen_subset_bitarray=BitVector(undef,0)

totalMeanAndVariance=copy(emptyCustomVariance)
for idx=1:length(moments_per_pdaclass)
    @inbounds merge!(totalMeanAndVariance,moments_per_pdaclass[idx])
end

weighttot=sum(sumweight)
weighttot_minw=weighttot-minweight

momentsl=copy(emptyCustomVariance) 
#momentsr needs to be the 'total' at the start of this.
momentsr=totalMeanAndVariance 

sumwl=0.0
varl=varr=val=-Inf
chosen_sumwl=NaN
  for i in subs
  #@show "vv=$(i)" #should be the element which switches
  #i is an index, it indicates which element of labellist will flip sides for the next calculation
  @inbounds switching_class=moments_per_pdaclass[i]
  @inbounds if elementsInLeftChildBV[i]
    #update elementsInLeftChildBV, i.e. toggle the value of the ith component, $=XOR
    @inbounds  elementsInLeftChildBV[i]=xor(elementsInLeftChildBV[i],true) #updating needs to occur here before we copy the array (in case it is a better split than what we have seen so far)
    #move class from left to right side
    unmerge!(momentsl,switching_class) #switching_class is 'removed' from the left variance/mean statistcs momentsl
    merge!(momentsr,switching_class) #switching_class is 'added' to the right variance/mean statistcs momentsl
    #sumnl-=sumnumerator[i]
    #sumdl-=sumdenominator[i]
    @inbounds  sumwl-=sumweight[i]
  else
      #update elementsInLeftChildBV, i.e. toggle the value of the ith component, $=XOR
      @inbounds  elementsInLeftChildBV[i]=xor(elementsInLeftChildBV[i],true)
      #move class from right to left side
      merge!(momentsl,switching_class)
      unmerge!(momentsr,switching_class)
      #sumnl+=sumnumerator[i]
    	#sumdl+=sumdenominator[i]
    	@inbounds  sumwl+=sumweight[i]
  end
    if (sumwl>minweight)&&(weighttot_minw>sumwl) #do we have enough exposure? is the split valid?
	  #vold=valnew
      varl=momentsl.var
      varr=momentsr.var
	  #@show fname,sumwl,weighttot-sumwl
      valnew = -(varl+varr)
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
  chosen_subset=Array{UInt}(undef,0)
end
return val,chosen_subset,chosen_sumwl,weighttot-chosen_sumwl
end


function calculateSplitValue(a::mseSplit,number_of_char_features::Int,labellist::Array{UInt8,1},sumnumerator::Array{Float64,1},sumdenominator::Array{Float64,1},sumweight::Array{Float64,1},countlistfloat::Array{Float64,1},minweight::Float64,subs::DTSubsets,feature_column_id::Int,moments_per_pdaclass)
error("221x this will not work yet.....")
#here randomweight>0
#for subsets, exhaustive search with flipping members (gray code) or "increasing" subset search ({1}, {1,2}, {1,2,3}, .... {1,2,3, ....., n-1,2})
#all input lists (labellist,sumnumerator,sumdenominator,sumweight,countlistfloat) need to be sorted in the same manner
#convention, in the beginning everything is on the right side
end

function aggregate_data_mse(f,numerator::Array{Float64,1},denominator::Array{Float64,1},weight::Array{Float64,1})
    #this is the core function of the modelling process
    #besides copying of the data, the vast majority of time is spent in here!
    #most of the time is spent here, if we can improve the for loop below, that would improve performance greatly!
    #one possibility would be to introduce parallelization here (which is not straightforward, I think....)
    lo=one(eltype(f.parent.refs))
    hi=convert(eltype(f.parent.refs),length(f.parent.pool)) 
    ooo=one(lo)-lo
    vecsize=hi+ooo
    cnt = zeros(Int, vecsize)
    sumnumerator = zeros(Float64, vecsize)
    sumdenominator = zeros(Float64, vecsize)
    sumweight= zeros(Float64, vecsize)
    moments_per_pdaclass=Array{CustomVariance}(undef,vecsize)
    for i=1:vecsize
        @inbounds moments_per_pdaclass[i]=copy(emptyCustomVariance)
    end
    @inbounds for count in f.indices[1]
        @inbounds idx=f.parent.refs[count] + ooo
        @inbounds cnt[idx] += 1
        @inbounds this_num = numerator[count]
        @inbounds sumnumerator[idx] += this_num
        @inbounds this_denom = denominator[count]
        @inbounds sumdenominator[idx] += this_denom
        @inbounds sumweight[idx] += weight[count]
        @inbounds addSinglePoint!(moments_per_pdaclass[idx],this_num,this_denom)
    end
  return cnt,sumnumerator,sumdenominator,sumweight,moments_per_pdaclass
end


mutable struct CustomVariance
    meanratio::Float64
    sn::Float64 #sum nominator
    sd::Float64 #sum denominator    
    n::Int #count
    m2::Float64 # var*(n-1)
    var::Float64 #weigthedVariance
    function CustomVariance(num,denom)
        n=length(num)
        @assert n==length(denom) # length(w)
        sn=sum(num)
        sd=sum(denom)        
        meanratio=sn/sd
        m2=zero(eltype(num))           
        for i=1:n
            @inbounds m2 += denom[i]*(num[i]/denom[i]-meanratio)^2
        end        
        if n>1
            var = m2 / (n-1)
        else 
            var=zero(eltype(num))
        end
        return new(meanratio,sn,sd,n,m2,var)
    end
    function CustomVariance()
        return new(0.0,0.0,0.0,0,0.0,0.0)
    end
end

function removeSinglePoint!(a::CustomVariance,num,denom)
    sn=a.sn-num
    sd=a.sd-denom
    meanratio=sn/sd
    #merge: m2 = a.m2 + b.m2 + 2*(a.sn*(a.meanratio-meanratio)+§*(num/denom-meanratio)) + a.sd*(meanratio^2-a.meanratio^2) + denom*(meanratio^2-num/denom^2)
    #unmerge:    
    a.m2 = a.m2 + 2*(num*(a.meanratio-num/denom) + sn*(a.meanratio-meanratio)) + denom*((num/denom)^2-a.meanratio^2) + sd*(meanratio^2-a.meanratio^2)
    a.sn=sn
    a.sd=sd
    a.n -= 1
    a.meanratio=meanratio
    a.var=a.m2/(a.n-1)
    return nothing 
end

function addSinglePoint!(a::CustomVariance,num,denom)
    sn=a.sn+num
    sd=a.sd+denom   
    a.n=a.n+1
    meanratio=sn/sd
    a.m2 = a.m2 + 2*(a.sn*(a.meanratio-meanratio)+num*(num/denom-meanratio)) + a.sd*(meanratio^2-a.meanratio^2) + denom*(meanratio^2-(num/denom)^2)
    a.sn=sn
    a.sd=sd
    a.meanratio=meanratio
    a.var=a.m2/(a.n-1)
    return nothing 
end

function Base.merge!(a::CustomVariance,b::CustomVariance)
    #merges b into a
    #i.e. the points form b are 'added' to a
    sn=a.sn+b.sn
    sd=a.sd+b.sd   
    n=a.n+b.n
    meanratio=sn/sd
    m2 = a.m2 + b.m2 + 2*(a.sn*(a.meanratio-meanratio)+b.sn*(b.meanratio-meanratio)) + a.sd*(meanratio^2-a.meanratio^2) + b.sd*(meanratio^2-b.meanratio^2)
    var = m2/(n-1)
    a.sn=sn
    a.sd=sd
    a.n=n
    a.meanratio=meanratio
    a.m2=m2
    a.var=var
    return nothing
end

function Base.merge(a::CustomVariance,b::CustomVariance)
    c=copy(a)
    merge!(c,b)
    return c
end 

"""
function unmerge!(a::CustomVariance,b::CustomVariance)
this function updates the mean and variance of a
The assumption is that the underlying data of b is a subset of a.
If this is not the case, the result of this function is most likely meaningless
That is to say "subtract" the series b from the series a and update the mean and variance of a accordingly
"""
function unmerge!(a::CustomVariance,b::CustomVariance)
    sn=a.sn-b.sn
    sd=a.sd-b.sd   
    n=a.n-b.n
    meanratio=sn/sd
    #merge: m2 = a.m2 + b.m2 + 2*(a.sn*(a.meanratio-meanratio)+b.sn*(b.meanratio-meanratio)) + a.sd*(meanratio^2-a.meanratio^2) + b.sd*(meanratio^2-b.meanratio^2)
    #unmerge:    
    m2 = a.m2 - b.m2 + 2*(b.sn*(a.meanratio-b.meanratio) + sn*(a.meanratio-meanratio)) + b.sd*(b.meanratio^2-a.meanratio^2) + sd*(meanratio^2-a.meanratio^2)
    var = m2/(n-1)
    a.sn=sn
    a.sd=sd
    a.n=n
    a.meanratio=meanratio
    a.m2=m2
    a.var=var    
    return nothing
end

function Base.copy(a::CustomVariance)
    b=CustomVariance()
    b.sn=a.sn
    b.sd=a.sd
    b.n=a.n
    b.meanratio=a.meanratio
    b.m2=a.m2
    b.var=a.var    
    return b
end