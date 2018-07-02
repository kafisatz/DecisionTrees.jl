#calculateSplitValue(::DecisionTrees.NormalDevianceSplit, ::Symbol, ::Int64, ::Array{UInt8,1}, ::Array{Float64,1}, ::Array{Float64,1}, ::Array{Float64,1}, ::Array{Float64,1}, ::Float64, ::DecisionTrees.MyIncreasingSubsets, ::Array{OnlineStats.Series,1})
function calculateSplitValue(a::NormalDevianceSplit,fname::Symbol,number_of_char_features::Int,labellist::Vector{T},sumnumerator::Array{Float64,1},sumdenominator::Array{Float64,1},sumweight::Array{Float64,1},countlistfloat::Array{Float64,1},minweight::Float64,subs::DTSubsets,moments_per_pdaclass) where T<:Unsigned
#here randomweight==0
#for subsets, exhaustive search with flipping members (gray code) or "increasing" subset search ({1}, {1,2}, {1,2,3}, .... {1,2,3, ....., n-1,2})
#all input lists (labellist,sumnumerator,sumdenominator,sumweight,countlistfloat) need to be sorted in the same manner
#convention: in the beginning everything is on the right side
elementsInLeftChildBV=BitVector(undef,length(labellist));fill!(elementsInLeftChildBV,false) #indicates which classes belong to right child
chosen_subset_bitarray=BitVector(undef,0)

totalMeanAndVariance=copy(emptyMeanVarSeries)
for idx=1:length(moments_per_pdaclass)
    @inbounds  mi=moments_per_pdaclass[idx]
    @inbounds statsI=mi.stats[1]
    if statsI.n >0 #this is needed because of https://github.com/joshday/OnlineStats.jl/issues/125
        @inbounds merge!(totalMeanAndVariance,mi)
    end
end

weighttot=sum(sumweight)
weighttot_minw=weighttot-minweight

momentsl=copy(emptyMeanVarSeries) #Series(OnlineStats.Mean(),OnlineStats.Variance())
#momentsr needs to be the 'total' at the start of this.
momentsr=totalMeanAndVariance #copy(emptyMeanVarSeries) #Series(OnlineStats.Mean(),OnlineStats.Variance())

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
      @inbounds varl=momentsl.stats[2].σ2
      @inbounds varr=momentsr.stats[2].σ2
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


function calculateSplitValue(a::NormalDevianceSplit,number_of_char_features::Int,labellist::Array{UInt8,1},sumnumerator::Array{Float64,1},sumdenominator::Array{Float64,1},sumweight::Array{Float64,1},countlistfloat::Array{Float64,1},minweight::Float64,subs::DTSubsets,feature_column_id::Int,moments_per_pdaclass)
@warn("221 this will not work yet.....")
#here randomweight>0
#for subsets, exhaustive search with flipping members (gray code) or "increasing" subset search ({1}, {1,2}, {1,2,3}, .... {1,2,3, ....., n-1,2})
#all input lists (labellist,sumnumerator,sumdenominator,sumweight,countlistfloat) need to be sorted in the same manner
#convention, in the beginning everything is on the right side
elementsInLeftChildBV=BitVector(undef,length(labellist));fill!(elementsInLeftChildBV,false) #indicates which classes belong to right child

totalMeanAndVariance=copy(emptyMeanVarSeries)
for idx=1:length(moments_per_pdaclass)
    @inbounds  mi=moments_per_pdaclass[idx]
    @inbounds statsI=mi.stats[1]
    if statsI.n >0 #this is needed because of https://github.com/joshday/OnlineStats.jl/issues/125
        @inbounds merge!(totalMeanAndVariance,mi)
    end
end

weighttot=sum(sumweight)
weighttot_minw=weighttot-minweight

momentsl=copy(emptyMeanVarSeries) #Series(OnlineStats.Mean(),OnlineStats.Variance())
#momentsr needs to be the 'total' at the start of this.
momentsr=totalMeanAndVariance #copy(emptyMeanVarSeries) #Series(OnlineStats.Mean(),OnlineStats.Variance())

sumwl=0.0
this_splitlist=Array{Splitdef}(undef,0)

for i in subs
    #i is an index, it indicates which element of labellist will flip sides for the next calculation
 if elementsInLeftChildBV[i]
    #update elementsInLeftChildBV, i.e. toggle the value of the ith component, $=XOR
    @inbounds elementsInLeftChildBV[i]=xor(elementsInLeftChildBV[i],true) #updating needs to occur here before we copy the array (in case it is a better split than what we have seen so far)
    #move class from left to right side
    unmerge!(momentsl,switching_class)
    merge!(momentsr,switching_class)
	@inbounds   sumwl-=sumweight[i]
  else
    #update elementsInLeftChildBV, i.e. toggle the value of the ith component, $=XOR
    @inbounds elementsInLeftChildBV[i]=xor(elementsInLeftChildBV[i],true)
    #move class from right to left side
    merge!(momentsl,switching_class)
    unmerge!(momentsr,switching_class)
	@inbounds   sumwl+=sumweight[i]
  end
    if (sumwl>minweight)&&(weighttot_minw>sumwl) #do we have enough exposure? is the split valid?
        @inbounds varl=momentsl.stats[2].σ2
        @inbounds varr=momentsr.stats[2].σ2
	  	valnew=-(varl+varr)
        @inbounds push!(this_splitlist,Splitdef{T}(feature_column_id,feature_column_id2,fname,labellist[elementsInLeftChildBV],valnew,sumwl,weighttot-sumwl))
    end
end
return this_splitlist
end

function aggregate_data_normal_deviance(f,numerator::Array{Float64,1},denominator::Array{Float64,1},weight::Array{Float64,1})
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
    moments_per_pdaclass=Array{OnlineStats.Series}(undef,vecsize)
    for i=1:vecsize
        @inbounds moments_per_pdaclass[i]=copy(emptyMeanVarSeries)
    end
    @inbounds for count in f.indices[1]
        @inbounds idx=f.parent.refs[count] + ooo
        @inbounds cnt[idx] += 1
        @inbounds this_num = numerator[count]
        @inbounds sumnumerator[idx] += this_num
        @inbounds this_denom = denominator[count]
        @inbounds sumdenominator[idx] += this_denom
        @inbounds sumweight[idx] += weight[count]
        @inbounds OnlineStats.fit!(moments_per_pdaclass[idx],this_num/this_denom)
    end
  return cnt,sumnumerator,sumdenominator,sumweight,moments_per_pdaclass
end

"""
function unmerge!(s1::T, s2::T) where {T <: OnlineStats.Series}
this function updates the mean and variance of s1
The assumption is that the underlying data of s2 is a subset of s1.
If this is not the case, the result of this function is most likely meaningless
That is to say "subtract" the series s2 from the series s1 and update the mean and variance of s1 accordingly
both s1 and s2 must be of the 'form' Series(OnlineStats.Mean(),OnlineStats.Variance())
Mean must be the fist statistic and Variance the second
Example
vec1=[1.2,3.4,99.1];
vec2=[11.1,3.2];
vecCombined=vcat(vec1,vec2);
#@show isapprox(StatsBase.var(vecCombined,correcte=true),thisVar)

Consider also the merge method: merge!{T <: Series}(s1::T, s2::T, method::Symbol = :append) of OnlineStats
https://github.com/joshday/OnlineStats.jl
"""
 function unmerge!(s1::T, s2::T) where {T <: OnlineStats.Series}
    n2 = OnlineStats.nobs(s2)
    n2 == 0 && return s1
    #OnlineStats.updatecounter!(s1, -n2)
    s1.stats[1].n-=n2
    s1.stats[2].n-=n2

    sizetot=OnlineStats.nobs(s1)+n2
    sizewo=sizetot-n2
    gamma=n2/sizetot
    #we currently access the statistics (i.e. variance and mean) by index which is suboptimal
    #this will throw a typeassert error if the index is not right or the series consists of different types
     o=s1.stats[2]::OnlineStats.Variance #WARNING if this is not of type variance, then all this will fail miserably
     o2=s2.stats[2]::OnlineStats.Variance
    #update variance
        one_minus_gamma=1.0-gamma
        one_minus_gamma_inverted=1/one_minus_gamma
        muwo = o.μ= 1/one_minus_gamma* (o.μ-gamma*o2.μ)
        o.σ2 = one_minus_gamma_inverted*(o.σ2-gamma*(o2.σ2+(o2.μ-muwo)^2*one_minus_gamma))
        o.μ = muwo #note the variance "struct" has a mu parameter as well (in addition to the mean of the series)
    #update mean
    s1.stats[1].μ=muwo
    return s1
end

