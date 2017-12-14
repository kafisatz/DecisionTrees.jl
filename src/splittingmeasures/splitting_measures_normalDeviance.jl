#=
	info("""
	Splitting Function = NormalDeviance
	for the normal deviance
	the concept should be like this
	calc var for each indidvidual 'subset' -> using COUNTSORT similar to the difference preparation
	add up subsets however you like (for num vars this is already "it")
	substract subsets as needed (gray code loop) with the function below
	"""
	)
	warn("TBD: check allocation and time of aggregate_data_normal_deviance")

	#this runs a small version to compile all functions
	loc="C:\\jdata\\datafolder\\mortgage\\"
	ARGS=[string(loc,"jmortgage.settings.csv") string(loc,"jmortgage.CSV") string("out_jmortgage")];
	ARGS=[string(loc,"jmortgage.settings_profiling.csv") string(loc,"jmortgage.CSV") string("out_jmortgage")];
	resbool=DecisionTrees.run_model(ARGS) #small version to compile functions

	using Iterators
	using DTM
	#import DecisionTrees.pdaMod
	using OnlineStats,SQLite,DataFrames, ProgressMeter, Iterators,PyCall,JLD2
	dix=load(string(loc,"jmortgage.jld2"));


	numerator=dix["numeratortrn"];
	denominator=dix["denominatortrn"];
	weight=dix["weighttrn"];
	trn_charfeatures_PDA=dix["trn_charfeatures_PDA"];
	f=trn_charfeatures_PDA[1]
	 cnt,sumnumerator,sumdenominator,sumweight,moments_per_pdaclass=aggregate_data_normal_deviance(f,numerator,denominator,weight)
	features=f
	minweight=100.0
	crit=DifferenceSplit()
	feature_column_id=-1
	randomweight=0.05
	catSortByThreshold=7
	catSortBy=SortByMean()

=#

function calculateSplitValue(a::NormalDevianceSplit,number_of_char_features::Int,labellist::Array{UInt8,1},sumnumerator::Array{Float64,1},sumdenominator::Array{Float64,1},sumweight::Array{Float64,1},countlistfloat::Array{Float64,1},minweight::Float64,subs::DTSubsets,moments_per_pdaclass)
  warn("229 this will not work yet.....")
  warn("randomw should be zero here")
  save("C:\\temp\\myf.jld2","a",a,"labellist",labellist,"sumnumerator",sumnumerator,"sumdenominator",sumdenominator,"sumweight",sumweight,"countlistfloat",countlistfloat,"minweight",minweight,"subs",subs,"moments_per_pdaclass",moments_per_pdaclass)

#here randomweight==0
#for subsets, exhaustive search with flipping members (gray code) or "increasing" subset search ({1}, {1,2}, {1,2,3}, .... {1,2,3, ....., n-1,2})
#all input lists (labellist,sumnumerator,sumdenominator,sumweight,countlistfloat) need to be sorted in the same manner
#convention, in the beginning everything is on the right side
elementsInLeftChildBV=BitVector(length(labellist));fill!(elementsInLeftChildBV,false) #indicates which classes belong to right child
chosen_subset_bitarray=BitVector(0)
#momentstot=deepcopy(moments_per_pdaclass[1])
#for idx=2:length(moments_per_pdaclass)
#  merge!(momentstot,moments_per_pdaclass[idx])
#end
#vartot=value(momentstot.stats[2])
#numtot=sum(sumnumerator)
#denomtot=sum(sumdenominator)
weighttot=sum(sumweight)
weighttot_minw=weighttot-minweight

momentsl=Series(OnlineStats.Mean(),OnlineStats.Variance())
momentsr=Series(OnlineStats.Mean(),OnlineStats.Variance())
sumwl=0.0
varl=varr=val=-Inf
chosen_sumwl=NaN
  for i in subs
  #@show "vv=$(i)" #should be the element which switches
  #i is an index, it indicates which element of labellist will flip sides for the next calculation
  switching_class=moments_per_pdaclass[i]
  if elementsInLeftChildBV[i]
    #update elementsInLeftChildBV, i.e. toggle the value of the ith component, $=XOR
    elementsInLeftChildBV[i]=xor(elementsInLeftChildBV[i],true) #updating needs to occur here before we copy the array (in case it is a better split than what we have seen so far)
    #move class from left to right side
    unmerge!(momentsl,switching_class)
    merge!(momentsr,switching_class)
    #sumnl-=sumnumerator[i]
	  #sumdl-=sumdenominator[i]
	  sumwl-=sumweight[i]
  else
      #update elementsInLeftChildBV, i.e. toggle the value of the ith component, $=XOR
      elementsInLeftChildBV[i]=xor(elementsInLeftChildBV[i],true)
      #move class from right to left side
      merge!(momentsl,switching_class)
      unmerge!(momentsr,switching_class)
      #sumnl+=sumnumerator[i]
    	#sumdl+=sumdenominator[i]
    	sumwl+=sumweight[i]
  end
    if (sumwl>minweight)&&(weighttot_minw>sumwl) #do we have enough exposure? is the split valid?
	  #vold=valnew
      varl=value(momentsl.stats[2])
      varr=value(momentsl.stats[2])
	     valnew=varl+varr #abs(sumnl/sumdl-(numtot-sumnl)/(denomtot-sumdl))
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
  chosen_subset=Array{UInt}(0)
end
return val,chosen_subset,chosen_sumwl,weighttot-chosen_sumwl
end



function aggregate_data_normal_deviance(f::pdaMod,numerator::Array{Float64,1},denominator::Array{Float64,1},weight::Array{Float64,1})
    #this is the core function of the modelling process
	#besides copying of the data, the vast majority of time is spent in here!
	#most of the time is spent here, if we can improve the for loop below, that would improve performance greatly!
	#one possibility would be to introduce parallelization here (which is not straightforward, I think....)
	a= f.pda.refs
	(lo::UInt8, hi::UInt8) = extrema(levels(f))  
  ooo=one(lo)-lo
	vecsize=hi+ooo
  cnt = zeros(Int, vecsize)
  sumnumerator = zeros(Float64, vecsize)
	sumdenominator = zeros(Float64, vecsize)
	sumweight= zeros(Float64, vecsize)
	warn("check if we need sumdenom and sumnom at all")
	moments_per_pdaclass=Array{OnlineStats.Series}(vecsize)
	for i=1:vecsize
		moments_per_pdaclass[i]=Series(OnlineStats.Mean(),OnlineStats.Variance())
	end
	#note: inbounds increases efficiency here (about a factor of 2), however if the bounds are violated something nasty might happen (quote: If the subscripts are ever out of bounds, you may suffer crashes or silent corruption.)
	for count=1:length(a)
        #i=a[count] #i in a
		#idx=i - lo + 1
		@inbounds idx=a[count] + ooo
		@inbounds cnt[idx] += 1
        @inbounds this_num = numerator[count]
		@inbounds sumnumerator[idx] += this_num
		@inbounds this_denom = denominator[count]
		@inbounds sumdenominator[idx] += this_denom
		@inbounds sumweight[idx] += weight[count]
		@inbounds fit!(moments_per_pdaclass[idx],this_num/this_denom)
    end
  return cnt,sumnumerator,sumdenominator,sumweight,moments_per_pdaclass
end

# merge method merge!{T <: Series}(s1::T, s2::T, method::Symbol = :append)
# #this extends the functionality of the package
# https://github.com/joshday/OnlineStats.jl
 function unmerge!{T <: OnlineStats.Series}(s1::T, s2::T)
  #this function updates the variance of s1
  #we assume that the underlying data of s2 is a subset of s1
  #thus we "substract" the series s2 from the series s1 and update the mean and variance of s1 accordingly

  #IMPORTANT NOTE, this will only work for mean and variance (OnlineStats provides a wealth of additional metrics)
    n2 = OnlineStats.nobs(s2)
    n2 == 0 && return s1
	OnlineStats.updatecounter!(s1, -n2)

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
		#muwo = o.μ= 1/(1.0-gamma)* (o.μ-gamma*o2.μ)
		muwo = o.μ= 1/one_minus_gamma* (o.μ-gamma*o2.μ)
		#o.σ2 = 1/(1.0-gamma)*(o.σ2-gamma*(o2.σ2+(o2.μ-muwo)^2*(1.0-gamma)))
		o.σ2 = one_minus_gamma_inverted*(o.σ2-gamma*(o2.σ2+(o2.μ-muwo)^2*one_minus_gamma))
		o.μ = muwo #note the variance "struct" has a mu parameter as well (in addition to the mean of the series)
	#update mean
	s1.stats[1].μ=muwo
	#update nobs
	o.nobs -= o2.nobs
    return s1
end


function calculateSplitValue(a::NormalDevianceSplit,number_of_char_features::Int,labellist::Array{UInt8,1},sumnumerator::Array{Float64,1},sumdenominator::Array{Float64,1},sumweight::Array{Float64,1},countlistfloat::Array{Float64,1},minweight::Float64,subs::DTSubsets,feature_column_id::Int64,moments_per_pdaclass)
warn("221 this will not work yet.....")
warn("randomw should be postive here")

save("C:\\temp\\myf.jld2","a",a,"labellist",labellist,"sumnumerator",sumnumerator,"sumdenominator",sumdenominator,"sumweight",sumweight,"countlistfloat",countlistfloat,"minweight",minweight,"subs",subs,"feature_column_id",feature_column_id,"moments_per_pdaclass",moments_per_pdaclass)
#here randomweight>0
#for subsets, exhaustive search with flipping members (gray code) or "increasing" subset search ({1}, {1,2}, {1,2,3}, .... {1,2,3, ....., n-1,2})
#all input lists (labellist,sumnumerator,sumdenominator,sumweight,countlistfloat) need to be sorted in the same manner
#convention, in the beginning everything is on the right side
elementsInLeftChildBV=BitVector(length(labellist));fill!(elementsInLeftChildBV,false) #indicates which classes belong to right child

#momentstot=deepcopy(moments_per_pdaclass[1])
#for idx=2:length(moments_per_pdaclass)
#  merge!(momentstot,moments_per_pdaclass[idx])
#end
#vartot=value(momentstot.stats[2])
weighttot=sum(sumweight)
weighttot_minw=weighttot-minweight
momentsl=Series(OnlineStats.Mean(),OnlineStats.Variance())
momentsr=Series(OnlineStats.Mean(),OnlineStats.Variance())
sumwl=0.0
this_splitlist=Array{Splitdef}(0)

for i in subs
    #i is an index, it indicates which element of labellist will flip sides for the next calculation
 if elementsInLeftChildBV[i]
    #update elementsInLeftChildBV, i.e. toggle the value of the ith component, $=XOR
    elementsInLeftChildBV[i]=xor(elementsInLeftChildBV[i],true) #updating needs to occur here before we copy the array (in case it is a better split than what we have seen so far)
    #move class from left to right side
    unmerge!(momentsl,switching_class)
    merge!(momentsr,switching_class)
	  sumwl-=sumweight[i]
  else
    #update elementsInLeftChildBV, i.e. toggle the value of the ith component, $=XOR
    elementsInLeftChildBV[i]=xor(elementsInLeftChildBV[i],true)
    #move class from right to left side
    merge!(momentsl,switching_class)
    unmerge!(momentsr,switching_class)
	  sumwl+=sumweight[i]
  end
    if (sumwl>minweight)&&(weighttot_minw>sumwl) #do we have enough exposure? is the split valid?
      varl=value(momentsl.stats[2])
      varr=value(momentsl.stats[2])
	    valnew=varl+varr
	  push!(this_splitlist,Splitdef(feature_column_id,labellist[elementsInLeftChildBV],valnew,sumwl,weighttot-sumwl))
    end
end
return this_splitlist
end
