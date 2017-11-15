this_dir="C:\\Users\\bernhard.konig\\Documents\\ASync\\irobot\\TariffWatch New Repo\\algorithms\\Julia\\Code\\dev\\src";cd(this_dir);push!(LOAD_PATH,this_dir)
using Compat
using BenchmarkTools
using DTM
using DataFrames, StatsBase
n=2_000_000;
r=rand(n);
rint=convert(Array{String},string.(ceil.(r.*200)));
pda=PooledDataArray(rint);

idx=sample(1:n,div(n,4),replace=false,ordered=true)
#idx=[200,2000,5000005];
v=view(pda,idx);
pda_sub=pda[idx]
pdaM=DTM.pdaMod(pda_sub);
szz=size(v,1)
numeratortot=210.*rand(n);
denominatortot=map(x->x<1.2 ? 1.1 : x,((5.*rand(n)).^3));
weighttot=rand(n);
numerator=numeratortot[idx] #210.*rand(szz);
denominator=denominatortot[idx] #210.*rand(szz);
weight=weighttot[idx] #210.*rand(szz);

function aggregate_data_diff_new{T<:AbstractArray{String,1}}(f::T,numerator::Array{Float64,1},denominator::Array{Float64,1},weight::Array{Float64,1})  # = true
    warn("BK i am not sure if this is accurate....")
	lvl_mod=findin(levels(f),f.parent.pool)
	(lo::UInt8, hi::UInt8) = extrema(lvl_mod) #(lo, hi) = extrema(levels(f))
	#@show size(f.parent)
	@assert length(f.indexes)==1 #this must hold
	maho=0.0

	ooo=one(lo)-lo
	vecsize=hi+ooo
    cnt = zeros(Int, vecsize)
    sumnumerator = zeros(Float64, vecsize)
	sumdenominator = zeros(Float64, vecsize)
	sumweight= zeros(Float64, vecsize)

	warn("BK: need to add inbounds here as soon as things work.... ")
	for count = 1:length(v) #in f.indexes[1]
	#for count=1:length(a)
		#@inbounds idx=a[count] + ooo
		idx=f.parent.refs[count] + ooo
		cnt[idx] += 1
		sumnumerator[idx] += numerator[count]
		sumdenominator[idx] += denominator[count]
		sumweight[idx] += weight[count]
    end
  return cnt,sumnumerator,sumdenominator,sumweight

#	for x in f.indexes[1]
		#@show x
		#@show f.parent.refs[x]
#		maho+= float(f.parent.refs[x])
	#end
	#@show f.parent.pool
	#return maho
end

#aggregate_data_diff_new2{T<:AbstractArray{String,1}}(f::T,numerator::Array{Float64,1},denominator::Array{Float64,1},weight::Array{Float64,1})
function aggregate_data_diff_new2(f::SubArray{String,1,DataArrays.PooledDataArray{String,UInt32,1},Tuple{Array{Int64,1}},false},numerator::Array{Float64,1},denominator::Array{Float64,1},weight::Array{Float64,1})  # = true
    #warn("BK i am not sure if this is accurate....")
	#lvl_mod=findin(levels(f),f.parent.pool)
	#@show (lo::UInt8, hi::UInt8) = extrema(lvl_mod) #(lo, hi) = extrema(levels(f))
    #if we never subset the pda and always work on a view of the original data, then lo, hi are always 1: size(pool)!
    lo=one(eltype(f.parent.refs)) #UInt8(1)
    #hi=UInt8(length(f.parent.pool))
    hi=convert(eltype(f.parent.refs),length(f.parent.pool)) #hi=convert(UInt8,length(f.parent.pool))
	#@show size(f.parent)
#	@assert length(f.indexes)==1 #this must hold
	maho=0.0

	ooo=one(lo)-lo
	vecsize=hi+ooo
    cnt = zeros(Int, vecsize)
    sumnumerator = zeros(Float64, vecsize)
	sumdenominator = zeros(Float64, vecsize)
	sumweight= zeros(Float64, vecsize)

#	warn("BK: need to add inbounds here as soon as things work.... ")
	#@inbounds for count::Int64 = 1:length(v) #
    #@inbounds
    @inbounds for count in f.indexes[1]
	#for count=1:length(a)
		#@inbounds idx=a[count] + ooo
		idx=f.parent.refs[count] + ooo
		cnt[idx] += 1
		sumnumerator[idx] += numerator[count]
		sumdenominator[idx] += denominator[count]
		sumweight[idx] += weight[count]
    end
  return cnt,sumnumerator,sumdenominator,sumweight
end

function aggregate_data_diff_new32(f::T,numerator::Array{Float64,1},denominator::Array{Float64,1},weight::Array{Float64,1})  where T<:SubArray
    #warn("BK i am not sure if this is accurate....")
	#lvl_mod=findin(levels(f),f.parent.pool)
	#@show (lo::UInt8, hi::UInt8) = extrema(lvl_mod) #(lo, hi) = extrema(levels(f))
    #if we never subset the pda and always work on a view of the original data, then lo, hi are always 1: size(pool)!
    lo=one(eltype(f.parent.refs)) #UInt8(1)
    #hi=UInt8(length(f.parent.pool))
    hi=convert(eltype(f.parent.refs),length(f.parent.pool)) #hi=convert(UInt8,length(f.parent.pool))
	#@show size(f.parent)
#	@assert length(f.indexes)==1 #this must hold
	maho=0.0

	ooo=one(lo)-lo
	vecsize=hi+ooo
    cnt = zeros(Int, vecsize)
    sumnumerator = zeros(Float64, vecsize)
	sumdenominator = zeros(Float64, vecsize)
	sumweight= zeros(Float64, vecsize)

#	warn("BK: need to add inbounds here as soon as things work.... ")
	#@inbounds for count::Int64 = 1:length(v) #
    #@inbounds
    @inbounds for count in f.indexes[1]
	#for count=1:length(a)
		#@inbounds idx=a[count] + ooo
		idx=f.parent.refs[count] + ooo
		cnt[idx] += 1
		sumnumerator[idx] += numerator[count]
		sumdenominator[idx] += denominator[count]
		sumweight[idx] += weight[count]
    end
  return cnt,sumnumerator,sumdenominator,sumweight
end




function aggregate_data_diff(f::DTM.pdaMod,numerator::Array{Float64,1},denominator::Array{Float64,1},weight::Array{Float64,1})
    #this is the core function of the modelling process
	#besides copying of the data, the vast majority of time is spent in here!
	#most of the time is spent here, if we can improve the for loop below, that would improve performance greatly!
	#one possibility would be to introduce parallelization here (which is not straightforward, I think....)
	a= f.pda.refs
	#this should hold by definition/construction
		#@assert extrema(a)==extrema(levels(f))

	(lo::UInt8, hi::UInt8) = extrema(levels(f)) #(lo, hi) = extrema(levels(f))
	ooo=one(lo)-lo
	vecsize=hi+ooo
    cnt = zeros(Int, vecsize)
    sumnumerator = zeros(Float64, vecsize)
	sumdenominator = zeros(Float64, vecsize)
	sumweight= zeros(Float64, vecsize)
	#note: inbounds increases efficiency here (about a factor of 2), however if the bounds are violated something nasty might happen (quote: If the subscripts are ever out of bounds, you may suffer crashes or silent corruption.)
	for count=1:length(a)
        #i=a[count] #i in a
		#idx=i - lo + 1
		@inbounds idx=a[count] + ooo
		@inbounds cnt[idx] += 1
        @inbounds sumnumerator[idx] += numerator[count]
		@inbounds sumdenominator[idx] += denominator[count]
		@inbounds sumweight[idx] += weight[count]
    end
  return cnt,sumnumerator,sumdenominator,sumweight
end


#c,sn,sd,sw=aggregate_data_diff_new(v,numerator,denominator,weight)
c,sn,sd,sw=aggregate_data_diff_new2(v,numeratortot,denominatortot,weighttot);
c,sn,sd,sw=aggregate_data_diff_new32(v,numeratortot,denominatortot,weighttot);
#@time c,sn,sd,sw=aggregate_data_diff_new(v,numerator,denominator,weight)
@btime c,sn,sd,sw=aggregate_data_diff_new2(v,numeratortot,denominatortot,weighttot);
@btime c,sn,sd,sw=aggregate_data_diff_new32(v,numeratortot,denominatortot,weighttot);


rint2=deepcopy(rint);
rint2[1:end-200].="3.0";
pda2=PooledDataArray(rint2)

idx=collect(40:1232)
unshift!(idx,3)
v2=view(pda2,idx)
pda_trn=deepcopy(v2.parent[v2.indexes[1]])
idx_negation=collect(1:length(pda2))
deleteat!(idx_negation,idx)
pda_val=deepcopy(v2.parent[idx_negation])
@btime v2.parent[idx_negation]

#some checks
@assert size(rint2,1)==size(pda_val,1)+size(pda_trn,1)
@assert isequal(pda_val.pool,pda_trn.pool)
@assert isequal(sort(pda_trn.pool),sort(unique(rint2)))
testv=sort(parse.(Float64,pda_trn.pool))
@assert eltype(pda_val.refs)==eltype(pda_trn.refs)
eltype(pda_trn.refs)


@code_warntype aggregate_data_diff_new2(v,numeratortot,denominatortot,weighttot)

@code_warntype aggregate_data_diff(v,numerator,denominator,weight)
c2,sn2,sd2,sw2=aggregate_data_diff(pdaM,numerator,denominator,weight)
@time c2,sn2,sd2,sw2=aggregate_data_diff(pdaM,numerator,denominator,weight)


@btime aggregate_data_diff(pdaM,numerator,denominator,weight);
@btime aggregate_data_diff_new2(v,numeratortot,denominatortot,weighttot);


@benchmark aggregate_data_diff(pdaM,numerator,denominator,weight)
@benchmark aggregate_data_diff_new2(v,numeratortot,denominatortot,weighttot)
