using BenchmarkTools 
using StatsBase
#function lrIndices(trnidx::Vector{Int},f::DecisionTrees.PooledArraysDTM.PooledArray{String,T,1},subset::Array) where T<:Unsigned
include("c:\\temp\\2.jl") #PDA module
using .PooledArraysDTM


"""
returns l and r indices (l corresponds to all elements of trnidx that 'match' subset)
This verison of the function is for String / Categorical variables
"""
function lrIndices(trnidx::Vector{Int},f::PooledArray{String,T,1},subset::Array) where T<:Unsigned
	l=Vector{Int}(undef,0)
	r=Vector{Int}(undef,0)
	sizehint!(l,length(trnidx))
	sizehint!(r,length(trnidx))
	for i in trnidx
		@inbounds thisref=f.refs[i]
		if in(thisref,subset)
			push!(l,i)
		else
			push!(r,i)
		end
	end
	return l,r
end

"""
returns l and r indices (l corresponds to all elements of trnidx that 'match' subset)
This verison of the function is for String / Categorical variables
"""
function lrIndices2(trnidx::Vector{Int},f::PooledArray{String,T,1},subset::Array) where T<:Unsigned
	l=Vector{Int}(undef,0)
	r=Vector{Int}(undef,0)
	sizehint!(l,length(trnidx))
	sizehint!(r,length(trnidx))
    subsetEnd=subset[end]
	for i in trnidx
		@inbounds thisref=f.refs[i]
		if thisref<=subsetEnd
			push!(l,i)
		else
			push!(r,i)
		end
	end
	return l,r
end

import Random


for klU=1:4

nn=300000
trni=sample(1:nn,floor(Int,.4*nn),replace=false)
sort!(trni)
@assert trni==unique(trni)

strArr=[randstring(1) for i=1:nn]
pd=PooledArray(strArr)

subsetStr=sample(pd.pool,floor(Int,.35*length(pd.pool)),replace=false)
subset=convert(Vector{eltype(pd.refs)},findall((in)(subsetStr),pd.pool))
sort!(subset)

subsetContiguous=convert(Vector{UInt8},1:ceil(Int,rand()*length(pd.pool)))
isC=true #rand()<.5
if isC
    subset=subsetContiguous
end

lll,rrr=lrIndices(trni,pd,subset)
ll2,rr2=lrIndices2(trni,pd,subset)
@assert lrIndices2(trni,pd,subset)==lrIndices(trni,pd,subset)

@show isC
@btime lrIndices($trni,$pd,$subset)
@btime lrIndices2($trni,$pd,$subset)
end

"""
isContiguous(v) returns true if the vector v is of the form collect(n:m)
"""
function isContiguous(subset::Vector)
    mi=subset[1]
    ma=subset[end]
    return (length(subset)==(ma-mi+1))&&issorted(subset)
end

@assert isContiguous(collect(1:23))
@assert isContiguous(collect(22:23))
@assert isContiguous(collect(23:23))
@assert !isContiguous([1,9,3])
@assert !isContiguous([1,3,2])


