using BenchmarkTools 
using StatsBase
#function lrIndices(trnidx::Vector{Int},f::DecisionTrees.PooledArraysDTM.PooledArray{String,T,1},subset::Array) where T<:Unsigned
include("c:\\temp\\2.jl") #PDA module
using .PooledArraysDTM


"""
returns l and r indices (l corresponds to all elements of trnidx that 'match' subset)
Here subset needs to be of the form collect(m:n)
"""
function lrIndicesForContiguousSubset(trnidx::Vector{Int},f,subset::Array)
	if (length(subset)>0)&&(isone(subset[1]))
    	return lrIndicesForContiguousSubsetStartingAtONE(trnidx,f,subset)
	else 
		return lrIndicesForContiguousSubsetNOTStartingAtONE(trnidx,f,subset)
	end
end

"""
returns l and r indices (l corresponds to all elements of trnidx that 'match' subset)
Here subset needs to be of the form collect(1:n)
"""
function lrIndicesForContiguousSubsetStartingAtONE(trnidx::Vector{Int},f,subset::Array)
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


"""
returns l and r indices (l corresponds to all elements of trnidx that 'match' subset)
Here subset needs to be of the form collect(m:n)
"""
function lrIndicesForContiguousSubsetNOTStartingAtONE(trnidx::Vector{Int},f,subset::Array)
	l=Vector{Int}(undef,0)
	r=Vector{Int}(undef,0)
	sizehint!(l,length(trnidx))
	sizehint!(r,length(trnidx))
	subsetEnd=subset[end]
	subsetStart=subset[1]
	for i in trnidx
		@inbounds thisref=f.refs[i]
		if (thisref<=subsetEnd)&&(thisref>=subsetStart)
			push!(l,i)
		else
			push!(r,i)
		end
	end
	return l,r
end


"""
isContiguous(v) returns true if the vector v is of the form collect(n:m)
"""
function isContiguous(subset::Vector)
    mi=subset[1]
    ma=subset[end]
    return (length(subset)==(ma-mi+1))&&issorted(subset)
end

import Random

ntest=0
whereStartsAtONE=0
nn=300
trni=sample(1:nn,floor(Int,.4*nn),replace=false)
sort!(trni)
strArr=[randstring(1) for i=1:nn]
pd=PooledArray(strArr)

for klU=1:4000000

subsetStr=sample(pd.pool,4,replace=false)
subset=convert(Vector{eltype(pd.refs)},findall((in)(subsetStr),pd.pool))
sort!(subset)

#subsetContiguous=convert(Vector{UInt8},5:ceil(Int,rand()*length(pd.pool)))
#subset=subsetContiguous 
    
#isC=true #rand()<.5
if isContiguous(subset)
   
whereStartsAtONE+=isone(subset[1])
ntest+=1


lll,rrr=lrIndices(trni,pd,subset)
ll2,rr2=lrIndicesForContiguousSubset(trni,pd,subset)
@assert lrIndicesForContiguousSubset(trni,pd,subset)==lrIndices(trni,pd,subset)

#@btime lrIndices($trni,$pd,$subset)
#@btime lrIndices2($trni,$pd,$subset)
end
end
@show ntest,whereStartsAtONE

@assert isContiguous(collect(1:23))
@assert isContiguous(collect(22:23))
@assert isContiguous(collect(23:23))
@assert !isContiguous([1,9,3])
@assert !isContiguous([1,3,2])


