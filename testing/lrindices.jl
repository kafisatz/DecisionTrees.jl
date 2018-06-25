using BenchmarkTools 
using StatsBase
#function lrIndices(trnidx::Vector{Int},f::DecisionTrees.PooledArraysDTM.PooledArray{String,T,1},subset::Array) where T<:Unsigned
include("c:\\temp\\2.jl") #PDA module
using .PooledArraysDTM


"""
function lrIndices!(trnidx::Vector{Int},f,subset::Array)
modifies trnidx to 'become' the left child index, returns r
l and r are defined by subset. l corresponds to the observations 'in' subset
"""
function lrIndices!(trnidx::Vector{Int},f,subset::Array)
if eltype(f.pool)<:Number
    #NOTE: subset may not necessarily start at one here (because after a few splits, certain 'middle' values of the pool may not exist in the current node)
    #still, the definition of the split should be of the form x<subset[end] -> go left (otherwise go right)
    #therefore we can consider this function for numerical variables
    return lrIndicesForContiguousSubsetStartingAtONE!(trnidx,f,subset)
else
    if isContiguous(subset)
        if (length(subset)>0)&&(isone(subset[1]))
            return lrIndicesForContiguousSubsetStartingAtONE!(trnidx,f,subset)
        else 
            return lrIndicesForContiguousSubsetNOTStartingAtONE!(trnidx,f,subset)
        end
    else
        return lrIndicesDefault!(trnidx,f,subset)
    end
end
end

"""
function lrIndices!(trnidx::Vector{Int},f,subset::Array)
modifies trnidx to 'become' the left child index, returns r
l and r are defined by subset. l corresponds to the observations 'in' subset
"""
function lrIndicesDefault!(trnidx::Vector{Int},f,subset::Array)
	r=Vector{Int}(undef,0)
	sizehint!(r,length(trnidx))
	indexSize=length(trnidx)
    for jj=indexSize:-1:1
        @inbounds i=trnidx[jj]
		@inbounds thisref=f.refs[i]
		if !in(thisref,subset)
			pushfirst!(r,i)
			deleteat!(trnidx,jj)
		end
	end
	return r
end

"""
function lrIndicesForContiguousSubsetStartingAtONE!(trnidx::Vector{Int},f,subset::Array)
this is a similar function as lrIndices! (->see the docs of lrIndices!)
Here subset needs to be of the form collect(1:n)
"""
function lrIndicesForContiguousSubsetStartingAtONE!(trnidx::Vector{Int},f,subset::Array)
	r=Vector{Int}(undef,0)
	sizehint!(r,length(trnidx))
    subsetEnd=subset[end]
    indexSize=length(trnidx)
    for jj=indexSize:-1:1
        @inbounds i=trnidx[jj]
		@inbounds thisref=f.refs[i]
		if !(thisref<=subsetEnd)
			pushfirst!(r,i)
			deleteat!(trnidx,jj)
		end
	end
	return r
end


"""
this is a similar function as lrIndices! (->see the docs of lrIndices!)
Here subset needs to be of the form collect(m:n)
"""
function lrIndicesForContiguousSubsetNOTStartingAtONE!(trnidx::Vector{Int},f,subset::Array)
	r=Vector{Int}(undef,0)
	sizehint!(r,length(trnidx))
    subsetEnd=subset[end]
    indexSize=length(trnidx)    
	subsetEnd=subset[end]
	subsetStart=subset[1]
    for jj=indexSize:-1:1
        @inbounds i=trnidx[jj]
		@inbounds thisref=f.refs[i]
		#if !((thisref<=subsetEnd)&&(thisref>=subsetStart))
        if (thisref>subsetEnd)||(thisref<subsetStart)
			pushfirst!(r,i)
			deleteat!(trnidx,jj)
		end
	end
	return r
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
nn=15
trni=sample(1:nn,floor(Int,.4*nn),replace=false)
sort!(trni)
strArr=[randstring(1) for i=1:nn]
floatArr=float.(abs.(rand(Int,nn).%9))
pdStr=PooledArray(strArr)
pdFloat=PooledArray(floatArr)

for klU=1:400000
    #@show klU
    pd = rand()<.5 ? pdStr : pdFloat

    subsetStrOrFloat=sample(pd.pool,4,replace=false)
    subset=convert(Vector{eltype(pd.refs)},findall((in)(subsetStrOrFloat),pd.pool))
    sort!(subset)

    #if isContiguous(subset)
           
        whereStartsAtONE+=isone(subset[1])
        ntest+=1

        trniWillBeModified=deepcopy(trni)
        rrr=lrIndices!(trniWillBeModified,pd,subset)
        if eltype(pd.pool)==Float64
            ll2,rr2=DecisionTrees.lrIndicesForContiguousSubsetStartingAtONE(trni,pd,subset)
        else
            ll2,rr2=DecisionTrees.lrIndices(trni,pd,subset)
        end
        #=
        @show ll2,rr2
        @show subset
        @show pd
        @show pd.pool
        @show sort!(unique(pd.refs))
        @show ll2
        @show rrr
        @show trni
        @show trniWillBeModified
        =#
        @assert ll2==trniWillBeModified 
        @assert rr2==rrr

        #@btime lrIndices($trni,$pd,$subset)
        #@btime lrIndices2($trni,$pd,$subset)
    #end
end
