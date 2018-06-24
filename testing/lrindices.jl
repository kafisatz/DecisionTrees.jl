using BenchmarkTools 
using StatsBase
#function lrIndices(trnidx::Vector{Int},f::DecisionTrees.PooledArraysDTM.PooledArray{String,T,1},subset::Array) where T<:Unsigned
include("c:\\temp\\2.jl") #PDA module
using .PooledArraysDTM


function lrIndices(idx::Vector{Int},f,subset::Array) 
	#todo tbd, we could restrict the input to U<:Number and T<:Unsigned for the pda f
	l=Vector{Int}(undef,0)
	r=Vector{Int}(undef,0)
	sizehint!(l,length(idx))
	sizehint!(r,length(idx))
    subsetEnd=subset[end]
	for i in idx
		@inbounds thisref=f.refs[i]
		if  thisref<=subsetEnd # this should be equivalent to 'in(thisref,subset)' by construction
			push!(l,i)
		else
			push!(r,i)
		end
	end
	return l,r
end

function lrIndices2(trnidx::Vector{Int},f::PooledArray{String,T,1},subset::Array) where T<:Unsigned
	l=Vector{Int}(undef,0)
	r=Vector{Int}(undef,0)
	sizehint!(l,length(trnidx))
	sizehint!(r,length(trnidx))
    subsetEnd=subset[end]
    a1=find((in)(trnidx)&&(x->f.refs[x]<subsetEnd),1:length(f.refs))
    @show a1
    @show boolIdx
    l=find(boolIdx)
    r=find(.!boolIdx)
	return l,r
end

import Random
nn=30
trni=sample(1:nn,floor(Int,.4*nn),replace=false)
sort!(trni)
@assert trni==unique(trni)

strArr=[randstring(1) for i=1:nn]
pd=PooledArray(strArr)

subsetStr=sample(pd.pool,floor(Int,.35*length(pd.pool)),replace=false)
subset=convert(Vector{eltype(pd.refs)},findall((in)(subsetStr),pd.pool))
sort!(subset)

lll,rrr=lrIndices(trni,pd,subset)
ll2,rr2=lrIndices2(trni,pd,subset)
@assert lrIndices2(trni,pd,subset)==lrIndices(trni,pd,subset)

trniU=convert(Vector{UInt32},trni)
@benchmark lrIndices($trni,$pd,$subset)
@benchmark lrIndices2($trniU,$pd,$subset)


a=rand(nn)
b=rand(nn)
@btime c=$a.*$b
@btime c.=$a.*$b