__precompile__()

module PooledArraysDTM
# BK this is commit (the newer versions use Dicts instead of arrays for the pool which is not compatible with DTM)
# https://github.com/JuliaComputing/PooledArrays.jl/tree/324263b39be8dc1982b0971f6031cbde92222268

#Here is the copy of the MIT license of PooledArrays

#=
"""
The MIT License (MIT)

Copyright (c) 2016 Julia Computing, Inc.
Some portions based on DataArrays.jl, Copyright (c) 2012: Harlan Harris, EPRI (Tom Short's code), Chris DuBois, John Myles White, and other contributors.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
""" 
=#



export PooledArray

##############################################################################
##
## PooledArray type definition
##
##############################################################################

const DEFAULT_POOLED_REF_TYPE = UInt32

# This is used as a wrapper during PooledArray construction only, to distinguish
# arrays of pool indices from normal arrays
mutable struct RefArray{R}
    a::R
end

struct PooledArray{T, R<:Integer, N, RA} <: AbstractArray{T, N}
    refs::RA
    pool::Vector{T}

    function PooledArray(rs::RefArray{RA}, p::Vector{T}) where {T,R,N,RA<:AbstractArray{R, N}}
        # refs mustn't overflow pool
        if length(rs.a) > 0 && maximum(rs.a) > prod(size(p))
            throw(ArgumentError("Reference array points beyond the end of the pool"))
        end
        new{T,R,N,RA}(rs.a,p)
    end
end
const PooledVector{T,R} = PooledArray{T,R,1}
const PooledMatrix{T,R} = PooledArray{T,R,2}

##############################################################################
##
## PooledArray constructors
##
# Algorithm:
# * Start with:
#   * A null pool
#   * A pre-allocated refs
#   * A hash from T to Int
# * Iterate over d
#   * If value of d in pool already, set the refs accordingly
#   * If value is new, add it to the pool, then set refs
##############################################################################

# Echo inner constructor as an outer constructor
function PooledArray(refs::RefArray{R}, pool::Vector{T}) where {T,R}
    PooledArray{T,eltype(R),ndims(R),R}(refs, pool)
end

# A no-op constructor
PooledArray(d::PooledArray) = d

# Constructor from array, pool, and ref type
function PooledArray(d::AbstractArray, pool::Vector{T},
                     r::Type{R} = DEFAULT_POOLED_REF_TYPE) where {T,R<:Integer}
    if length(pool) > typemax(R)
        throw(ArgumentError("Cannot construct a PooledVector with type $R with a pool of size $(length(pool))"))
    end

    newrefs = Array{R}(size(d))
    poolref = Dict{T, R}()

    # loop through once to fill the poolref dict
    for i = 1:length(pool)
        poolref[pool[i]] = i
    end

    # fill in newrefs
    for i = 1:length(d)
        newrefs[i] = get(poolref, d[i], 0)
    end
    return PooledArray(RefArray(newrefs), pool)
end

"""
    PooledArray(array, [reftype])

Convert the given array to a PooledArray where each element will be referenced
as an integer of the given type. If no `reftype` is specified one is chosen
automatically based on the number of unique elements.
"""
PooledArray

function PooledArray{T}(d::AbstractArray, r::Type{R}) where {T,R<:Integer}
    pool = convert(Vector{T}, unique(d))
    if hasmethod(isless, (T, T))
        sort!(pool)
    end
    PooledArray(d, pool, r)
end

function PooledArray{T}(d::AbstractArray) where T
    pool = convert(Vector{T}, unique(d))
    u = length(pool)
    R = u < typemax(UInt8) ? UInt8 :
        u < typemax(UInt16) ? UInt16 : UInt32
    if hasmethod(isless, (T, T))
        sort!(pool)
    end
    PooledArray(d, pool, R)
end

PooledArray(d::AbstractArray{T}, r::Type{R}) where {T,R<:Integer} = PooledArray{T}(d, r)
PooledArray(d::AbstractArray{T}) where {T} = PooledArray{T}(d)

# Construct an empty PooledVector of a specific type
PooledArray(t::Type) = PooledArray(Array(t,0))
PooledArray(t::Type, r::Type{R}) where {R<:Integer} = PooledArray(Array(t,0), r)

##############################################################################
##
## Basic interface functions
##
##############################################################################

Base.size(pa::PooledArray) = size(pa.refs)
Base.length(pa::PooledArray) = length(pa.refs)
Base.endof(pa::PooledArray) = endof(pa.refs)

Base.copy(pa::PooledArray) = PooledArray(RefArray(copy(pa.refs)), copy(pa.pool))
# TODO: Implement copy_to()

function Base.resize!(pa::PooledArray{T,R,1}, n::Integer) where {T,R}
    oldn = length(pa.refs)
    resize!(pa.refs, n)
    pa.refs[oldn+1:n] = zero(R)
    pa
end

Base.reverse(x::PooledArray) = PooledArray(RefArray(reverse(x.refs)), x.pool)

function Base.permute!!(x::PooledArray, p::AbstractVector{T}) where T<:Integer
    Base.permute!!(x.refs, p)
    x
end

function Base.invpermute!(x::PooledArray, p::AbstractVector{T}) where T<:Integer
    Base.invpermute!(x.refs, p)
    x
end

function Base.similar(pa::PooledArray{T,R}, S::Type, dims::Dims) where {T,R}
    PooledArray(RefArray(zeros(R, dims)), S[])
end

Base.find(pdv::PooledVector{Bool}) = findall(convert(Vector{Bool}, pdv, false))

##############################################################################
##
## map
## Calls `f` only once per pool entry, plus preserves sortedness.
##
##############################################################################

function Base.map(f, x::PooledArray{T,R}) where {T,R<:Integer}
    newpool = map(f, x.pool)
    uniquedpool = similar(newpool, 0)
    lookup = similar(x.pool, Int)
    n = 1
    d = Dict{eltype(newpool), Int}()
    @inbounds for i = 1:length(newpool)
        el = newpool[i]
        idx = get(d, el, 0)
        if idx == 0
            d[el] = n
            idx = n
            n += 1
            push!(uniquedpool, el)
        end
        lookup[i] = idx
    end
    sp = sortperm(uniquedpool)
    isp::Vector{Int} = invperm(sp)
    PooledArray(RefArray(R[ isp[lookup[i]] for i in x.refs ]),
                uniquedpool[sp])
end

##############################################################################
##
## Sorting can use the pool to speed things up
##
##############################################################################

function groupsort_indexer(x::AbstractVector, ngroups::Integer)
    # translated from Wes McKinney's groupsort_indexer in pandas (file: src/groupby.pyx).

    # count group sizes, location 0 for NA
    n = length(x)
    # counts = x.pool
    counts = fill(0, ngroups + 1)
    @inbounds for i = 1:n
        counts[x[i] + 1] += 1
    end

    # mark the start of each contiguous group of like-indexed data
    where = fill(1, ngroups + 1)
    @inbounds for i = 2:ngroups+1
        where[i] = where[i - 1] + counts[i - 1]
    end

    # this is our indexer
    result = fill(0, n)
    @inbounds for i = 1:n
        label = x[i] + 1
        result[where[label]] = i
        where[label] += 1
    end
    result, where, counts
end

groupsort_indexer(pv::PooledVector) = groupsort_indexer(pv.refs, length(pv.pool))

function Base.sortperm(pa::PooledArray; alg::Base.Sort.Algorithm=Base.Sort.DEFAULT_UNSTABLE,
                       lt::Function=isless, by::Function=identity,
                       rev::Bool=false, order=Base.Sort.Forward)
    order = Base.ord(lt, by, rev, order)

    # TODO handle custom ordering efficiently
    if !isa(order, Base.Order.ForwardOrdering) && !isa(order, Base.Order.ReverseOrdering)
        return sort!(collect(1:length(pa)), alg, Base.Order.Perm(order,pa))
    end

    perm = groupsort_indexer(pa)[1]
    isa(order, Base.Order.ReverseOrdering) && reverse!(perm)
    perm
end

Base.sort(pa::PooledArray; kw...) = pa[sortperm(pa; kw...)]

#type FastPerm{O<:Base.Sort.Ordering,V<:AbstractVector} <: Base.Sort.Ordering
#    ord::O
#    vec::V
#end
#Base.sortperm{V}(x::AbstractVector, a::Base.Sort.Algorithm, o::FastPerm{Base.Sort.ForwardOrdering,V}) = x[sortperm(o.vec)]
#Base.sortperm{V}(x::AbstractVector, a::Base.Sort.Algorithm, o::FastPerm{Base.Sort.ReverseOrdering,V}) = x[reverse(sortperm(o.vec))]
#Perm{O<:Base.Sort.Ordering}(o::O, v::PooledVector) = FastPerm(o, v)

##############################################################################
##
## conversions
##
##############################################################################

Base.convert(::Type{PooledArray{S,R1,N}}, pa::PooledArray{T,R2,N}) where {S,T,R1<:Integer,R2<:Integer,N} =
    PooledArray(RefArray(convert(Array{R1,N}, pa.refs)), convert(Vector{S}, pa.pool))
Base.convert(::Type{PooledArray{S,R,N}}, pa::PooledArray{T,R,N}) where {S,T,R<:Integer,N} =
    PooledArray(RefArray(copy(pa.refs)), convert(Vector{S}, pa.pool))
Base.convert(::Type{PooledArray{T,R,N}}, pa::PooledArray{T,R,N}) where {T,R<:Integer,N} = pa
Base.convert(::Type{PooledArray{S,R1}}, pa::PooledArray{T,R2,N}) where {S,T,R1<:Integer,R2<:Integer,N} =
    convert(PooledArray{S,R1,N}, pa)
Base.convert(::Type{PooledArray{S}}, pa::PooledArray{T,R,N}) where {S,T,R<:Integer,N} =
    convert(PooledArray{S,R,N}, pa)
Base.convert(::Type{PooledArray}, pa::PooledArray{T,R,N}) where {T,R<:Integer,N} = pa

Base.convert(::Type{PooledArray{S,R,N}}, a::AbstractArray{T,N}) where {S,T,R<:Integer,N} =
    PooledArray(convert(Array{S,N}, a), R)
Base.convert(::Type{PooledArray{S,R}}, a::AbstractArray{T,N}) where {S,T,R<:Integer,N} =
    PooledArray(convert(Array{S,N}, a), R)
Base.convert(::Type{PooledArray{S}}, a::AbstractArray{T,N}) where {S,T,N} =
    PooledArray(convert(Array{S,N}, a))
Base.convert(::Type{PooledArray}, a::AbstractArray) =
    PooledArray(a)

function Base.convert(::Type{Array{S, N}}, pa::PooledArray{T, R, N}) where {S, T, R, N}
    res = Array{S}(size(pa))
    for i in 1:length(pa)
        if pa.refs[i] != 0
            res[i] = pa.pool[pa.refs[i]]
        end
    end
    return res
end

Base.convert(::Type{Vector}, pv::PooledVector{T, R}) where {T, R} = convert(Array{T, 1}, pv)

Base.convert(::Type{Matrix}, pm::PooledMatrix{T, R}) where {T, R} = convert(Array{T, 2}, pm)

Base.convert(::Type{Array}, pa::PooledArray{T, R, N}) where {T, R, N} = convert(Array{T, N}, pa)

##############################################################################
##
## indexing
##
##############################################################################

# Scalar case
function Base.getindex(pa::PooledArray, I::Integer)
    return pa.pool[getindex(pa.refs, I)]
end
function Base.getindex(pa::PooledArray, I::Integer...)
    return pa.pool[getindex(pa.refs, I...)]
end

# Vector case
function Base.getindex(A::PooledArray, I::Union{Real,AbstractVector}...)
    PooledArray(RefArray(getindex(A.refs, I...)), copy(A.pool))
end

# Dispatch our implementation for these cases instead of Base
Base.getindex(A::PooledArray, I::AbstractVector) =
    PooledArray(RefArray(getindex(A.refs, I)), copy(A.pool))
Base.getindex(A::PooledArray, I::AbstractArray) =
    PooledArray(RefArray(getindex(A.refs, I)), copy(A.pool))

##############################################################################
##
## setindex!() definitions
##
##############################################################################

function getpoolidx(pa::PooledArray{T,R}, val::Any) where {T,R}
    val::T = convert(T,val)
    pool_idx = searchsortedfirst(pa.pool, val)
    if pool_idx > length(pa.pool) || pa.pool[pool_idx] != val
        pool_idx = unsafe_pool_push!(pa, val)
    end
    return pool_idx
end

function unsafe_pool_push!(pa::PooledArray{T,R}, val) where {T,R}
    push!(pa.pool, val)
    pool_idx = length(pa.pool)
    if pool_idx > typemax(R)
        throw(ErrorException(string(
            "You're using a PooledArray with ref type $R, which can only hold $(Int(typemax(R))) values,\n",
            "and you just tried to add the $(typemax(R)+1)th reference.  Please change the ref type\n",
            "to a larger int type, or use the default ref type ($DEFAULT_POOLED_REF_TYPE)."
           )))
    end
    if pool_idx > 1 && isless(val, pa.pool[pool_idx-1])
        # maintain sorted order
        sp = sortperm(pa.pool)
        isp = invperm(sp)
        refs = pa.refs
        for i = 1:length(refs)
            # after resize we might have some 0s
            if refs[i] != 0
                @inbounds refs[i] = isp[refs[i]]
            end
        end
        pool_idx = isp[pool_idx]
        copy!(pa.pool, pa.pool[sp])
    end
    pool_idx
end

function Base.setindex!(x::PooledArray, val, ind::Integer)
    x.refs[ind] = getpoolidx(x, val)
    return x
end

##############################################################################
##
## growing and shrinking
##
##############################################################################

function Base.push!(pv::PooledVector{S,R}, v::T) where {S,R,T}
    v = convert(S,v)
    push!(pv.refs, getpoolidx(pv, v))
    return v
end

Base.pop!(pv::PooledVector) = pv.pool[pop!(pv.refs)]

function Base.unshift!(pv::PooledVector{S,R}, v::T) where {S,R,T}
    v = convert(S,v)
    pushfirst!(pv.refs, getpoolidx(pv, v))
    return v
end

Base.popfirst!(pv::PooledVector) = pv.pool[popfirst!(pv.refs)]

Base.empty!(pv::PooledVector) = (empty!(pv.refs); pv)

function Base.vcat(a::PooledArray{S}, b::Array{T}) where {S,T}
    output = Array{promote_type(S, T)}(length(b) + length(a))
    copy!(output, 1, a, 1, length(a))
    copy!(output, length(a)+1, b, 1, length(b))
end

function Base.vcat(a::Array{T}, b::PooledArray{S}) where {T,S}
    output = Array{promote_type(S, T)}(length(b) + length(a))
    copy!(output, 1, a, 1, length(a))
    copy!(output, length(a)+1, b, 1, length(b))
end

function Base.vcat(a::PooledArray, b::PooledArray)
    ap = a.pool
    bp = b.pool

    poolmap = Dict{Int, Int}()
    l = length(ap)
    for (i, x) in enumerate(bp)
        if x in ap
            poolmap[i] = findfirst(ap, x)
        else
            poolmap[i] = (l+=1)
        end
    end
    newpool = Array{promote_type(eltype(a), eltype(b))}(l)
    for i = 1:length(ap)
        newpool[i] = ap[i]
    end
    invmap = Dict(v=>k for (k, v) in poolmap)
    for i = length(ap)+1:l
        newpool[i] = bp[invmap[i]]
    end
    refs2 = map(r->poolmap[r], b.refs)
    types = [UInt8, UInt16, UInt32, UInt]
    tidx = findfirst(t->l < typemax(t), types)
    T = types[tidx]
    newrefs = Base.typed_vcat(T, a.refs, refs2)
    return PooledArray(RefArray(newrefs), newpool)
end

end
