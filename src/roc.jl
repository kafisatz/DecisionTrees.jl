# from https://github.com/davidavdav/ROCAnalysis.jl/blob/master/src/main.jl

# Note: this file may cause issues when copy/pasted to julia with autohotkey
# it is best to include it include("thisfile") and then test the functions

#Copyright (c) 2014--2015 David A. van Leeuwen
#The Julia package "ROC.jl" is licenced under the MIT License.
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## roc.jl  Receiver Operating Characteristics computations
## (c) 2013--2015 David A. van Leeuwen
##
## Licensed under the MIT software license, see LICENSE.md
## Routines to compute Receiver Operating Characteristics

## These routines need a re-write, we can probably unroll a few loops
## and be more efficient with memory and sorting.

## Roc: essential data to store Receiver Operating Characteristics
"""
`Roc` is a type that stores the essential performance information that can be extracted from a 
set of supervised trials, i.e., target and non-target scores from a two-class classifier.  Apart from
the (minimalized) arrays for probabilities of false-alarm and miss---the coordinates of the ROC 
curve---, they are the threshold for these, a boolean whether-or-not this point lies on the 
convex hull, and the associated optimal log-likelihood-ratio associated to the line segment. 
"""    
mutable struct Roc{T<:Real}
    pfa::Vector{Float64}        # probability of false alarm
    pmiss::Vector{Float64}      # probability of miss
    ch::BitArray                # point lies on convex hull
    theta::Vector{T}                # threshold
    llr::Vector{Float64}        # optimal log likelihood ratio
end

"""
`roc(tar, non; laplace, collapse)` computes the essential statistics for evaluation of the 
performance of a two-class classifier.  The true class of the scores is encoded in the 
array in which they appear, i.e., `tar`get or `non`-target scores. 

The results are given in a type `Roc`, containing the receiver-operating-characteristics (ROC) data
for this test.  

Optional arguments:
 
 - `collapse=true`, indicating whether the resulting arrays for the probability of 
false positives (alarms) and negatives (misses) should be collapsed, i.e., consecutive points 
on the ROC that have the same false positive or false negatove rate are removed, retaining 
only the `corner' points of the ROC. 

 - `laplace=false`, indicating whether two additional data points at either end of the score scale 
should be added to the target an dnon-target scores.  This corresponds to the Laplace prior, and
has a result of limiting the magnitude of the optimal log-likelihood-ratio associated with the
line segments of the convex hull of the ROC. 
"""
function roc(tar::Vector{T}, non::Vector{T}; laplace::Bool=false, collapse=true) where {T <: Real} 
    xo, tc = sortscores(tar, non)
    ## first collect point of the same threshold (for discrete score data)
    theta, tc, nc, = binscores(xo, tc)
    ## now we can compute pmiss and pfa
    pfa = [1; 1 - cumsum(nc)/length(non)]
    pmiss = [0; cumsum(tc)/length(tar)]
    ## convex hull and llr
    ch, llr = chllr(tc, nc, xo, laplace=laplace)
    if collapse
        ## remove points on the ROC in the middle of a horizontal or vertical segment
        changes = changepoints(pfa, pmiss)
        (pfa, pmiss, ch) = map(a -> a[changes], (pfa, pmiss, ch))
        (llr, theta) = map(a -> a[changes[1:end-1]], (llr, theta))
    end
    Roc(pfa, pmiss, ch, theta, llr)
end

## this returns the target count (1's an 0's for targets and non targets) and scores, 
## in the order of the scores. 
function sortscores(tar::Vector{T}, non::Vector{T}) where {T <: Real}
    x = vcat(tar, non)                     # order targets before non-targets
    so = sortperm(x)                       # sort order
    xo = x[so]                             # scores ordered
    tc = vcat(ones(Int, length(tar)), zeros(Int, length(non)))[so] # target count
    return xo, tc
end    

## bins scores that are the same into aggredated score counts
function binscores(xo::Vector{T},tc::Vector{Int}) where {T <: Real}
    nc = 1 - tc
    changes = findall([true; diff(xo) .!= 0; true]) # points where threshold change
    theta = xo                                       # threshold
    keep = trues(length(tc))
    @inbounds for i=1:length(changes)-1
        start = changes[i]
        stop = changes[i+1]-1
        if stop>start
            tc[start] = sum(tc[start:stop])
            nc[start] = sum(nc[start:stop])
            keep[start+1:stop]= false
        end
    end
    (tc, nc, theta) = map(a -> a[keep], (tc, nc, theta))
    return theta, tc, nc, keep
end
## This finds the change points on the ROC (the corner points)
function changepoints(pfa::Vector{Float64}, pmiss::Vector{Float64})
    (const_pfa, const_pmiss) = map(a -> diff(a) .== 0, (pfa, pmiss)) # no changes
    return [true; ! (const_pfa[1:end-1]  &  const_pfa[2:end] | 
                     const_pmiss[1:end-1] & const_pmiss[2:end]); true]
end

## compute convex hull and optimal llr from target count, nontarget count, and ordered scores
function chllr(tc::Vector{Int}, nc::Vector{Int}, xo::Vector{T}; laplace::Bool=true) where {T}
    if laplace
        tc = [1; 0; tc; 1; 0]
        nc = [0; 1; nc; 0; 1]
        xo = [-Inf; -Inf; xo; Inf; Inf]
    end
    ntar = sum(tc)
    nnon = sum(nc)
#    if fast
        pfa = [1; 1 - cumsum(nc)/nnon] # use rationals for accuracy of the convex hull
        pmiss = [0; cumsum(tc)/ntar]
        index = rochull(pfa, pmiss)
#    else
#        pfa = [1, 1 - cumsum(nc)/nnon]
#        pmiss = [0, cumsum(tc)/ntar]
        ## convex hull
#        hull = chull(vcat(hcat(pfa, pmiss), [2 2])) # convex hull points
#        index = sort(hull.vertices[:,1])[1:end-1]   # indices of the points on the CH
#    end
    ch = falses(length(pfa))
    ch[index] = true
    ## LLR
    (dpfa, dpmiss) = map(a -> diff(a[ch]), (pfa, pmiss))
    llr = log(abs(dpmiss ./ dpfa))[cumsum(ch)[1:end-1]]
    if laplace
        return [true; ch[4:end-3]; true], llr[3:end-2]
    else
        return ch, llr
    end
end

## returns positive number if test point is "left" of line (x1,y1) -- (x2,y2)
## positive: left, negative: right, zero: on
## This works best with rationals, but that makes it slow 
function isleft(x1::T, y1::T, x2::T, y2::T, xt::T, yt::T) where {T <: Real}
    (x2 - x1)*(yt - y1) - (xt - x1)*(y2 - y1) 
end

## Andrew's Monotone Chain Algorithm, from http://geomalgorithms.com/a10-_hull-1.html
## (ignoring the C-code:-)
## points are (pmiss, pfa), pmiss (x) is increasing, pfa (y) is decreasing
## we implictly add a point at (2,2) that is on the convex hull
function rochull(pfa::Vector{T}, pmiss::Vector{T}) where {T <: Real}
    two = 2one(T)
    n = length(pmiss)
    ## minmin and minmax
    i = 1
    while i < n
        if pmiss[i+1] == pmiss[1]
            i += 1
        else
            break
        end
    end
    minmin = i
    minmax = 1
    ## maxmin and maxmax
    maxmin = maxmax = n+1       # virtual point for now
    stack = [minmin]
    for i=minmin+1:n
        if isleft(pmiss[minmin], pfa[minmin], two, two, pmiss[i], pfa[i]) >= zero(T)
            continue
        end
        while length(stack) >= 2
            if isleft(pmiss[stack[end-1]], pfa[stack[end-1]],
                      pmiss[stack[end]], pfa[stack[end]], pmiss[i], pfa[i]) > zero(T)
                break
            else
                pop!(stack)
            end
        end
        push!(stack, i)
    end
    pushfirst!(stack, minmax)
    stack
end

## some support for standard Julia functions
Base.length(r::Roc) = length(r.pfa)
Base.size(r::Roc) = (length(r),)
function Base.size(r::Roc, d::Integer)
    d==1 || error("Roc is a 1-dimensional structure")
    return length(r)
end
start(r::Roc) = 1
done(r::Roc, state) = state>length(r)
endof(r::Roc) = length(r)
Base.next(r::Roc, state) = (r[state], state+1)
@inline iterate(x::Roc) = next(x, start(x))
@inline iterate(x::Roc, i) = done(x, i) ? done : next(x, i)

function Base.getindex(r::Roc, i::Integer)
    if i == length(r)
        return (r.pfa[i], r.pmiss[i], Inf, r.ch[i], Inf)
    else
        return (r.pfa[i], r.pmiss[i], r.theta[i], r.ch[i], r.llr[i])
    end
end

## R terminology, quantile function qnorm() and cumulative distribution pnorm()
# qnorm(x) = v2 * erfinv(2x-1)
# pnorm(x) = (1 + erf(x/v2)) / 2 

## The probability that a random non-target score is higher than a target score
auc(r::Roc) = -dot(r.pmiss[1:end-1] + r.pmiss[2:end], diff(r.pfa)) / 2
auc(tar::Vector, non::Vector) = auc(roc(tar, non))

#=
	#test
	x=[.1 .91 0.1 0 0][:]
	y=[1.0 1 0 1 0][:]
	rr=roc(x,y)
	a=auc(x,y)
=#