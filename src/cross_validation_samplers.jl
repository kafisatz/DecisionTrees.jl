#from https://github.com/JuliaStats/MLBase.jl/blob/master/src/crossval.jl 
import Base: start, next, done,length

abstract type CrossValGenerator end

struct Kfold <: CrossValGenerator
    permseq::Vector{Int}
    k::Int
    coeff::Float64

    function Kfold(n::Int, k::Int)
        2 <= k <= n || error("DTM: (Kfold CV Sampler): The value of k must be in [2, length(a)]. You provided k=$(k)")
        new(randperm(n), k, n / k)
    end
end

length(c::Kfold) = c.k

struct KfoldState
    i::Int      # the i-th of the subset 
    s::Int      # starting index
    e::Int      # ending index
end

start(c::Kfold) = KfoldState(1, 1, round.(Integer,c.coeff))
next(c::Kfold, s::KfoldState) = (i = s.i+1; (setdiff(1:length(c.permseq), c.permseq[s.s:s.e]), KfoldState(i, s.e+1, round.(Integer,c.coeff * i))))
done(c::Kfold, s::KfoldState) = (s.i > c.k) 
@inline iterate(x::Kfold) = next(x, start(x))
@inline iterate(x::Kfold, i) = done(x, i) ? done : next(x, i)

# Repeated random sub-sampling

struct RandomSub <: CrossValGenerator
    n::Int    # total length
    sn::Int   # length of each subset
    k::Int    # number of subsets
end

length(c::RandomSub) = c.k

start(c::RandomSub) = 1
next(c::RandomSub, s::Int) = (sort!(StatsBase.sample(1:c.n, c.sn; replace=false)), s+1)
done(c::RandomSub, s::Int) = (s > c.k)
@inline iterate(x::RandomSub) = next(x, start(x))
@inline iterate(x::RandomSub, i) = done(x, i) ? done : next(x, i)

#Kfold, but disjoing sets ('negation' of Kfold sets)
struct KfoldDisjoint <: CrossValGenerator
    permseq::Vector{Int}
    k::Int
    coeff::Float64

    function KfoldDisjoint(n::Int, k::Int)
        2 <= k <= n || error("DTM: (KfoldDisjoint CV Sampler): The value of k must be in [2, length(a)]. You provided k=$(k)")
        #2 <= k <= n || error("The value of k must be in [2, length(a)].")
        new(randperm(n), k, n / k)
    end
end

length(c::KfoldDisjoint) = c.k

start(c::KfoldDisjoint) = KfoldState(1, 1, round.(Integer,c.coeff))
next(c::KfoldDisjoint, s::KfoldState) =
    (i = s.i+1; (sort(c.permseq[s.s:s.e]), KfoldState(i, s.e+1, round.(Integer,c.coeff * i))))
done(c::KfoldDisjoint, s::KfoldState) = (s.i > c.k) 
@inline iterate(x::KfoldDisjoint) = next(x, start(x))
@inline iterate(x::KfoldDisjoint, i) = done(x, i) ? done : next(x, i)

#=

#testing 

collect(KfoldDisjoint(10,3))
collect(RandomSub(10,4,20))
collect(Kfold(10,3))

    nn=10
    kk=3
    sampled_kfold=Kfold(nn,kk)
    for i in sampled_kfold
        @show i    
        @show setdiff(1:nn,i)
    end
=#