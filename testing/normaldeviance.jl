#todo

#test the merge and unmerge functions, are these working correctly for all scenarios?
#check if OnlineStats was 'materially' upgraded (is there something new/faster than merge/unmerge)
#eventually create a tree and check if it is correct (compare to rpart maybe?)
using Test
using BenchmarkTools 
import StatsBase
import OnlineStats
#function lrIndices(trnidx::Vector{Int},f::DecisionTrees.PooledArraysDTM.PooledArray{String,T,1},subset::Array) where T<:Unsigned
include("c:\\temp\\2.jl") #PDA module
using .PooledArraysDTM

import Random


function test_mean_and_var_merge(nloops,sizeofVectors)
    for loopV=1:nloops
        #nn=150
        nn=sizeofVectors
        weight=rand(nn)
        denom=rand(nn)
        num=rand(nn)
        trni=StatsBase.sample(1:nn,floor(Int,.4*nn),replace=false)
        sort!(trni)
        strArr=[randstring(1) for i=1:nn]
        floatArr=float.(abs.(rand(Int,nn).%9))
        pdStr=PooledArray(strArr)
        pdFloat=PooledArray(floatArr)

        if rand()<.5
            pd=pdStr
        else    
            pd=pdFloat
        end 

        thispd=view(pd,trni)
        cnt,sumnumerator,sumdenominator,sumweight,moments_per_pdaclass=aggregate_data_msePointwise(thispd,num,denom,weight)
        #check if the 'merge!' of all compoments matches the total
        total=Series(OnlineStats.Mean(),OnlineStats.Variance())
        for i=1:length(moments_per_pdaclass)
            if moments_per_pdaclass[i].stats[1].n >0 #this is needed because of https://github.com/joshday/OnlineStats.jl/issues/125
                merge!(total,moments_per_pdaclass[i])
            end
        end
        ratioOfND=num./denom
        #@show ratioOfND[trni]
        #@show total 
        #@show mean(ratioOfND[trni]),total.stats[1].μ
        #@show StatsBase.var(ratioOfND[trni]),total.stats[2].σ2
        if length(trni)>0
            @test isapprox(mean(view(ratioOfND,trni)),total.stats[1].μ)
            @test isapprox(StatsBase.var(view(ratioOfND,trni),corrected=false),total.stats[2].σ2)
        end
    end
return nothing 
end

function test_mean_and_var_UNmerge(nloops,sizeofVectors)
        #sizeofVectors=10
        for jj=1:nloops
        
        nn=sizeofVectors
        weight=rand(nn)
        denom=rand(nn)
        num=rand(nn)
        trni=StatsBase.sample(1:nn,floor(Int,.4*nn),replace=false)
        sort!(trni)
        strArr=[randstring(1) for i=1:nn]
        floatArr=float.(abs.(rand(Int,nn).%9))
        pdStr=PooledArray(strArr)
        pdFloat=PooledArray(floatArr)

        if rand()<.5
            pd=pdStr
        else    
            pd=pdFloat
        end 

        thispd=view(pd,trni)
        cnt,sumnumerator,sumdenominator,sumweight,moments_per_pdaclass=aggregate_data_msePointwise(thispd,num,denom,weight)
        #check if the 'merge!' of all compoments matches the total
        total=Series(OnlineStats.Mean(),OnlineStats.Variance())
        for i=1:length(moments_per_pdaclass)
            if moments_per_pdaclass[i].stats[1].n >0 #this is needed because of https://github.com/joshday/OnlineStats.jl/issues/125
                merge!(total,moments_per_pdaclass[i])
            end
        end
        
        if length(moments_per_pdaclass)>2
        #remove one or two random elements of moments_per_pdaclass
        ratioOfND=num./denom
        pick=rand()>.5
        reducedTotal=deepcopy(total)
        if pick
            el1=StatsBase.sample(1:length(moments_per_pdaclass))
            unmerge!(reducedTotal,moments_per_pdaclass[el1])
            valuesRemaining=ratioOfND[trni][(pd[trni].!=pd.pool[el1])]
        else 
            (el1,el2)=StatsBase.sample(1:length(moments_per_pdaclass),2,replace=false)
            unmerge!(reducedTotal,moments_per_pdaclass[el1])
            unmerge!(reducedTotal,moments_per_pdaclass[el2])  
            idx1=(pd.refs.!=pd.pool[2]).&(pd.refs.!=pd.pool[4])
            valuesRemaining=ratioOfND[trni][(pd[trni].!=pd.pool[el1]).&(pd[trni].!=pd.pool[el2])]
        end   
            #=
            @show valuesRemaining
            @show reducedTotal 
            @show mean(valuesRemaining),reducedTotal.stats[1].μ
            @show StatsBase.var(valuesRemaining),reducedTotal.stats[2].σ2
            =#
            if length(trni)>0
                @test isapprox(mean(valuesRemaining),reducedTotal.stats[1].μ)
                @test isapprox(StatsBase.var(valuesRemaining,corrected=false),reducedTotal.stats[2].σ2,rtol=1e-4)
            #    @show eltype(pd),pick
            end
        end
      end
      return nothing 

      end
Random.srand(99231)
test_mean_and_var_merge(1000,1000)
test_mean_and_var_UNmerge(1000,1000)
#@benchmark test_mean_and_var_UNmerge(10,100)
#@benchmark test_mean_and_var_merge(10,1000)
@benchmark cnt,sumnumerator,sumdenominator,sumweight,moments_per_pdaclass=aggregate_data_msePointwise(thispd,num,denom,weight)