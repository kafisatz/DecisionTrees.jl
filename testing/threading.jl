# https://docs.julialang.org/en/latest/manual/parallel-computing/#The-@threads-Macro-1

# july 5, 2018
# threading does not seem to improve performance when applied to the innermost get_deviances loop

# ramp up
include(joinpath(Pkg.dir("DecisionTrees"), "testing", "prepare_a_dtmtable.jl"))

using StatsBase

function get_deviances_threaded(a::DecisionTrees.PoissonDevianceSplit, current_meanl::Float64, current_meanr::Float64, lo, ooo, f, numerator::Array{Float64,1}, denominator::Array{Float64,1}, weight::Array{Float64,1}, elementsInLeftChildBV)
    # for the poisson deviance consider the derivative of the poisson loss (or google, the reacfin paper or the 'axa' master thesis)
	# this is the core function of the modelling process
	# besides copying of the data, the vast majority of time is spent in here!
	# most of the time is spent here, if we can improve the for loop below, that would improve performance greatly!
	# one possibility would be to introduce parallelization here (which is not straightforward, though
	
	# dr and dl are 'reused' by each iteration
	dr = Threads.Atomic{Float64}(0.0)
	dl = Threads.Atomic{Float64}(0.0)
	# note: inbounds increases efficiency here (about a factor of 2), however if the bounds are violated something nasty might happen (quote: If the subscripts are ever out of bounds, you may suffer crashes or silent corruption.)
	Threads.@threads for count in 1:length(f)	
		@inbounds idx = f.parent.refs[count] + ooo		
		@inbounds ni = numerator[count]
        @inbounds wi = weight[count]
		@inbounds eli = elementsInLeftChildBV[idx]    
        if eli
            if iszero(ni)
                Threads.atomic_add!(dl, -(ni - wi * current_meanl))
            else 
                Threads.atomic_add!(dl, -(ni - wi * current_meanl) - ni * log(wi * current_meanl / ni))
            end
        else            
            if iszero(ni)
                Threads.atomic_add!(dr, -(ni - wi * current_meanr))
            else 
                Threads.atomic_add!(dr, -(ni - wi * current_meanr) - ni * log(wi * current_meanr / ni))
            end        
        end
	end
	return dl[], dr[]
end


function get_deviancesALTERNATIVE(a::DecisionTrees.PoissonDevianceSplit, current_meanl::Float64, current_meanr::Float64, lo, ooo, f, numerator::Array{Float64,1}, denominator::Array{Float64,1}, weight::Array{Float64,1}, elementsInLeftChildBV)
    # for the poisson deviance consider the derivative of the poisson loss (or google, the reacfin paper or the 'axa' master thesis)
	# this is the core function of the modelling process
	# besides copying of the data, the vast majority of time is spent in here!
	# most of the time is spent here, if we can improve the for loop below, that would improve performance greatly!
	# one possibility would be to introduce parallelization here (which is not straightforward, though
	
	# dr and dl are 'reused' by each iteration
	dr = 0.0 # Threads.Atomic{Float64}(0.0)
	dl = 0.0 # Threads.Atomic{Float64}(0.0)
	# note: inbounds increases efficiency here (about a factor of 2), however if the bounds are violated something nasty might happen (quote: If the subscripts are ever out of bounds, you may suffer crashes or silent corruption.)
	for count in 1:length(f)	
		@inbounds idx = f.parent.refs[count] + ooo	
		@inbounds ni = numerator[count]
        @inbounds wi = weight[count]
		@inbounds eli = elementsInLeftChildBV[idx]       
        if eli
            if iszero(ni)
                dl += -(ni - wi * current_meanl)
            else 
                dl += -(ni - wi * current_meanl) - ni * log(wi * current_meanl / ni)
            end
        else            
            if iszero(ni)
                dr += -(ni - wi * current_meanr)
            else 
                dr += -(ni - wi * current_meanr) - ni * log(wi * current_meanr / ni)
            end        
        end
	end
	return dl, dr
end


function check_equality(dataloops, nloops, nsizeNumerator)
# nloops=2#00

    for i = 1:dataloops
        nsize = nsizeNumerator # 15
        denominatoR = rand(nsize)
        numeratoR = rand(nsize)
        weight = rand(nsize)
        trni = sample(1:nsize, floor(Int, .4 * nsize), replace=false)
        sort!(trni)
        psplit = DecisionTrees.PoissonDevianceSplit()

        for thisloop2 = 1:nloops
            strArr = [randstring(1) for i = 1:nsize]
            floatArr = float.(abs.(rand(Int, nsize) .% 99))
            pdStr = DecisionTrees.PooledArray(strArr)
            pdFloat = DecisionTrees.PooledArray(floatArr)

            if rand() < .5
                pd = pdStr 
            else 
                pd = pdFloat
            end

            fpd = view(pd, trni)
            lo = one(eltype(fpd.parent.refs)) # UInt8(1)
            ooo = one(lo) - lo
            uqrefsmod = sort(unique(fpd.parent.refs)) .+ ooo
            elLBV = BitVector(rand(length(uqrefsmod)) .< rand())
            current_meanl = rand()
            current_meanr = rand()

            resl, resr = DecisionTrees.get_deviances(psplit, current_meanl, current_meanr, lo, ooo, fpd, numeratoR, denominatoR, weight, elLBV)
            resl2, resr2 =    get_deviances_threaded(psplit, current_meanl, current_meanr, lo, ooo, fpd, numeratoR, denominatoR, weight, elLBV)
            resl3, resr3 = get_deviancesALTERNATIVE(psplit, current_meanl, current_meanr, lo, ooo, fpd, numeratoR, denominatoR, weight, elLBV)
            @assert isapprox(resl, resl2)
            @assert isapprox(resl, resl3)
            @assert isapprox(resr, resr2)
            @assert isapprox(resr, resr3)
        end 
    end 

    return nothing

end


function timeOLD(dataloops, nloops, nsizeNumerator)
    # nloops=2#00
    
    for i = 1:dataloops
        nsize = nsizeNumerator # 15
        denominatoR = rand(nsize)
        numeratoR = rand(nsize)
        weight = rand(nsize)
        trni = sample(1:nsize, floor(Int, .4 * nsize), replace=false)
        sort!(trni)
        psplit = DecisionTrees.PoissonDevianceSplit()
    
        for thisloop2 = 1:nloops
            strArr = [randstring(1) for i = 1:nsize]
            floatArr = float.(abs.(rand(Int, nsize) .% 99))
            pdStr = DecisionTrees.PooledArray(strArr)
            pdFloat = DecisionTrees.PooledArray(floatArr)
    
            if rand() < .5
                pd = pdStr 
            else 
                pd = pdFloat
            end
    
            fpd = view(pd, trni)
            lo = one(eltype(fpd.parent.refs)) # UInt8(1)
            ooo = one(lo) - lo
            uqrefsmod = sort(unique(fpd.parent.refs)) .+ ooo
            elLBV = BitVector(rand(length(uqrefsmod)) .< rand())
            current_meanl = rand()
            current_meanr = rand()
    
            resl, resr = DecisionTrees.get_deviances(psplit, current_meanl, current_meanr, lo, ooo, fpd, numeratoR, denominatoR, weight, elLBV)
            # resl2,resr2=    get_deviances_threaded(psplit,current_meanl,current_meanr,lo,ooo,fpd,numeratoR,denominatoR,weight,elLBV)
            # @assert isapprox(resl,resl2)
            # @assert isapprox(resr,resr2)
        end 
    end 
    
    return nothing
    
end

    
    
function timeNEW(dataloops, nloops, nsizeNumerator)
    # nloops=2#00
    
    for i = 1:dataloops
        nsize = nsizeNumerator # 15
        denominatoR = rand(nsize)
        numeratoR = rand(nsize)
        weight = rand(nsize)
        trni = sample(1:nsize, floor(Int, .4 * nsize), replace=false)
        sort!(trni)
        psplit = DecisionTrees.PoissonDevianceSplit()
    
        for thisloop2 = 1:nloops
            strArr = [randstring(1) for i = 1:nsize]
            floatArr = float.(abs.(rand(Int, nsize) .% 99))
            pdStr = DecisionTrees.PooledArray(strArr)
            pdFloat = DecisionTrees.PooledArray(floatArr)
    
            if rand() < .5
                pd = pdStr 
            else 
                pd = pdFloat
            end
    
            fpd = view(pd, trni)
            lo = one(eltype(fpd.parent.refs)) # UInt8(1)
            ooo = one(lo) - lo
            uqrefsmod = sort(unique(fpd.parent.refs)) .+ ooo
            elLBV = BitVector(rand(length(uqrefsmod)) .< rand())
            current_meanl = rand()
            current_meanr = rand()
    
            # resl,resr=DecisionTrees.get_deviances(psplit,current_meanl,current_meanr,lo,ooo,fpd,numeratoR,denominatoR,weight,elLBV)
            resl2, resr2 =    get_deviances_threaded(psplit, current_meanl, current_meanr, lo, ooo, fpd, numeratoR, denominatoR, weight, elLBV)
            # @assert isapprox(resl,resl2)
            # @assert isapprox(resr,resr2)
        end 
    end 
    
    return nothing
    
end


    
function timeALT(dataloops, nloops, nsizeNumerator)
        # nloops=2#00
        
    for i = 1:dataloops
        nsize = nsizeNumerator # 15
        denominatoR = rand(nsize)
        numeratoR = rand(nsize)
        weight = rand(nsize)
        trni = sample(1:nsize, floor(Int, .4 * nsize), replace=false)
        sort!(trni)
        psplit = DecisionTrees.PoissonDevianceSplit()
        
        for thisloop2 = 1:nloops
            strArr = [randstring(1) for i = 1:nsize]
            floatArr = float.(abs.(rand(Int, nsize) .% 99))
            pdStr = DecisionTrees.PooledArray(strArr)
            pdFloat = DecisionTrees.PooledArray(floatArr)
        
            if rand() < .5
                pd = pdStr 
            else 
                pd = pdFloat
            end
        
            fpd = view(pd, trni)
            lo = one(eltype(fpd.parent.refs)) # UInt8(1)
            ooo = one(lo) - lo
            uqrefsmod = sort(unique(fpd.parent.refs)) .+ ooo
            elLBV = BitVector(rand(length(uqrefsmod)) .< rand())
            current_meanl = rand()
            current_meanr = rand()
        
                # resl,resr=DecisionTrees.get_deviances(psplit,current_meanl,current_meanr,lo,ooo,fpd,numeratoR,denominatoR,weight,elLBV)
            resl2, resr2 =    get_deviancesALTERNATIVE(psplit, current_meanl, current_meanr, lo, ooo, fpd, numeratoR, denominatoR, weight, elLBV)
                # @assert isapprox(resl,resl2)
                # @assert isapprox(resr,resr2)
        end 
    end 
        
    return nothing
        
end
    

        
@time check_equality(10, 20, 200000);

@time timeOLD(10, 20, 200000);
@time timeNEW(10, 20, 200000); # 50% more than old (with 1 thread only!)
@time timeALT(10, 20, 200000); # same as old

using BenchmarkTools

@benchmark get_deviances_threaded(psplit, current_meanl, current_meanr, lo, ooo, fpd, numeratoR, denominatoR, weight, elLBV)
@benchmark get_deviancesALTERNATIVE($psplit, $current_meanl, $current_meanr, $lo, $ooo, $fpd, $numeratoR, $denominatoR, $weight, $elLBV)
@benchmark DecisionTrees.get_deviances($psplit, $current_meanl, $current_meanr, $lo, $ooo, $fpd, $numeratoR, $denominatoR, $weight, $elLBV)

@code_warntype get_deviancesALTERNATIVE(psplit, current_meanl, current_meanr, lo, ooo, fpd, numeratoR, denominatoR, weight, elLBV)
@code_warntype DecisionTrees.get_deviances(psplit, current_meanl, current_meanr, lo, ooo, fpd, numeratoR, denominatoR, weight, elLBV)
@code_warntype get_deviances_threaded(psplit, current_meanl, current_meanr, lo, ooo, fpd, numeratoR, denominatoR, weight, elLBV)