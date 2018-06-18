function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(Tuple{typeof(Random.shuffle!), Random.MersenneTwister, Array{Symbol, 1}})
    precompile(Tuple{typeof(Random.srand), Random.MersenneTwister, Array{UInt32, 1}})
    precompile(Tuple{typeof(Random.rand!), Random.MersenneTwister, Array{Int64, 1}, Base.UnitRange{Int64}})
    precompile(Tuple{typeof(Random.srand), Array{UInt32, 1}})
end
