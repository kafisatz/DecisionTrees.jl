function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(Tuple{typeof(Missings.levels), DecisionTrees.PooledArraysDTM.PooledArray{String, UInt8, 1, Array{UInt8, 1}}})
    precompile(Tuple{typeof(Missings.T), Type{Int64}})
    precompile(Tuple{typeof(Missings.T), Type{Float64}})
    precompile(Tuple{typeof(Missings.T), Type{String}})
    precompile(Tuple{typeof(Missings.levels), DecisionTrees.PooledArraysDTM.PooledArray{Float64, UInt8, 1, Array{UInt8, 1}}})
end
