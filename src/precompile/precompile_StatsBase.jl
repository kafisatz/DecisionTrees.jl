function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(Tuple{getfield(StatsBase, Symbol("##s174#130")), Type{Int}, Type{Int}, Tuple{DataType}})
    precompile(Tuple{typeof(StatsBase.corm), Array{Float64, 1}, Float64, Array{Float64, 1}, Float64})
end
