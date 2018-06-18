function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(Tuple{getfield(JLD2, Symbol("##s36#29")), Type{Int}, Type{Int}, Type{Int}})
    precompile(Tuple{typeof(JLD2.hasdata), DataType})
    precompile(Tuple{typeof(JLD2.odr_sizeof), JLD2.OnDiskRepresentation{(0, 8, 16, 24, 32, 40, 48, 56, 64), Tuple{Array{String, 1}, Array{Int64, 1}, Array{Int64, 1}, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}, DataFrames.DataFrame, Array{Array{Float64, 1}, 1}, Array{Array{String, 1}, 1}}, Tuple{JLD2.RelOffset, JLD2.RelOffset, JLD2.RelOffset, JLD2.RelOffset, JLD2.RelOffset, JLD2.RelOffset, JLD2.RelOffset, JLD2.RelOffset, JLD2.RelOffset}}})
    precompile(Tuple{typeof(JLD2.samelayout), DataType})
    precompile(Tuple{typeof(JLD2.odr_sizeof), DataType})
    precompile(Tuple{typeof(JLD2.hasfielddata), Type{Int}})
    precompile(Tuple{getfield(JLD2, Symbol("##s36#30")), Type{Int}, Type{Int}, Type{Int}})
    precompile(Tuple{typeof(JLD2.odr_sizeof), JLD2.OnDiskRepresentation{(0, 16), Tuple{String, Array{Any, 1}}, Tuple{JLD2.Vlen{String}, JLD2.Vlen{JLD2.RelOffset}}}})
    precompile(Tuple{getfield(JLD2, Symbol("##s36#21")), Type{Int}, Type{Int}, Type{Int}, Type{Int}})
    precompile(Tuple{typeof(JLD2.ismutabletype), DataType})
    precompile(Tuple{typeof(JLD2.writeas), Type{Int}})
end
