function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(Tuple{getfield(DelimitedFiles, Symbol("##writedlm#15")), (Base.Iterators).Pairs{Union{}, Union{}, Tuple{}, NamedTuple{(), Tuple{}}}, typeof(identity), Base.IOStream, Array{Any, 2}, Char})
    precompile(Tuple{typeof(DelimitedFiles.writedlm_cell), Base.GenericIOBuffer{Array{UInt8, 1}}, String, Char, Bool})
    precompile(Tuple{getfield(DelimitedFiles, Symbol("##writecsv#3")), (Base.Iterators).Pairs{Union{}, Union{}, Tuple{}, NamedTuple{(), Tuple{}}}, typeof(identity), String, Array{Any, 2}})
    precompile(Tuple{typeof(DelimitedFiles.writecsv), String, Array{Any, 2}})
    precompile(Tuple{typeof(DelimitedFiles.writedlm_cell), Base.GenericIOBuffer{Array{UInt8, 1}}, Int64, Char, Bool})
    precompile(Tuple{typeof(DelimitedFiles.writedlm_cell), Base.GenericIOBuffer{Array{UInt8, 1}}, Float64, Char, Bool})
end
