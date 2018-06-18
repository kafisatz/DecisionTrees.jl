function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(Tuple{typeof(Mmap.gethandle), Base.IOStream})
    precompile(Tuple{getfield(Mmap, Symbol("##mmap#1")), Bool, Bool, typeof(identity), Base.IOStream, Type{Array{UInt8, 1}}, Tuple{Int64}, Int64})
    precompile(Tuple{typeof(Mmap.mmap), Base.IOStream, Type{Array{UInt8, 1}}})
    precompile(Tuple{getfield(Mmap, Symbol("##3#5")), Array{UInt8, 1}})
end
