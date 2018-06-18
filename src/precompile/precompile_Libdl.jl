function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(Tuple{typeof(Libdl.dlsym), Ptr{Nothing}, Symbol})
end
