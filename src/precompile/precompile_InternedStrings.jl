function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(Tuple{typeof(InternedStrings.intern), Type{String}, WeakRefStrings.WeakRefString{UInt8}})
    precompile(Tuple{typeof(InternedStrings.getvalue), Type{String}, WeakRef})
end
