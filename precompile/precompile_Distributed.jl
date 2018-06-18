function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(Tuple{typeof(Distributed.nprocs)})
    precompile(Tuple{typeof(Distributed.send_del_client), Distributed.Future})
    precompile(Tuple{typeof(Distributed.procs)})
    precompile(Tuple{typeof(Distributed.procs)})
    precompile(Tuple{typeof(Distributed.nprocs)})
    precompile(Tuple{typeof(Distributed._require_callback), Base.PkgId})
    precompile(Tuple{typeof(Distributed._require_callback), Base.PkgId})
    precompile(Tuple{typeof(Distributed.fetch_ref), Distributed.RRID})
end
