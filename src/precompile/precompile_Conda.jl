function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(Tuple{typeof(Conda.szExeFile), Conda.PROCESSENTRY32})
    precompile(Tuple{typeof(Conda.isrunning), String})
    precompile(Tuple{typeof(Conda._install_conda), String, Bool})
    precompile(Tuple{typeof(Conda._set_conda_env), Base.Cmd, String})
    precompile(Tuple{typeof(Conda.add_channel), String})
    precompile(Tuple{typeof(Conda._installer_url)})
    precompile(Tuple{typeof(Conda.runconda), Base.Cmd})
    precompile(Tuple{typeof(Conda._quiet)})
    precompile(Tuple{typeof(Conda.runconda), Base.Cmd, String})
    precompile(Tuple{typeof(Conda.add_channel), String, String})
    precompile(Tuple{typeof(Conda._install_conda), String})
end
