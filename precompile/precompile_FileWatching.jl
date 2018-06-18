function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(Tuple{typeof(FileWatching.start_watching), FileWatching.FileMonitor})
    precompile(Tuple{typeof(FileWatching.start_watching), FileWatching.PollingFileWatcher})
    precompile(Tuple{typeof(FileWatching.poll_file), String, Float64, Int64})
    precompile(Tuple{typeof(FileWatching.poll_file), String, Float64, Int64})
    precompile(Tuple{typeof(FileWatching.watch_file), String, Int64})
    precompile(Tuple{typeof(FileWatching.watch_file), String, Int64})
end
