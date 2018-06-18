function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(Tuple{getfield(Logging, Symbol("##handle_message#2")), Nothing, (Base.Iterators).Pairs{Union{}, Union{}, Tuple{}, NamedTuple{(), Tuple{}}}, typeof(identity), Logging.ConsoleLogger, (Base.CoreLogging).LogLevel, String, Module, Symbol, Symbol, String, Int64})
    precompile(Tuple{getfield(Logging, Symbol("##handle_message#2")), Int64, (Base.Iterators).Pairs{Union{}, Union{}, Tuple{}, NamedTuple{(), Tuple{}}}, typeof(identity), Logging.ConsoleLogger, (Base.CoreLogging).LogLevel, String, Module, Symbol, Symbol, String, Int64})
    precompile(Tuple{getfield(Logging, Symbol("##handle_message#2")), Nothing, (Base.Iterators).Pairs{Union{}, Union{}, Tuple{}, NamedTuple{(), Tuple{}}}, typeof(identity), Logging.ConsoleLogger, (Base.CoreLogging).LogLevel, String, Nothing, Symbol, Symbol, String, Int64})
    precompile(Tuple{getfield(Logging, Symbol("##handle_message#2")), Int64, (Base.Iterators).Pairs{Symbol, (Base.StackTraces).StackFrame, Tuple{Symbol}, NamedTuple{(:caller,), Tuple{(Base.StackTraces).StackFrame}}}, typeof(identity), Logging.ConsoleLogger, (Base.CoreLogging).LogLevel, String, Module, Symbol, Tuple{Ptr{Nothing}, Symbol}, String, Int64})
    precompile(Tuple{typeof(Logging.default_metafmt), Base.CoreLogging.LogLevel, Module, Symbol, Symbol, String, Int64})
    precompile(Tuple{typeof(Logging.default_metafmt), Base.CoreLogging.LogLevel, Nothing, Symbol, Symbol, String, Int64})
    precompile(Tuple{typeof(Logging.default_metafmt), Base.CoreLogging.LogLevel, Module, Symbol, Tuple{Ptr{Nothing}, Symbol}, String, Int64})
end
