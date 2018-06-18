function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(Tuple{typeof(Dates.format), Base.GenericIOBuffer{Array{UInt8, 1}}, Dates.DatePart{Char(0x59000000)}, Dates.DateTime})
    precompile(Tuple{typeof(Dates.format), Base.GenericIOBuffer{Array{UInt8, 1}}, Dates.DatePart{Char(0x73000000)}, Dates.DateTime})
    precompile(Tuple{typeof(Dates.character_codes), Core.SimpleVector})
    precompile(Tuple{getfield(Dates, Symbol("##s623#35")), Int, Int, Int, Int, Int, Int, Int, Int})
    precompile(Tuple{getfield(Dates, Symbol("##s624#32")), Int, Int, Int, Int, Int, Int})
    precompile(Tuple{typeof(Dates.format), Base.GenericIOBuffer{Array{UInt8, 1}}, Dates.DatePart{Char(0x6d000000)}, Dates.DateTime})
    precompile(Tuple{typeof(Dates.format), Base.GenericIOBuffer{Array{UInt8, 1}}, Dates.DatePart{Char(0x64000000)}, Dates.DateTime})
    precompile(Tuple{typeof(Dates.character_codes), Type{Dates.DateFormat{Symbol("yyyy-mm-dd"), Tuple{Dates.DatePart{Char(0x79000000)}, Dates.Delim{Char, 1}, Dates.DatePart{Char(0x6d000000)}, Dates.Delim{Char, 1}, Dates.DatePart{Char(0x64000000)}}}}})
    precompile(Tuple{typeof(Dates._directives), Type{Dates.DateFormat{Symbol("yyyy-mm-dd"), Tuple{Dates.DatePart{Char(0x79000000)}, Dates.Delim{Char, 1}, Dates.DatePart{Char(0x6d000000)}, Dates.Delim{Char, 1}, Dates.DatePart{Char(0x64000000)}}}}})
    precompile(Tuple{typeof(Dates.genvar), DataType})
end
