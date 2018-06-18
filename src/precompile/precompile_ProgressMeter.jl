function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(Tuple{typeof(ProgressMeter.durationstring), Float64})
    precompile(Tuple{typeof(ProgressMeter.durationstring), Int64})
    precompile(Tuple{getfield(ProgressMeter, Symbol("##barstring#17")), ProgressMeter.BarGlyphs, typeof(identity), Int64, Float64})
    precompile(Tuple{getfield(ProgressMeter, Symbol("##updateProgress!#5")), Array{Any, 1}, Symbol, Int, ProgressMeter.Progress})
    precompile(Tuple{typeof(ProgressMeter.move_cursor_up_while_clearing_lines), Base.TTY, Int64})
    precompile(Tuple{typeof(ProgressMeter.printover), Base.TTY, String, Symbol})
    precompile(Tuple{getfield(ProgressMeter, Symbol("##printvalues!#16")), Symbol, typeof(identity), ProgressMeter.Progress, Array{Any, 1}})
    precompile(Tuple{typeof(ProgressMeter.tty_width), String})
    precompile(Tuple{getfield(ProgressMeter, Symbol("##Progress#1#2")), Int64, String, Symbol, Base.TTY, Int64, ProgressMeter.BarGlyphs, Type{Int}, Int64})
end
