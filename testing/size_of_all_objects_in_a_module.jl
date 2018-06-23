
objs=names(Main)
objs_sizes=Vector{Int}(length(objs))
for x=1:length(objs)
    @show x
    objs_sizes[x]= Base.summarysize(Core.eval(objs[x]))
end
srt=sortperm(objs_sizes)
objs_sizes_srt=objs_sizes[srt]
objs_sizes_srt_mb=objs_sizes_srt./1024./1024
objs_srt=objs[srt]
rs=hcat(string.(objs_srt),string.(objs_sizes_srt_mb))

import Core 
for x in fieldnames(typeof(dtm1model))
   zi= Core.eval("Base.summarysize(dtm2model.$(string(x)))))
   @show x,zi
end)))"