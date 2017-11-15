
objs=names(Main)
objs_sizes=Vector{Int}(length(objs))
for x=1:length(objs)
    @show x
    objs_sizes[x]= Base.summarysize(eval(objs[x]))
end
srt=sortperm(objs_sizes)
objs_sizes_srt=objs_sizes[srt]
objs_sizes_srt_mb=objs_sizes_srt./1024./1024
objs_srt=objs[srt]
rs=hcat(string.(objs_srt),string.(objs_sizes_srt_mb))

#objs_sizes=map(Base.summarysize(eval(x)) for x in objs)
#Base.summarysize(eval(objs[end]))

for x in fieldnames(dtm1model)
   zi= eval("Base.summarysize(dtm2model.$(string(x))   )))
   @show x,zi
end)))"