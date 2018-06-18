#SnoopComplie

using SnoopCompile

### Log the compiles
# This only needs to be run once (to generate "/tmp/images_compiles.csv")


snoopfile="C:\\temp\\_snoop_dtm_compiles.csv"
SnoopCompile.@snoop snoopfile begin
    fn=joinpath(dirname(@__FILE__),"..","test","runtests.jl")
    #fn=Pkg.dir("DecisionTrees", "test","runtests.jl")
    @assert isfile(fn)
    include(fn)
end

### Parse the compiles and generate precompilation scripts
# This can be run repeatedly to tweak the scripts

data = SnoopCompile.read(snoopfile)

pc = SnoopCompile.parcel(reverse!(data[2]))
SnoopCompile.write("c:\\temp\\precompile", pc)