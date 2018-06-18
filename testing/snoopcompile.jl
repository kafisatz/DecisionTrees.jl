#SnoopComplie

using SnoopCompile

### Log the compiles
# This only needs to be run once (to generate "/tmp/images_compiles.csv")



SnoopCompile.@snoop "dtm_compiles.csv" begin
    fn=joinpath(dirname(@__FILE__),"..","test","runtests.jl")
    #fn=Pkg.dir("DecisionTrees", "test","runtests.jl")
    @assert isfile(fn)
    include(fn)
end

### Parse the compiles and generate precompilation scripts
# This can be run repeatedly to tweak the scripts

data = SnoopCompile.read("dtm_compiles.csv")

pc = SnoopCompile.parcel(reverse!(data[2]))
SnoopCompile.write("precompile", pc)