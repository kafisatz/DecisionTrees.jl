#SnoopComplie

using SnoopCompile

### Log the compiles
# This only needs to be run once (to generate "/tmp/images_compiles.csv")


#snoopfile="C:\\temp\\__snoop_dtm_compiles.csv"
#@assert isdir(splitdir(snoopfile)[1])
#fn=joinpath(dirname(@__FILE__),"..","test","runtests.jl")
#if isfile(fn)
#else 
#    fn=Pkg.dir("DecisionTrees", "test","runtests.jl")    
#end
#@assert isfile(fn)

SnoopCompile.@snoop "C:\\temp\\__snoop_dtm_compiles.csv" begin
    include(Pkg.dir("DecisionTrees", "test","runtests.jl")    )
end

### Parse the compiles and generate precompilation scripts
# This can be run repeatedly to tweak the scripts

data = SnoopCompile.read("C:\\temp\\__snoop_dtm_compiles.csv")

pc = SnoopCompile.parcel(reverse!(data[2]))
SnoopCompile.write("c:\\temp\\precompile", pc)