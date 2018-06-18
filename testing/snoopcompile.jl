#SnoopComplie

using SnoopCompile

### Log the compiles
# This only needs to be run once (to generate "/tmp/images_compiles.csv")
SnoopCompile.@snoop "C:\\temp\\__snoop_dtm_compiles.csv" begin
    include(Pkg.dir("DecisionTrees", "test","runtests.jl")    )
end


### Parse the compiles and generate precompilation scripts
# This can be run repeatedly to tweak the scripts
data = SnoopCompile.read("C:\\temp\\__snoop_dtm_compiles.csv")

pc = SnoopCompile.parcel(reverse!(data[2]))
SnoopCompile.write("c:\\temp\\precompile", pc)