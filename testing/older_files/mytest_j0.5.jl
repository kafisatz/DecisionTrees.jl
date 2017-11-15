this_dir=string(dictGlobalSettings["Julia Code Folder"],"\\src");cd(this_dir);push!(LOAD_PATH,this_dir)
using DTM

#this runs a small version to compile all functions
chosenj="boosting_small";ARGS=[string("a:\\JuliaTestData\\",chosenj,".settings.csv") string("a:\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];
resbool=DTM.run_model(ARGS) #small version to compile functions


using Gallium
include(Pkg.dir("Gallium/examples/testprograms/misc.jl"))
Gallium.breakpoint(sinthesin,Tuple{Int64})
inaloop(2)

include(string(dictGlobalSettings["Julia Code Folder"],"src\\DTM.jl"))
import DTM.derive_scores
Gallium.breakpoint(derive_scores) #,Array{Float64,1},Array{Float64,1})

#use ? for help
#use `res (or `size(res)) to evaluate "something" in the current scope

chosenj="julia200k_100s";ARGS=[string("a:\\JuliaTestData\\",chosenj,".settings.csv") string("a:\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];
resbool=DTM.run_model(ARGS)

include("R:\\TariffWatch New Repo\\algorithms\\Julia\\Code\\dev_julia_0.5\\testing\\perf_test_boosting.jl")
