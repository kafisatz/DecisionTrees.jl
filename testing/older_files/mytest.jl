include("C:\\julia\\Algorithms\\Julia\\Code\\dev_julia_0.5\\testing\\dict_local_w530.jl")

#=
ENV["PYTHON"]="C:\\Anaconda3\\python.exe"
Pkg.build("PyCall")
=#
this_dir=string(dictGlobalSettings["Julia Code Folder"],"\\src");cd(this_dir);push!(LOAD_PATH,this_dir)
using DTM

#this runs a small version to compile all functions
chosenj="boosting_small";ARGS=[string("a:\\JuliaTestData\\",chosenj,".settings.csv") string("a:\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];
resbool=DTM.run_model(ARGS) #small version to compile functions

chosenj="mx"
ARGS=[string("e:\\temp\\",chosenj,".settings.csv") string("e:\\temp\\",chosenj,".jld") string("out_MMXNE")];

resbool=DTM.run_model(ARGS) #small version to compile functions
