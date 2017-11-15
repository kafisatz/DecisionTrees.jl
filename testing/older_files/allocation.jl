this_dir=string(dictGlobalSettings["Julia Code Folder"],"\\src");cd(this_dir);push!(LOAD_PATH,this_dir)
chosenj="tree_small" #this runs a small version to compile all functions
include(string(this_dir,"\\DTM.jl"))
using DTM
using Graphics,ProfileView

#this runs a small version to compile all functions
loc="C:\\Users\\tbd\\Documents\\Arbeitsordner (sync)\\Public NL\\Personal\\Bernhard\\Resources\\Milliman Software\\Predictive Analytics Paris\\Training\\data\\"
ARGS=[string(loc,"jmortgage.settings_profiling.csv") string(loc,"jmortgage.jld") string("out_jmortgage")];
resbool=DTM.run_model(ARGS) #small version to compile functions	


Profile.clear()          # in case we have any previous profiling data
Profile.init()
#Profile.init(10^7, 0.01) #n,delay #Profiling settings

Profile.clear_malloc_data()
@profile juliaresult=DTM.run_model(ARGS)   # Profile the program


#exit()
#=
this_dir=string(dictGlobalSettings["Julia Code Folder"],"\\src");cd(this_dir);push!(LOAD_PATH,this_dir)
using DTM
#using Graphics,ProfileView

#this runs a small version to compile all functions
chosenj="boosting_small";ARGS=[string("a:\\JuliaTestData\\",chosenj,".settings.csv") string("a:\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];
resbool=DTM.run_model(ARGS) #small version to compile functions	

chosenj="julia200k";ARGS=[string("a:\\JuliaTestData\\",chosenj,".settings.csv") string("a:\\JuliaTestData\\",chosenj,".jld") string("out_",chosenj)];
juliaresult=DTM.run_model(ARGS)   

Profile.clear_malloc_data()
chosenj="julia200k";ARGS=[string("a:\\JuliaTestData\\",chosenj,".settings.csv") string("a:\\JuliaTestData\\",chosenj,".jld") string("out_",chosenj)];
juliaresult=DTM.run_model(ARGS)   
=#

