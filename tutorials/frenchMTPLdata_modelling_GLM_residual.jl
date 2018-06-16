###########


t0=time_ns()
cd(joinpath(GLOBAL_julia_code_folder,"DecisionTrees.jl"))
#@info("You may want to run 'pkg> instantiate' when you first run this. Use ] to enter the package mode.")

@everywhere using Revise
@everywhere using CSV,DataFrames

@everywhere using DecisionTrees  

tela = (time_ns()-t0)/1e9
@show tela #precompilationcan take considerable time under Julia 0.7alpha (up to 15 minutes in some cases (as ofJune 11, 2018)), this should improve soon though 


##############################
#Read the data
##############################
datafile=string("data\\freMTPL2\\freMTPL2.csv")
@assert isfile(datafile);
@time fullData=CSV.read(datafile,rows_for_type_detect=100000,allowmissing=:none,categorical=false);
