#=
	#in shell 
	sudo update-ca-certificates
	
	sudo apt-get install python3-pip 
	pip3 install pandas
	pip3 install numpy
	pip3 install xlsxwriter
	
=#

./julia1.2/bin/julia

cd("/home/bernhard/DecisionTrees.jl") 
using Pkg
Pkg.activate(pwd())

using Revise 

ENV["PYTHON"]="python3"
import Pkg;println("");println("Building PyCall...");
#Pkg.build("Conda");
try 
	using PyCall;;global const pyModPandas = PyNULL();global const pyModxlsxwriter = PyNULL();copy!(pyModPandas, pyimport_conda("pandas","pandas"));copy!(pyModxlsxwriter, pyimport_conda("xlsxwriter","xlsxwriter"));writer=pyModPandas[:ExcelWriter]("blabla_deletme.xlsx", engine = "xlsxwriter");@show typeof(writer)
catch
	try
	Pkg.build("PyCall");using PyCall;global const pyModPandas = PyNULL();global const pyModxlsxwriter = PyNULL();copy!(pyModPandas, pyimport_conda("pandas","pandas"));copy!(pyModxlsxwriter, pyimport_conda("xlsxwriter","xlsxwriter"));writer=pyModPandas[:ExcelWriter]("blabla_deletme.xlsx", engine = "xlsxwriter");@show typeof(writer);println("If this worked, then pandas and xlsxwriter should be installed");println("Done building PyCall.Restarting Julia...");	
	catch e;println("this did not work");@show e;end; 
	#exit();
end 
using DecisionTrees

t0=time_ns()
import CSV
import DataFrames
import DataFrames: DataFrame
import Random
import OnlineStats
import StatsBase

using Distributed #only for Tests?
using Statistics
using Pkg #only for Tests?
using Test #only for Tests? 
using DecisionTrees

println("Time for using & import: ",round((time_ns()-t0)/1e9,digits=3),"s")

#set directories
global pkgdir=""
global testdir = joinpath(dirname(@__FILE__),"") #joinpath(@__DIR__,"..","test",filename)
#the following if clause is for "direct execution" (e.g. in VScode where @__FILE__ does not seem to work)
if testdir==""
    try 
        global pkgdir=Pkg.dir("DecisionTrees") #how long will this remain supported?
        global testdir=joinpath(pkgdir,"test")
    catch
    end
else 
    global pkgdir=normpath(joinpath(testdir,".."))
end
datadir=joinpath(pkgdir,"data")

Random.seed!(1239946)
#DecisionTreesTests = @testset "DecisionTrees" begin
    t0runtests=time_ns()
    println("DTM: Running tests...") # Time = $(t0runtests)")
        
    global eltypesData1=[Int,	Float64,	Float64,	Float64,	Float64,	Float64,	Int,	String,	String,	String,	Int,	String,	String,	String,	String,	String,	String,	String,	String,	String,	String,	Int,	String,	Int,	String,	String,	Int,	String,	String,	Int,	String,	String,	Int,	String,	String,	Int,	String,	Int,	String,	String,	Int,	Int,	String,	Int,	Int,	String,	Int,	Int,	String,	Int,	String,	String,	String,	Int,	String,	String,	String,	String,	Int,	Int,	Int,	String,	String,	String,	String,	String,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Float64,	Int,	Float64,	Float64,	Float64,	Float64,	Int,	Int,	Int,	Int,	Int,	Int]
    
    @info("DTM: Note that the tests produce a relatively long log and write a number of files in temporary directories. You can disregard the log, if all tests pass.")
    @warn("DTM: So far only very rudimentary tests have been implemented!")
    
datafile=joinpath(datadir,"freMTPL2","freMTPL2.csv")

import Random: rand,randstring
tolForTheseTests=1e-15

    include("test/smoketests.jl")
