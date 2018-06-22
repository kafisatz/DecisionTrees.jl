t0=time_ns()

import CSV
import DataFrames
import DataFrames: DataFrame
using Test 
using DecisionTrees 

@show tela = (time_ns()-t0)/1e9

#set directories
global pkgdir=""
global testdir = joinpath(dirname(@__FILE__),"")
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

@testset "DecisionTrees" begin
    println("DTM: Running tests:")
    
    global eltypesData1=[Int,	Float64,	Float64,	Float64,	Float64,	Float64,	Int,	String,	String,	String,	Int,	String,	String,	String,	String,	String,	String,	String,	String,	String,	String,	Int,	String,	Int,	String,	String,	Int,	String,	String,	Int,	String,	String,	Int,	String,	String,	Int,	String,	Int,	String,	String,	Int,	Int,	String,	Int,	Int,	String,	Int,	Int,	String,	Int,	String,	String,	String,	Int,	String,	String,	String,	String,	Int,	Int,	Int,	String,	String,	String,	String,	String,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Float64,	Int,	Float64,	Float64,	Float64,	Float64,	Int,	Int,	Int,	Int,	Int,	Int]
    
    @warn("So far only very rudimentary tests have been implemented!")
    
    include("smoketests.jl")
    include("errors_and_warnings.jl")

    println("DTM: Testing finished.")
end

@show telaTotal = (time_ns()-t0)/1e9