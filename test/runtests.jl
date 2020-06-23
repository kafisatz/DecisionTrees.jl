t0 = time_ns()

# test
import CSV
import DataFrames
import DataFrames: DataFrame
import Random
import OnlineStats
import StatsBase

using Distributed # only for Tests?
using Statistics
using Pkg # only for Tests?
using Test # only for Tests? 
# using Pkg; Pkg.activate(normpath(joinpath(pwd(),".."));Pkg.instantiate();#Pkg.test() 
using DecisionTrees

println("Time for using & import: ",round((time_ns() - t0) / 1e9, digits=3),"s")

# set directories
global pkgdir = ""
global testdir = joinpath(dirname(@__FILE__), "") # joinpath(@__DIR__,"..","test",filename)
# the following if clause is for "direct execution" (e.g. in VScode where @__FILE__ does not seem to work)
if testdir == ""
    try 
        global pkgdir = normpath(joinpath(dirname(pathof(DecisionTrees)), ".."))
        global testdir = joinpath(pkgdir, "test")
    catch
    end
else 
    global pkgdir = normpath(joinpath(testdir, ".."))
end
datadir = joinpath(pkgdir, "data")

Random.seed!(1239946)
DecisionTreesTests = @testset "DecisionTrees" begin
    t0runtests = time_ns()
    println("DTM: Running tests...") # Time = $(t0runtests)")
        
    global eltypesData1 = [Int,	Float64,	Float64,	Float64,	Float64,	Float64,	Int,	String,	String,	String,	Int,	String,	String,	String,	String,	String,	String,	String,	String,	String,	String,	Int,	String,	Int,	String,	String,	Int,	String,	String,	Int,	String,	String,	Int,	String,	String,	Int,	String,	Int,	String,	String,	Int,	Int,	String,	Int,	Int,	String,	Int,	Int,	String,	Int,	String,	String,	String,	Int,	String,	String,	String,	String,	Int,	Int,	Int,	String,	String,	String,	String,	String,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Float64,	Int,	Float64,	Float64,	Float64,	Float64,	Int,	Int,	Int,	Int,	Int,	Int]
    
    @info("DTM: Note that the tests produce a relatively long log and write a number of files in temporary directories. You can disregard the log, if all tests pass.")
    @warn("DTM: So far only very rudimentary tests have been implemented!")
    
    include("smoketests.jl")
    include("errors_and_warnings.jl")
    include("functionalityTests.jl")
    include("helperfunctions.jl")

    @testset "Other Tests" begin
        @test 1 == 1
    end
    println("DTM: Testing finished. Time needed:")
    TimeForTesting = (time_ns() - t0runtests) / 1e9
    println((round(TimeForTesting, digits=0)), "s")
    println(round(TimeForTesting / 60, digits=1), "m")
    println(round(TimeForTesting / 3600, digits=3), "h")
end
# dump(DecisionTreesTests)