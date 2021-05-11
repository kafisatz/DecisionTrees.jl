t0 = time_ns()
#=
    cd("/home/bernhard/DecisionTrees.jl");using Pkg;Pkg.activate(pwd());
    Pkg.test();
=#

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
DecisionTreesTests = @testset verbose = true "DecisionTrees" begin
    t0runtests = time_ns()
    println("DTM: Running tests...") # Time = $(t0runtests)")
        
    global eltypesData1 = [Int,	Float64,	Float64,	Float64,	Float64,	Float64,	Int,	String,	String,	String,	Int,	String,	String,	String,	String,	String,	String,	String,	String,	String,	String,	Int,	String,	Int,	String,	String,	Int,	String,	String,	Int,	String,	String,	Int,	String,	String,	Int,	String,	Int,	String,	String,	Int,	Int,	String,	Int,	Int,	String,	Int,	Int,	String,	Int,	String,	String,	String,	Int,	String,	String,	String,	String,	Int,	Int,	Int,	String,	String,	String,	String,	String,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Float64,	Int,	Float64,	Float64,	Float64,	Float64,	Int,	Int,	Int,	Int,	Int,	Int]
    
    @info("DTM: Note that the tests produce a relatively long log and write a number of files in temporary directories. You can disregard the log, if all tests pass.")
    @warn("DTM: So far only very rudimentary tests have been implemented!")
    
    global graphvizexe
    global graphvizsetup_is_working
    graphvizsetup_is_working = false
    @testset verbose = true "GraphViz Setup" begin 
        if Sys.iswindows()
            #graphvizexe = raw"c:\program files\graphviz\bin\dot.exe" #as per github ci log of graphviz setup this seems to be the location of dot.exe
            graphvizexe = "dot"
            #if !isfile(graphvizexe)
                if haskey(ENV,"graphvizdot")
                    #if key exists this takes precedence
                    graphvizexe = ENV["graphvizdot"]
                    #=
                        ENV["graphvizdot"]="c:\\program files\\Graphviz\\bin\\dot.exe"
                    =#
                end 
            #end
        else 
            graphvizexe = "dot" #linux & appleOS assume it is in in path
        end
        
        #@test isfile(graphvizexe)
        gcmd = `$(graphvizexe) "-V"`
        try 
            gres = run(gcmd)
            @test gres.exitcode == 0
            @show gres.exitcode
            graphvizsetup_is_working = true
            @test true #graphviz setup works
        catch e
            @warn("DTM: unable to invoke gcmd")
            @show gcmd
            print(e)
            @test false
        end
    end

    include("smoketests.jl") #running a number of functions to see if anything 'smokes' (runs into an unexpected error)
    include("errors_and_warnings.jl") #testing a few warnings and errors
    include("functionalityTests.jl") #testing resulting values of models
    include("helperfunctions.jl")

    @testset verbose = true "Other Tests" begin
        @test 1 == 1
    end
    println("DTM: Testing finished. Time needed:")
    TimeForTesting = (time_ns() - t0runtests) / 1e9
    println((round(TimeForTesting, digits=0)), "s")
    println(round(TimeForTesting / 60, digits=1), "m")
    println(round(TimeForTesting / 3600, digits=3), "h")
end
# dump(DecisionTreesTests)