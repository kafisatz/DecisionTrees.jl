t0=time_ns()

using CSV
using DataFrames
using Test 
using DecisionTrees 

@show tela = (time_ns()-t0)/1e9

#set directories
global pkgdir=""
global testdir = joinpath(dirname(@__FILE__),"")
#the following if clause is for "direct execution" (e.g. in VScode where @__FILE__ does not seem to work)
if testdir==""
    try 
        pkgdir=Pkg.dir("DecisionTrees") #how long will this remain supported?
        testdir=joinpath(pkgdir,"test")
    catch
    end
else 
    pkgdir=normpath(joinpath(testdir,".."))
end
datadir=joinpath(pkgdir,"data")

@testset "DecisionTrees" begin
    println("DTM: Running tests:")

    warn("So far only very rudimentary tests have been implemented!")
    @test 1==1
    include("smoketests.jl")

    println("DTM: Testing finished.")
end

@show telaTotal = (time_ns()-t0)/1e9