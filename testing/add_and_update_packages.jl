#this is file is for testing purposes only.
add CSV Loess BenchmarkTools Compat OnlineStats StatsBase Distributions DataFrames ProgressMeter Iterators JLD2 ProfileView Coverage FileIO Random StatsFuns

up CSV Loess BenchmarkTools Compat OnlineStats StatsBase Distributions DataFrames ProgressMeter Iterators JLD2 ProfileView Coverage FileIO Random StatsFuns

#consider (in the Pkg manager)
add OnlineStats#master

using PyCall #works on master
using CSV #works on masters
using Compat
using OnlineStats #not working on June 7
using StatsBase
using Distributions #needed to Pkg.build("Rmath")  and Pkg.build("SpecialFunctions") by hand using DataFrames
using ProgressMeter
using Iterators
using JLD2#does it work on master?
using FileIO
using StatsFuns #works on master
