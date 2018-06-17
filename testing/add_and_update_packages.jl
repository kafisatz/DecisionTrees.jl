#this is file is only a helper file to add and update packages
add PyCall CSV BenchmarkTools OnlineStats StatsBase DataFrames ProgressMeter JLD2 FileIO StatsFuns

up PyCall CSV BenchmarkTools OnlineStats StatsBase DataFrames ProgressMeter JLD2 FileIO StatsFuns

#consider (in the Pkg manager)
add OnlineStats#master

using PyCall #works on master
using CSV #works on masters
using OnlineStats #not working on June 7
using StatsBase
using ProgressMeter
using JLD2#does it work on master?
using FileIO
using StatsFuns #works on master
