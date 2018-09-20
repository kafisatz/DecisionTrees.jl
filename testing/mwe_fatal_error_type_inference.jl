#sample script to reproduce
#ERROR: fatal error in type inference (type bound)

#it seems related to Union{mseSplit,sseSplit}
#(search the code for that expression)
#I am quite sure that this did not happen on 0.7beta2

using Distributed
using Pkg
pkg"develop https://github.com/kafisatz/DecisionTrees.jl#master"
pkg"activate DecisionTrees.jl"
Distributed.@everywhere using DecisionTrees

#if the above does not work, maybe you need to clone the git
#git clone --depth=50 --branch=master https://github.com/kafisatz/DecisionTrees.jl.git 
#cd kafisatz/DecisionTrees.jl

@assert VERSION>=v"0.7.0-rc2.0"

Distributed.@everywhere import CSV
Distributed.@everywhere import DataFrames
Distributed.@everywhere import DataFrames: DataFrame

##############################
#Read the data
##############################

dr="C:\\Users\\bernhard.konig\\Documents\\ASync\\home\\Code\\Julia\\DecisionTrees.jl"
isdir(dr)&&cd(dr)
pwd()

@warn("""This assumes pwd() == Pkg.dir("DecisionTrees.jl")""") 
datafile=joinpath("data","freMTPL2","freMTPL2.csv")
if Sys.islinux()
    #trying to unzip the file...
    pathToZip=string(splitext(datafile)[1],".zip") 
    unzipLoc=splitdir(pathToZip)[1]
    cmd=`unzip $(pathToZip) -d $(unzipLoc)`
    run(cmd)
end

@assert isfile(datafile) "You need to unzip the file freMTPL2.zip in the folder data/freMTPL2"

@time fullData=CSV.read(datafile,rows_for_type_detect=100000,allowmissing=:none,categorical=false);

#some data manipulation
    areasSorted=sort(unique(fullData[:Area]))
    AreaInteger=map(x->findall((in)([x]),areasSorted)[1],fullData[:Area])
    fullData[:AreaInteger]=AreaInteger
    #correct for unreasonable observations
    for i=1:size(fullData,1)
        fullData[:ClaimNb][i]=min(4,fullData[:ClaimNb][i])
        fullData[:Exposure][i]=min(1.0,fullData[:Exposure][i])
    end
    
#set independent variables
selected_explanatory_vars=["Area","AreaInteger","VehPower","VehAge","DrivAge","BonusMalus","VehBrand","VehGas","Density","Region"]

##############################
#Prepare the data
##############################
dtmtable,sett,dfprepped=prepare_dataframe_for_dtm!(fullData,keycol="IDpol",trnvalcol="trnTest",numcol="ClaimNb",denomcol="Exposure",weightcol="Exposure",independent_vars=selected_explanatory_vars);

#run a few models

#this runs as expected
    updateSettings!(sett,minWeight=-0.03,model_type="build_tree",calculatePoissonError=true)
    resultingFiles,resM=dtm(dtmtable,sett)

#this runs as expected
    updateSettings!(sett,crit="sse",model_type="build_tree")
    resultingFiles,resM=dtm(dtmtable,sett)

#this one fails!
    updateSettings!(sett,crit="sse",model_type="boosted_tree")
    resultingFiles,resM=dtm(dtmtable,sett)

#=
# this is the error I get

"""
ERROR: fatal error in type inference (type bound)
Stacktrace:
 [1] _typed_hcat(::Type{Float64}, ::Tuple{Array{Int64,1},Array{Int64,1},Array{Float64,1},Array{Float64,1},Array{Float64,1},Array{Float64,1},Array{Float64,1},Array{Float64,1}}) at .\abstractarray.jl:1255
 [2] typed_hcat at .\abstractarray.jl:1243 [inlined]
 [3] hcat at C:\cygwin\home\Administrator\buildbot\worker\package_win64\build\usr\share\julia\stdlib\v0.7\SparseArrays\src\sparsevector.jl:1063 [inlined]
 [4] macro expansion at .\logging.jl:307 [inlined]
 [5] macro expansion at .\util.jl:156 [inlined]
 [6] macro expansion at C:\Users\bernhard.konig\Documents\ASync\home\Code\Julia\DecisionTrees.jl\src\boosting.jl:89 [inlined]
 [7] macro expansion at .\util.jl:156 [inlined]
 [8] macro expansion at C:\Users\bernhard.konig\Documents\ASync\home\Code\Julia\DecisionTrees.jl\src\helper_functions.jl:1469 [inlined]
 [9] boosted_tree(::DecisionTrees.DTMTable, ::ModelSettings) at C:\Users\bernhard.konig\Documents\ASync\home\Code\Julia\DecisionTrees.jl\src\boosting.jl:82
 [10] run_model_actual(::DecisionTrees.DTMTable, ::ModelSettings, ::String) at C:\Users\bernhard.konig\Documents\ASync\home\Code\Julia\DecisionTrees.jl\src\runmodel.jl:467
 [11] #dtm#110 at C:\Users\bernhard.konig\Documents\ASync\home\Code\Julia\DecisionTrees.jl\src\runmodel.jl:798 [inlined]
 [12] dtm(::DecisionTrees.DTMTable, ::ModelSettings) at C:\Users\bernhard.konig\Documents\ASync\home\Code\Julia\DecisionTrees.jl\src\runmodel.jl:798
 [13] top-level scope at none:0
 """
 =#