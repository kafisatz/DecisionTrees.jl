cd(Pkg.dir("DecisionTrees"))
import Core 
import Dates
println(Dates.now())
import Random
import Distributed
import DelimitedFiles
import PyCall
import DataFrames
import DataFrames: DataFrame
import JLD2
import OnlineStats
import ProgressMeter
import StatsBase
import DataStreams
import CSV

#try PyCall
global const pyAB = PyCall.PyNULL()
global const pyCD = PyCall.PyNULL()
copy!(pyAB, PyCall.pyimport_conda("pandas","pandas"))
copy!(pyCD, PyCall.pyimport_conda("xlsxwriter","xlsxwriter"))

#using dt:
using DecisionTrees
println(Dates.now())