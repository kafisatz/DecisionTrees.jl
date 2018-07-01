__precompile__()
true

#DTM is used as an abbreviation for DecisionTreeModels
#note on v0.7 alpha precompilation through 'using' (i.e. when the code was not precomipled before) takes about 800 seconds
#once precompilation is done, the using statement might need around 500 seconds. 
@warn "BK: DTM may have a performance regression of roughly 50% in Julia 0.7alpha versus 0.6.3. \nYou may want to start julia with --depwarn=no option (in case there are still deprecation warnings in this package or its dependencies). I expect that this will vanish as 0.7 matures (or the bottleneck will hopefully be identified through profiling)"

#=
@info "DTM: BK Possibly add trnidx and validx to the resulting ensemble. This is relevant in case of a CV sampling which is performed. Otherwise it is not possible to reconstruct the Excel statistics after the model has run."
@info "DTM: BK 'Remove' variables which are not used by a model from the SAS/VBA/CSharp code (e.g. dim command in VBA)"
@info "DTM: BK Add 'time_finished' to model result (and possibly the time needed for the modelling)"
@info "DTM: BK Need to ensure that the code runs smoothly even when no split is found (e.g. minw too big)"
@info "DTM: BK Possibly introduce an option not to calculate certain intermediate metrics (scoring, sorting, ...) -> faster runtime"
@info "DTM: BK To check if sortperm! performance has improved in newer Julia versions"
@info "DTM: BK To fix this: the code should not rely on dtmtable.mappings and dtmtable.candmatwomax anymore (if possible). The pools have the equivalent information"
@info "DTM: BK To add purity improvement for each split! Stop splitting if it is > some_specified_threshold"
@info "DTM: BK Do not use df_name_vector anymore, but the symbol in the subset or simply names(features)"
@info "DTM: BK Add yes/no to the top of each printed DOT graph (in the PDF output) indicating that 'yes goes left'"
@info "DTM: BK Try to avoid branching in critical functions; use ifelse or similar constructs (this might improve the performance)"
=#

module DecisionTrees

export run_model,run_model_actual,define_eltypevector,prepare_dataframe_for_dtm!,updateSettingsMod!,prep_data_from_df_modernized,run_legacy,prep_data_from_df

#include PooledArraysDTM module
include("PooledArraysDTM.jl")
import .PooledArraysDTM: PooledArray

import Core 
import Dates
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

import DataStreams #We do not really need DataStreams (explicitly) but CSV keeps failing (in tests) if this is not here (on 0.7alpha)

#import MySQL, # temporarily disabled 
#import SQLite #disabled as it uses DataFrames 0.11

Random.srand(1234)

#include("check_if_python_packages_are_installed.jl")
include("types.jl")
include("cross_validation_samplers.jl") #needs to be included early on as types are in here
include("sorting.jl")
include(joinpath("splittingmeasures","splittingmeasures.jl"))
include("helper_functions.jl")
include("apply_tree_fn.jl")
include("build_tree.jl")
include("boosting.jl")
include("bagging.jl")
include("runmodel.jl")
include("runmodel_multi.jl")
include("charts.jl")
include("dot_graphs.jl")
include("roc.jl")
include("cross_validation.jl")
include("show.jl")
include("partialDependence.jl")
include("apply_tree_fn_unseen.jl")

global const pyModPandas = PyCall.PyNULL()
global const pyModxlsxwriter = PyCall.PyNULL()

global const emptyMeanVarSeries=OnlineStats.Series(OnlineStats.Mean(),OnlineStats.Variance())
global const globalConstAllowableMethodsForDefineCandidates=("equalWeight","basedOnWeightVector")
global const UInt8VECTORemptySplitDef = Vector{Splitdef{UInt8}}(undef,0)
global const UInt16VECTORemptySplitDef = Vector{Splitdef{UInt16}}(undef,0)
global const UInt8emptySplitDef=Splitdef(0,0,Symbol(),Vector{UInt8}(),-Inf,0.0,0.0) #formerly const_default_splitdef
global const UInt16emptySplitDef=Splitdef(0,0,Symbol(),Vector{UInt16}(),-Inf,0.0,0.0) 
global const defaultModelName="dtmresult"
global const defaultModelNameWtihCSVext=string(defaultModelName,".csv")
global nLevelsThreshold=2002
global const global_number_of_num_f_warning_mandatory_field="Abort. Settings do not contain any value for number_of_num_features (which is mandatory)."
global const stars="*******************************************************************************"
global const UInt8emptyRulepath=Array{Rulepath{UInt8}}(undef,0)
global const UInt16emptyRulepath=Array{Rulepath{UInt16}}(undef,0)
global const UInt8emptyLeaf=Leaf(0,NaN,NaN,NaN,-1,UInt8emptyRulepath,NaN,NaN,0)
global const UInt16emptyLeaf=Leaf(0,NaN,NaN,NaN,-1,UInt16emptyRulepath,NaN,NaN,0)
global const UInt8emptyNode=Node(0,0,Array{UInt8}(undef,0),deepcopy(UInt8emptyLeaf),deepcopy(UInt8emptyLeaf),deepcopy(UInt8emptyRulepath))
global const UInt16emptyNode=Node(0,0,Array{UInt16}(undef,0),deepcopy(UInt16emptyLeaf),deepcopy(UInt16emptyLeaf),deepcopy(UInt16emptyRulepath))
global const global_const_shift_cols=5
global const global_pldamod_valid_types=[String,Float64]	

global const global_nameOfSettingsSheet="ModelSettings"
global const global_nameOfModelStatisticsSheet="ModelStatistics"
global const global_nameOfTWOWayValidationSheet="2Way_ValidationGraphs"

global const DifferenceSplitConst      =  DifferenceSplit() #todo/tbd check if we need the constants like Difference at all (I think we do not need these)
global const MaxValue				   =  MaxValueSplit()
global const MaxAbsValue			   =  MaxAbsValueSplit()
global const MaxMinusValue			   =  MaxMinusValueSplit()
global const MSE                       =  MSESplit()
global const RankOpt                   =  RankOptSplit()
global const ROptMinRLostPct           =  ROptMinRLostPctSplit()
global const ROptMinRLost              =  ROptMinRLostSplit()
global const SORTBYMEAN				   =  SortByMean()

global const globalValidfitForStatsAndCharts=["rawRelativities","unsmoothedPerScore","smoothedPerScore"]	
global const global_statsperiter_header=["Iteration" "Correlation" "Std Dev Trn" "Std Dev Val" "Std Dev Ratio" "Lift Trn" "Lift Val" "Reversals Trn" "Reversals Val" "RelDiff TrnObs vs ValObs" "RelDiff ValObs vs TrnFitted" "Min Rel Trn" "Max Rel Trn" "Min Observed Ratio Trn" "Max Observed Ratio Trn"  "Min Rel Val" "Max Rel Val" "Min Observed Ratio Val" "Max Observed Ratio Val" "Gini Num (Two Arg Fn) Unweighted Trn" "Gini Num (Two Arg Fn) Unweighted Val" "Gini Num Trn (Single Agrument Fn)" "Gini Num Val (Single Argument Fn)" "RSquared of Ratio Trn" "RSquared of Ratio Val" "RSquared of Numerator Trn" "RSquared of Numerator Val" "RSS of Ratio Trn" "RSS of Ratio Val" "RSS of Numerator Trn" "RSS of Numerator Val" "Average Poisson Error Trn" "Average Poisson Error Val"]
global const GLOBALperfMeasursesWhereMaxIsBest=["Correlation","Lift Trn","Lift Val","Max Rel Trn","Max Observed Ratio Trn","Max Rel Val","Max Observed Ratio Val","Gini Num (Two Arg Fn) Unweighted Trn","Gini Num (Two Arg Fn) Unweighted Val","Gini Num Trn (Single Agrument Fn)","Gini Num Val (Single Argument Fn)","RSquared of Ratio Trn","RSquared of Ratio Val","RSquared of Numerator Trn","RSquared of Numerator Val","RSS of Ratio Trn","RSS of Ratio Val","RSS of Numerator Trn","RSS of Numerator Val"]

global const global_debug_csharp=false
global const global_byte_order_mark='\ufeff'
global const CSHARP_VALID_CHARS="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_123456789äëöüïÿâêîôûàèìòùáéíóúý" #the list is probably much longer, but this will do for now

function __init__()
	#the following lines may trigger the installation of the respective python packages
	#if true #temporarily disable this
        copy!(pyModPandas, PyCall.pyimport_conda("pandas","pandas"))
        copy!(pyModxlsxwriter, PyCall.pyimport_conda("xlsxwriter","xlsxwriter"))
    #end
end

function get_sha1()
	#which_res=@which DecisionTrees.__init__()
	#dtpath=splitdir(splitext(string(which_res.file))[1])[1][1:end-4]
	dtpath=dirname(@__FILE__)
	current_path=pwd()
	cd(dtpath)	
	sha1_of_dt="n/a"
	try
		cmd=`git rev-parse HEAD`
		sha1_of_dt=read(cmd,String)[1:end-1] #need to discard \n at the end
	catch 
		sha1_of_dt="unable to determine version. Git invocation possibly failed"
	end 
	cd(current_path)
	return sha1_of_dt
end

#Precompile files
    include(joinpath("precompile","precompile_DecisionTrees.jl"))
	#include("precompile.jl")

global const DoubleQuote='\"' #note, this is here at the end only because Notepad++ messes up the formatting after this definition
	
end #end Module DTM
