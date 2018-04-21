__precompile__()
VERSION >= v"0.7-"
@info "DTM: BK to add trnidx and validx to the resulting ensemble. This is relevant in case of a CV sampling which is performed. Otherwise it is not possible to reconstruct the Excel statistics after the model has run."
@info "DTM: BK to 'remove' variables which are not used by a model from the SAS/VBA/CSharp code (e.g. dim command in VBA)"
@info "DTM: BK add 'time_finished' to model result (and possibly the time needed for the modelling)"
@info "DTM: BK need to ensure that the code runs smoothly even when no split is found (e.g. minw too big)"
@info "DTM: BK possibly introduce an option not to calculate intermediate metrics/scoring/sorting -> faster runtime"
@info "DTM: BK to check if sortperm! performance has improved in newer Julia versions"
@info "DTM: BK to fix this - we should not rely on mappings and candmatwomax anymore (if possible)"
@info "DTM: BK to add purity improvement for each split! Stop splitting if it is > threshold"
@info "BK: do not use df_name_vector anymore, but the symbol in the subset or simply names(features) - possibly also get rid of mappings and candmat...."
@info "BK: tbd add yes/no to the top of each printed DOT graph indicating that 'yes goes left'"
@info "BK: avoid branching in critical functions / use ifelse or similar constructs"

module DecisionTrees

export run_model,run_model_actual,define_eltypevector,prepare_dataframe_for_dtm!,updateSettingsMod!,prep_data_from_df_modernized,run_legacy,prep_data_from_df

using Random,OnlineStats, StatsBase, StatsFuns, DataFrames, ProgressMeter, PyCall, JLD2, FileIO,PooledArrays # DataArrays is not used anymore
#as of March 4, JLD2 is not working due to this https://github.com/bicycle1885/CodecZlib.jl/issues/17
@warn "FileIO and PooledArrays are DIRTY for now (as they were not yet updated on master)"

#using MySQL, # temporarily disabled 
#using SQLite #disabled as it uses DataFrames 0.11

srand(1234)

#include("check_if_python_packages_are_installed.jl")
include("types.jl")
include("cross_validation_samplers.jl") #needs to be included early on as types are in here
include("sorting.jl")
include("splittingmeasures\\splittingmeasures.jl")
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
#include("onlineStats_addons.jl")
#include("sqlite_query.jl")



global const pyModPandas = PyNULL()
global const_default_splitdef=Splitdef(0,0,Symbol(),Vector{UInt8}(),-Inf,0.0,0.0)

function __init__()
	#we initialize Python here in the beginning, such that Julia would crash (prior to the modelling) in case there is an installation error with python
	#InitPython.initializePython() #todo/tbd check if and when this line is needed
	global const juliaStartedat=now() 
	#unshift!(PyVector(pyimport("sys")["path"]), "") #By default, PyCall doesn't include the current directory in the Python search path. If you want to do that (in order to load a Python module from the current directory) we need to run this
	#unshift!(PyVector(pyimport("sys")["path"]),"C:\\Anaconda3\\Lib\\site-packages")
	#global const pyModPandas=pywrap(pyimport("pandas")) #this line should be equivalent to  @pyimport pandas as pyModPandas
    copy!(pyModPandas, pyimport_conda("pandas","pandas"))

	global const global_number_of_num_f_warning_mandatory_field="Abort. Settings do not contain any value for number_of_num_features (which is mandatory)."
	global const stars="*******************************************************************************"
	global const emptyRulepath=Array{Rulepath}(0)
	global const emptyLeaf=Leaf(0,NaN,NaN,NaN,-1,emptyRulepath,NaN,NaN,0)
	global const emptyNode=Node(0,0,Array{UInt8}(0),deepcopy(emptyLeaf),deepcopy(emptyLeaf),deepcopy(emptyRulepath))
	global const global_const_shift_cols=5+1-1
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
	#const empty_expr=:("undefined"=="undef")
	#const STRING_WHICH_DOES_NOT_OCCUR_IN_DATA=convert(String,"Z.\???!?Z")
	#global const splittingcrits=["_mse" "_difference" "_difference_rankopt" "_minimize_lost_ranks_maximize_5pct" "_minimize_lost_ranks_maximize_minresidual" "_RankOpt" "_minimize_lost_ranks_difference" "_difference_vectorized" "_difference2" "_difference_devectorized_simd" "_mean_abs_error" "_mean_abs_error_rel" "_sum_abs_error" "_reduction_in_variance" "_reduction_in_standarddeviation" "_min_std" "_min_mean" "_max_mean"]

	global const global_statsperiter_header=["Iteration" "Correlation" "Std Dev Trn" "Std Dev Val" "Std Dev Ratio" "Lift Trn" "Lift Val" "Reversals Trn" "Reversals Val" "RelDiff TrnObs vs ValObs" "RelDiff ValObs vs TrnFitted" "Min Rel Trn" "Max Rel Trn" "Min Observed Ratio Trn" "Max Observed Ratio Trn"  "Min Rel Val" "Max Rel Val" "Min Observed Ratio Val" "Max Observed Ratio Val" "Gini Num (Two Arg Fn) Unweighted Trn" "Gini Num (Two Arg Fn) Unweighted Val" "Gini Num Trn (Single Agrument Fn)" "Gini Num Val (Single Argument Fn)" "RSquared of Ratio Trn" "RSquared of Ratio Val" "RSquared of Numerator Trn" "RSquared of Numerator Val" "RSS of Ratio Trn" "RSS of Ratio Val" "RSS of Numerator Trn" "RSS of Numerator Val"]

	global const global_debug_csharp=false
	global const global_byte_order_mark='\ufeff'
	global const CSHARP_VALID_CHARS="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_123456789äëöüïÿâêîôûàèìòùáéíóúý" #the list is probably much longer, but this will do for now
	global const DoubleQuote='\"' #note, this is here at the end only because Notepad++ messes up the formatting after this definition
end

function get_sha1()
	which_res=@which DecisionTrees.__init__()
	dtpath=splitdir(splitext(string(which_res.file))[1])[1][1:end-4]
	current_path=pwd()
	cd(dtpath)
	cmd=`git rev-parse HEAD`
	#txa=run(cmd)
	sha1_of_dt=readstring(cmd)[1:end-1] #need to discard \n at the end
	cd(current_path)
	return sha1_of_dt
end

end #end Module DTM
