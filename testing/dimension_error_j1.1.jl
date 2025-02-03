abstract type SplittingCriterion  end
struct DifferenceSplit      <: SplittingCriterion end

abstract type SortBy end
struct SortByMean 	<: SortBy end

mutable struct ModelSettings
	model_type::String # 1
	minWeight::Float64 # 2
	randomw::Float64 # 3
	crit::SplittingCriterion # fn version of 4
	maxSplittingPoints::Int # 5
	iterations::Int # 6
	learningRate::Float64 # 7 moderationfactor
	nScores::Int # 8
	adaptiveLearningRate::Float64 # 9
	prem_buffer::Int # 11
	startAtMean::Bool # 12
	writeTree::Bool # 13	
	number_of_num_features::Int # 14
	spawnsmaller::Bool # 26
    boolRankOptimization::Bool # 35  
	boolRandomizeOnlySplitAtTopNode::Bool # 37
	subsampling_prop::Float64 # 38 #NOTE Friedman suggests to sample without replacement see Friedman 2002 Stochastic Gradient Boosting,” Computational Statistics and Data Analysis 38(4):367-378
	subsampling_features_prop::Float64 # 39
	version::String # 40
	preppedJLDFileExists::Bool # 41
	catSortByThreshold::Int # 42
	catSortBy::SortBy # 43
	scorebandsstartingpoints::Array{Int,1} # 44
	showTimeUsedByEachIteration::Bool # 45
	smoothEstimates::String # 46
	deriveFitPerScoreFromObservedRatios::Bool # deriveFitPerScoreFromObservedRatios::Bool (if this is true then the fit per score will be based on the observed ratio per score (instead of the fitted ratio per score)) . Note that this does not influence the structure of the tree, but only the fitted ratio per score (and scoreband)
    roptForcedPremIncr::Bool # 47
    writeSasCode::Bool # 49
	writeIterationMatrix::Bool # 50
	writeResult::Bool # 51
	writeStatistics::Bool # 52
	boolCreateZipFile::Bool # 53
	writeCsharpCode::Bool
	writeVbaCode::Bool
	nDepthToStartParallelization::Int
	baggingWeightTreesError::String
	cBB_niterBoosting::Int
	cBB_niterBagging::Int
	fixedinds::Array{Int,1}
	boolNumeratorStats::Bool	
	statsByVariables::Array{Int,1}
	statsRandomByVariable::Int
	saveJLDFile::Bool # if this is set to false no *.jld2 file is saved #todo we should rename this to ensure it is clear that this refers to the prepped data only!
	saveResultAsJLDFile::Bool # this refers to the model output/result
	print_details::Bool # whether julia should print details to the log or not (note to developers: warnings/errors and certain important infos should always be printed!)
	seed::Int # optional entropy
	graphvizexecutable::String # e.g. C:\Program Files (x86)\Graphviz2.38\ if this
	showProgressBar_time::Bool # bar form iterators package (or some other package)
	prroduceEstAndLeafMatrices::Bool
	write_dot_graph::Bool
	calculateGini::Bool # whether or not to calculate the gini (which requires sorting of the data (which takes time!))
	calculatePoissonError::Bool
    performanceMeasure::String
	fitForStatsAndCharts::String
	ignoreZeroDenominatorValues::Bool

	# the following are treated specially
	df_name_vector::Array{String,1}
    number_of_char_features::Int # this is calculated by Julia from number_of_num_features and ncolsdfIndata (which is defined by the data)

	chosen_apply_tree_fn::String # defined below
	moderationvector::Array{Float64,1} # 50 onwards

    function ModelSettings()
	# initialize default settings
	model_type = "build_tree" # 1
	minWeight = -0.1 # 2
	randomw = 0.0 # 3
	crit = DifferenceSplit()	 # 4
	maxSplittingPoints = 250 # 5 #arbitrary choice. Preliminary testing has shown no performance difference for different values (say 10 or 200) for this setting
	iterations = 2 # 6
	learningRate = 0.1 # 7
	nScores = 1000 # 8
	adaptiveLearningRate = 1.0 # 9
	prem_buffer = 0 # 11
	startAtMean = true # 12
	writeTree = true # 13
	number_of_num_features = -1 # 14    
	spawnsmaller = true # 26
	boolRankOptimization = false # 35
	boolRandomizeOnlySplitAtTopNode = true # 37
	subsampling_prop = 1.0 # 38
	subsampling_features_prop = 1.0 # 39
	version = "not_initialized" # 40
	preppedJLDFileExists = false # 41
	catSortByThreshold = 7 # 42
	catSortBy = SortByMean() # 43
	scorebandsstartingpoints = [parse(Int, x) for x in split("1,100,200,300,400,500,600,700,800,900", ',')] # scorebandsstartingpoints=parse(Int,split("1,100,200,300,400,500,600,700,800,900",',')) #44
	showTimeUsedByEachIteration = false # true #45
	smoothEstimates = "1" # 46 this is a string for now (as there could be different smoothing methods specified by this string)
	deriveFitPerScoreFromObservedRatios = true
	roptForcedPremIncr = false # 47
        writeSasCode = false # 49
	writeIterationMatrix = false# 50
	writeResult = false # true #51
	writeStatistics = true # 52
	boolCreateZipFile = false # true #53
	writeCsharpCode = false # true
	writeVbaCode = false
	nDepthToStartParallelization = -1
	baggingWeightTreesError = "mae" # ::String
	# combined bagging&boosting model
	cBB_niterBoosting = 0
	cBB_niterBagging = 0
	fixedinds = Array{Int}(undef, 0)
	boolNumeratorStats = false	
	statsByVariables = Int[]
	statsRandomByVariable = 5
	saveJLDFile = true
	saveResultAsJLDFile = false
	print_details = true
	seed = 90210
	graphvizexecutable = "" # C:\\Program Files\\Graphviz\\bin\\dot.exe""
	showProgressBar_time = true
	prroduceEstAndLeafMatrices = false
	write_dot_graph = false
	calculateGini = false
	calculatePoissonError = false
        performanceMeasure = "Lift Val"
	fitForStatsAndCharts = "rawRelativities" # rawRelativities,unsmoothedPerScore,smoothedPerScore
	ignoreZeroDenominatorValues = false

	# the following are treated specially
	df_name_vector = Array{String}(undef, 0) #
        number_of_char_features = -1 #
	chosen_apply_tree_fn = "apply_tree_by_leaf" # #apply_tree_by_row does not seem to work for (certain?) boosting models
	moderationvector = [0.1] #

	return new(model_type,minWeight,randomw,crit,maxSplittingPoints,iterations,learningRate,nScores,adaptiveLearningRate,prem_buffer,startAtMean,writeTree,number_of_num_features,spawnsmaller,boolRankOptimization,boolRandomizeOnlySplitAtTopNode,subsampling_prop,subsampling_features_prop,version,preppedJLDFileExists,catSortByThreshold,catSortBy,scorebandsstartingpoints,showTimeUsedByEachIteration,smoothEstimates,deriveFitPerScoreFromObservedRatios,roptForcedPremIncr,writeSasCode,writeIterationMatrix,writeResult,writeStatistics,boolCreateZipFile,writeCsharpCode,writeVbaCode,nDepthToStartParallelization,baggingWeightTreesError,cBB_niterBoosting,cBB_niterBagging,fixedinds,boolNumeratorStats,statsByVariables,statsRandomByVariable,saveJLDFile,saveResultAsJLDFile,print_details,seed,graphvizexecutable,showProgressBar_time,prroduceEstAndLeafMatrices,write_dot_graph,calculateGini,calculatePoissonError,performanceMeasure,fitForStatsAndCharts,ignoreZeroDenominatorValues
	 ,df_name_vector,number_of_char_features,chosen_apply_tree_fn,moderationvector)
    end  # ModelSettings()

end # type definition


function writeAllFieldsToArray(s)	
	nn = fieldnames(typeof(s))
	res = Array{AbstractString}(undef, length(nn), 2)
	i = 0
	for x in nn
		i += 1
        @show x
        tmp = AbstractString[string(x) string(getfield(s, x))][:]
        @show tmp 
        @show size(tmp)
        @show typeof(tmp)
		res[i,:] = AbstractString[string(x) string(getfield(s, x))][:]
	end
	return res
end

mi = ModelSettings()
foo = writeAllFieldsToArray(mi)
