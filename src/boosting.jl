# Boosting
function boosted_tree(dtmtable::DTMTable, sett::ModelSettings)		
	trnidx = dtmtable.trnidx
	validx = dtmtable.validx
	features = dtmtable.features
	weight = dtmtable.weight
	actualNumerator = dtmtable.numerator
	denominator = dtmtable.denominator
	mappings = dtmtable.mappings
	candMatWOMaxValues = dtmtable.candMatWOMaxValues
	obstrn, obsval, trn_meanobservedvalue, val_meanobservedvalue, trn_numtot, val_numtot, trn_denomtot, val_denomtot, empty_rows_after_iteration_stats, showTimeUsedByEachIteration, chosen_apply_tree_fn, startAtMean, adaptiveLearningRate, moderationvector, iterations, nameOflistofpredictorsSheet, nameOfpredictorsSheet, nameOfModelStatisticsSheet, nameOfScoresSheet, nameOfOverviewSheet, nameOfSettingsSheet, nameOfValidationSheet, xlData, showProgressBar_time, prroduceEstAndLeafMatrices = initSettingsWhichAreTheSameForBoostingAndBagging(trnidx, validx, actualNumerator, denominator, sett)
	# trn_meanobservedvalue is the mean observed RATIO, i.e. sum(numerator)/sum(denominator) for the training data
	T_Uint8_or_UInt16 = find_max_type(features)
    # obs=obstrn+obsval
    obs = size(features, 1)
    # obs=obstrn+obsval  assumes that trn and val idx form a partition!
    if size(features, 1) != obstrn + obsval 
#        @warn("DTM: obstrn+obsval != size(features,1). This is experimental.")
    end
    
	current_error = 0.0
    moderationfactor = moderationvector[1]
	if size(moderationvector, 1) == 1
		actual_moderationvector = Array{Float64}(undef, iterations)
		fill!(actual_moderationvector, copy(moderationvector[1]))
	else
		actual_moderationvector = copy(moderationvector)
	end
	if startAtMean
		estimatedRatio = ones(length(weight)) .* trn_meanobservedvalue
	end
	
	estimatedNumerator = estimatedRatio .* denominator # these are for instance the estimated losses for a LR model
	
	estimatedNumeratorForStats = zeros(obs)
    indicatedRelativityForApplyTree_reused = zeros(obs)
	reused_fitted_leafnr_vector = zeros(Int, obs)
    sortvec_reused_trn_only = zeros(Int, length(trnidx))
    intVarsUsed = Array{Array{Array{Int,1},1}}(undef, 0);sizehint!(intVarsUsed, iterations)
	inds_considered = Array{Array{Int,1}}(undef, 0);sizehint!(inds_considered, iterations)
	
	res = Array{Union{Leaf{T_Uint8_or_UInt16},Node{UInt8},Node{UInt16}}}(undef, iterations)
	if prroduceEstAndLeafMatrices
		est_matrix = Array{Float64}(undef, obs, iterations + 1)
		est_matrixFromScores = copy(est_matrix)
		MatrixOfLeafNumbers = Array{Int}(undef, obs, iterations + 1)
		MatrixOfLeafNumbers[:,1] .= 0
		est_matrix[:,1] = copy(transpose(estimatedRatio))
		est_matrixFromScores[:,1] = copy(transpose(estimatedRatio))		
	else
		est_matrix = Array{Float64}(undef, 0, 0)
		est_matrixFromScores = copy(est_matrix)
		MatrixOfLeafNumbers = Array{Int}(undef, 0, 0)
	end
	
	scores = zeros(Int, obs)	
	estFromScores = zeros(obs)
	vectorOfLeafArrays = Array{Array{Leaf{T_Uint8_or_UInt16},1}}(undef, iterations + 1)
	vectorOfLeafArrays[1] = Array{Leaf{T_Uint8_or_UInt16}}(undef, 0) # the first entry is not defined		
	
    p = ProgressMeter.Progress(iterations, 2, "Progress of Boosting Model:") # minimum update interval: x seconds (2)
    ((adaptiveLearningRate >= 1.0) || (adaptiveLearningRate <= 0.0)) ? (BoolAdaptiveLearningRate = false) : (BoolAdaptiveLearningRate = true)
	# currently not supported
		# BoolAdaptiveLearningRate&&error("DTM: this is not yet updated to reflect numerator/denominator objectives!!")
	# this header needs to be in line with the output of the function createTrnValStatsForThisIteration
	statsPerIteration = global_statsperiter_header	
	currentRelativity = Array{Float64}(undef, obs)
	scoreBandLabels = createScorebandsUTF8List(sett.scorebandsstartingpoints, sett.nScores, addtotal=true)
	# this row is only here such that the variables are initialized outside the for loop
	local fristRowOfThisTable, obsPerScore, vectorWeightPerScore, numPerScore, denomPerScore, scores, scoresVAL, cumulativeStatsPerScoreBand, maxRawRelativityPerScoreSorted, rawObservedRatioPerScore, MAPPINGSmoothedEstimatePerScore

	# the following two lines are for the initialization in case any subsampling (boostrapping) of data is done
	# the idea is to preallocate many of the arrays which are used for sampling (and to save allocation/memory space)
		sampleSizeCanBeNEGATIVE = convert(Int, round(sett.subsampling_prop * obstrn))
		if abs(sampleSizeCanBeNEGATIVE) < 1
			sampleSizeCanBeNEGATIVE = convert(Int, sign(sett.subsampling_prop))
		end
		abssampleSize, sampleVector = initBootstrapSample(sampleSizeCanBeNEGATIVE) # ,trnidx,validx)

        # local stats,estimatedRatioUnsmoothed,estimatedRatioSmoothed,estimateRatioFromRelativities,resultingBT
        
    for iter = 1:iterations
        # currently disable timeConditional as it leads to a segfault (see https://github.com/JuliaLang/julia/issues/28536) \August 15, 2018
        # notably this segfault did not appear on 0.7beta2 (and earlier) where all tests ran perfectly fine (see travis)
        # @timeConditional(showTimeUsedByEachIteration, begin
		# Build Iteration iter
			current_mdf = moderationvector[min(iter, size(moderationvector, 1))] # mdf remains constant if a vector of size 1 was provided (instead of a vector of size iterations)			
			
			# @show iter,"here" 
			# @show size(trnidx)
			if showTimeUsedByEachIteration
			# todo tbd: this can be done more efficiently: We can avoid the tmpTree object (also we do not need to call "sometreesettings" every time!, improve this
				@time 	tmpTree = sample_data_and_build_tree!(trnidx, validx, indicatedRelativityForApplyTree_reused, candMatWOMaxValues, mappings, deepcopy(sett), actualNumerator, estimatedNumerator, weight, features, sampleSizeCanBeNEGATIVE, abssampleSize, sampleVector, T_Uint8_or_UInt16)
			else
						tmpTree = sample_data_and_build_tree!(trnidx, validx, indicatedRelativityForApplyTree_reused, candMatWOMaxValues, mappings, deepcopy(sett), actualNumerator, estimatedNumerator, weight, features, sampleSizeCanBeNEGATIVE, abssampleSize, sampleVector, T_Uint8_or_UInt16)
			end
			# indicatedRelativity
			if minimum(indicatedRelativityForApplyTree_reused) < 0.0
				# todo/tbd check this: if there are negative indicated relativites the boosting approach is not meaningful.
				@warn("The boosting model encountered negative indicated realtivities for some leaves.")
				@warn("This should not happen if you have only positive values for numerator and denominator.")
				error("Abort.")
			end
			# todo/tbd check this: what if a fitted value in a leaf is 0? what will the indicatedRelativityForApplyTree_reused be? will the boosting work as expected?
			push!(intVarsUsed, tmpTree.intVarsUsed)
			push!(inds_considered, tmpTree.inds_considered)
			res[iter] = tmpTree.rootnode # we do not need the full tree (with all its settings and so on) here, todo/tbd (someone else should) check if this is ok and creates no issues
			vectorOfLeafArrays[iter + 1] = create_leaves_array(res[iter])
			# vectorOfRulePathsToLeavesArrays[iter+1]=[x.rule_path for x in vectorOfLeafArrays[iter+1]]
			# Apply this iteration to the validation data set						
			apply_tree_by_leaf!(indicatedRelativityForApplyTree_reused, reused_fitted_leafnr_vector, validx, res[iter], features)
			# update_part_of_this_vector!(validx,indicatedRelativityForApplyTree_reused,indicatedRelativity)						
        if prroduceEstAndLeafMatrices
            # todo/tbd this could possibly be optimized. We are calling apply_tree_by_leaf! twice (if prroduceEstAndLeafMatrices==true), once for val and once for trn
            apply_tree_by_leaf!(indicatedRelativityForApplyTree_reused, reused_fitted_leafnr_vector, trnidx, res[iter], features)
            MatrixOfLeafNumbers[:,iter + 1] = reused_fitted_leafnr_vector # leaf_numbers(vectorOfLeafArrays[1+iter],obs) #TODO / TBD adjust this for subsampling. this may not work properly, as each tree will use a different amount of data AND DIFFERENT OBSERVATIONS! thus we cannot rely on the *.idx field of the leaves!			
        end			
			# Moderate the estimate, this is done for BOTH trn and val!
			_moderate!(estimatedRatio, indicatedRelativityForApplyTree_reused, current_mdf)
			if prroduceEstAndLeafMatrices
				write_column!(est_matrix, iter + 1, estimatedRatio)
			end
			update_current_rels!(currentRelativity, estimatedRatio, trn_meanobservedvalue)
			# update estimatedNumerator
			estimatedNumerator = estimatedRatio .* denominator
		# Derive Scores for this iteration #NOTE (tbd/todo: keep in mind) For large datasets (>5m rows) the sorting in construct scores may become a dominant part of the algorithm
			maxRawRelativityPerScoreSorted, MAPPINGSmoothedEstimatePerScore, vectorWeightPerScore, obsPerScore, rawObservedRatioPerScore, numPerScore, denomPerScore, nscoresPotentiallyReduced = constructANDderiveScores!(trnidx, validx, sortvec_reused_trn_only, estimatedRatio, currentRelativity, actualNumerator, denominator, weight, trn_meanobservedvalue, iter, sett)
			update_and_derive_scores!(scores, maxRawRelativityPerScoreSorted, currentRelativity)
			# todo, this could probably be improved
			# if we avoid creating small negative values in the function constructANDderiveScores altogether, then the following two function calls are obsolete
			correctSmallNegativeValues!(MAPPINGSmoothedEstimatePerScore)
			correctSmallNegativeValues!(rawObservedRatioPerScore)
		# Define estimated Numerator		
			sett.smoothEstimates == "0" ? referenceForNumEstimates = rawObservedRatioPerScore : referenceForNumEstimates = MAPPINGSmoothedEstimatePerScore
			update_mapped_estFromScores!(estFromScores, referenceForNumEstimates, scores)			
			
			if prroduceEstAndLeafMatrices
				write_column!(est_matrixFromScores, iter + 1, estFromScores)
			end
		# Derive Cumulative Statistics
			 # todo tbd there is potential for optimization here			 
			 # ["rawRelativities","unsmoothedPerScore","smoothedPerScore"]				
			if sett.fitForStatsAndCharts == "rawRelativities"				
				selectedEstimatedRatioForStats = currentRelativity .* trn_meanobservedvalue
			else 
				selectedEstimatedRatioForStats = estFromScores::Vector{Float64} # smooth or unsmooth has been defined/"selected" above already
			end
        for ij = 1:length(estimatedNumeratorForStats) # a.=b.*c seemed to allocate memory ..... (see memhist of June 24 2018)
            @inbounds  estimatedNumeratorForStats[ij] = selectedEstimatedRatioForStats[ij] * denominator[ij]
        end
			
			statsThisIteration, singleRowWithKeyMetrics, columnOfRelativityTrn = createTrnValStatsForThisIteration(scoreBandLabels, iter, sett.scorebandsstartingpoints, view(actualNumerator, trnidx), view(denominator, trnidx), view(weight, trnidx), view(selectedEstimatedRatioForStats, trnidx), view(scores, trnidx), view(actualNumerator, validx), view(denominator, validx), view(weight, validx), view(selectedEstimatedRatioForStats, validx), view(scores, validx), sett)			
			statsPerIteration = vcat(statsPerIteration, singleRowWithKeyMetrics)
			if iter == 1
				fristRowOfThisTable = sett.iterations + 2 + empty_rows_after_iteration_stats
				cumulativeStatsPerScoreBand = copy(statsThisIteration)
			else
				fristRowOfThisTable += size(statsThisIteration, 1)
				cumulativeStatsPerScoreBand = vcat(cumulativeStatsPerScoreBand, statsThisIteration)
			end
			# Add the respective chart
				loc = string(excelLetter(size(cumulativeStatsPerScoreBand, 2) + 2), fristRowOfThisTable)				
				thischart = defineRelativityChart(nameOfModelStatisticsSheet, nameOfModelStatisticsSheet, loc, size(statsThisIteration, 1) - 2, columnOfRelativityTrn, fristRowOfThisTable, xscale=2.0)
				push!(xlData.charts, deepcopy(thischart))
		# end	) #end of @timeConditional
		if showProgressBar_time # TimeUsedByEachIteration
			ProgressMeter.next!(p)
		end
    end
	# GC.enable(true)	#to check: is this needed? is there any benefit to it?
	# todo/tbd: note that we have at least three different estimates should we add all of them to the output?
	
	# one of these vectors is actually identical to estFromScores and would not need to be recalculated
	estimatedRatioUnsmoothed = map_these_values(rawObservedRatioPerScore, scores)	
	estimatedRatioSmoothed = map_these_values(MAPPINGSmoothedEstimatePerScore, scores)	
	estimateRatioFromRelativities = currentRelativity .* trn_meanobservedvalue		
	
	# fitForStatsAndCharts is an element of ["rawRelativities","unsmoothedPerScore","smoothedPerScore"]	
	if sett.fitForStatsAndCharts == "rawRelativities"		
		selectedEstimatedRatioForStats = estimateRatioFromRelativities
	end
	if sett.fitForStatsAndCharts == "unsmoothedPerScore"
		selectedEstimatedRatioForStats = estimatedRatioUnsmoothed
	end
	if sett.fitForStatsAndCharts == "smoothedPerScore"
		selectedEstimatedRatioForStats = estimatedRatioSmoothed
	end	
	estimatedNumeratorForStats = selectedEstimatedRatioForStats .* denominator
		
	# add empty lines
	statsPerIteration = vcat(statsPerIteration, repeat([""], empty_rows_after_iteration_stats, size(statsPerIteration, 2)))
	# Add Charts which show stats for all iterations
		add_iteration_charts!(xlData, sett, 1)
		size(statsPerIteration, 2) > size(cumulativeStatsPerScoreBand, 2) ? cumulativeStatsPerScoreBand = hcat(cumulativeStatsPerScoreBand, repeat([""], size(cumulativeStatsPerScoreBand, 1), size(statsPerIteration, 2) - size(cumulativeStatsPerScoreBand, 2))) : statsPerIteration = hcat(statsPerIteration, repeat([""], size(statsPerIteration, 1), -size(statsPerIteration, 2) + size(cumulativeStatsPerScoreBand, 2)))
		stats = vcat(statsPerIteration, cumulativeStatsPerScoreBand)
# Create Score List
	scoreheader = ["Score" "Row Count" "Weight" "Numerator" "Denominator" "Maximal Raw Relativity" "Observed Ratio" "Smoothed Fit"]	
	scoreMatrix = vcat(scoreheader, hcat(collect(1:size(maxRawRelativityPerScoreSorted, 1)), obsPerScore, vectorWeightPerScore, numPerScore, denomPerScore, maxRawRelativityPerScoreSorted, rawObservedRatioPerScore, MAPPINGSmoothedEstimatePerScore))
	scoreMatrix = hcat(scoreMatrix, repeat([""], size(scoreMatrix, 1), size(stats, 2) - size(scoreMatrix, 2)))
# Add Score Chart
	thischart = defineScoreChart(nameOfScoresSheet, nameOfScoresSheet, "B6", size(scoreMatrix, 1) - 1, 1, 3, 7, 8)
	push!(xlData.charts, deepcopy(thischart))
# Create univariate graphs
		predictorsData, predictorCharts = createPredictorData(trnidx, nameOfpredictorsSheet, mappings, candMatWOMaxValues, sett, scores, actualNumerator, denominator, weight, selectedEstimatedRatioForStats, features)
	# push predictorCharts into Excelobject
		push!(xlData.charts, predictorCharts...)
		predictorsSheet = ExcelSheet(nameOfpredictorsSheet,  DataFrame(predictorsData,:auto))
	# Populate ExcelSheet object with data
		statsSheet = ExcelSheet(nameOfModelStatisticsSheet, DataFrame(stats,:auto))
		overviewSheet = ExcelSheet(nameOfOverviewSheet, DataFrame())
		scoreMatrixDF = DataFrame(scoreMatrix,:auto)
		scoresSheet = ExcelSheet(nameOfScoresSheet, scoreMatrixDF)
		modelsettingsSheet = ExcelSheet(nameOfSettingsSheet, DataFrame( writeAllFieldsToArray(sett),:auto))
# Create two way validation charts
	twoWayValidation, validationCharts = createTwoWayValidationCharts(trnidx, validx, nameOfValidationSheet, scoreBandLabels, mappings, candMatWOMaxValues, sett, scores, actualNumerator, denominator, weight, selectedEstimatedRatioForStats, features)
	validationSheet = ExcelSheet(nameOfValidationSheet, DataFrame(twoWayValidation,:auto))
	push!(xlData.charts, validationCharts...)
# Define Excel
	xlData.sheets = [modelsettingsSheet,overviewSheet,statsSheet,scoresSheet,predictorsSheet,validationSheet]
# add boolNumeratorStats if requested
	if sett.boolNumeratorStats
		@warn("tariff est stats are not updated for DTM2: this needs review")
			errors_num_estimates = addTariffEstimationStatsAndGraphs!(xlData, trnidx, validx, actualNumerator, estimatedRatioUnsmoothed, estimatedRatioSmoothed, estimateRatioFromRelativities, est_matrix, est_matrixFromScores)
		# attach errors per iteration to modelStatistics== statsSheet (Excel file)
			statsSheetINDEX = findall(in([nameOfModelStatisticsSheet]), [x.name for x in xlData.sheets])
			@assert length(statsSheetINDEX) == 1 "Error no statsheet was found in Excel data (when trying to attach Numerator Errors to Modelstatistics)"
			statsSheetINDEX2 = statsSheetINDEX[1]
			errors_num_estimates = vcat(errors_num_estimates, repeat([""], size(stats, 1) - size(errors_num_estimates, 1), size(errors_num_estimates, 2)))
			stats = hcat(stats, errors_num_estimates)
			xlData.sheets[statsSheetINDEX2] = ExcelSheet(nameOfModelStatisticsSheet, DataFrame( stats,:auto))
	end
	z = deepcopy(stats[1:sett.iterations + 1,:])
	modelstats = DataFrame(convert(Array{Float64,2}, z[2:size(z, 1),:]),:auto)
	DataFrames.rename!(modelstats, Symbol[Symbol(x) for x in view(z, 1, :)])
# resulting BT
	# create trnidx
	@assert issorted(trnidx)
	trnidx_one_zero_full_length = map(x->UInt8(length(searchsorted(trnidx, x))), 1:length(scores))	
	fp = get_feature_pools(features)
	resultingBT = BoostedTree(res, sett, intVarsUsed, candMatWOMaxValues, mappings, inds_considered, actual_moderationvector, scores, currentRelativity, maxRawRelativityPerScoreSorted, trn_meanobservedvalue, startAtMean, MAPPINGSmoothedEstimatePerScore, rawObservedRatioPerScore, est_matrix, modelstats, xlData, trnidx_one_zero_full_length, fp)
	return xlData, estimatedRatio, MatrixOfLeafNumbers, vectorOfLeafArrays, rawObservedRatioPerScore, est_matrixFromScores, stats, estimatedRatioUnsmoothed, estimatedRatioSmoothed, estimateRatioFromRelativities, resultingBT
end

