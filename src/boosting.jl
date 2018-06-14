#Boosting
function boosted_tree(dtmtable::DTMTable,sett::ModelSettings)		
	trnidx=dtmtable.trnidx
	validx=dtmtable.validx
	features=dtmtable.features
	weight=dtmtable.weight
	actualNumerator=dtmtable.numerator
	denominator=dtmtable.denominator
	mappings=dtmtable.mappings
	candMatWOMaxValues=dtmtable.candMatWOMaxValues
	obstrn,obsval,trn_meanobservedvalue,val_meanobservedvalue,trn_numtot,val_numtot,trn_denomtot,val_denomtot,empty_rows_after_iteration_stats,showTimeUsedByEachIteration,chosen_apply_tree_fn,BoolStartAtMean,adaptiveLearningRate,moderationvector,iterations,nameOflistofpredictorsSheet,nameOfpredictorsSheet,nameOfModelStatisticsSheet,nameOfScoresSheet,nameOfOverviewSheet,nameOfSettingsSheet,nameOfValidationSheet,xlData,showProgressBar_time,boolProduceEstAndLeafMatrices=initSettingsWhichAreTheSameForBoostingAndBagging(trnidx,validx,actualNumerator,denominator,sett)
	obs=obstrn+obsval
	current_error=0.0
    moderationfactor=moderationvector[1]
	if (size(moderationvector,1)!=iterations && size(moderationvector,1)!=1)
		error("DTM: This should not have happend: size(mdf_vector) is neither 1 nor #iterations")
    end
	if size(moderationvector,1)==1
		actual_moderationvector=Array{Float64}(iterations)
		fill!(actual_moderationvector,deepcopy(moderationvector[1]))
	else
		actual_moderationvector=deepcopy(moderationvector)
	end
	if BoolStartAtMean
		estimatedRatio=ones(length(weight)).*trn_meanobservedvalue
	else
		error("DTM: This is currently not working: please set BoolStartAtMean=true)")
	end

	indicatedRelativityForApplyTree_reused=zeros(obs)
	reused_fitted_leafnr_vector=zeros(Int,obs)
    sortvec_reused_trn_only=zeros(Int,length(trnidx))
    intVarsUsed=Array{Array{Array{Int,1},1}}(0);sizehint!(intVarsUsed,iterations)
	inds_considered=Array{Array{Int64,1}}(0);sizehint!(inds_considered,iterations)
	#intCharVarsUsed=Array{Array{Array{Int,1},1}}(0);sizehint!(intCharVarsUsed,iterations)
	#char_inds_considered=Array{Array{Int64,1}}(0);sizehint!(char_inds_considered,iterations)
	
	res=Array{Union{Leaf,Node{UInt8},Node{UInt16}}}(iterations)
	if boolProduceEstAndLeafMatrices
		est_matrix=Array{Float64}(obs,iterations+1)
		est_matrixFromScores=deepcopy(est_matrix)
		MatrixOfLeafNumbers=Array{Int}(obs,iterations+1)
		MatrixOfLeafNumbers[:,1]=0
		est_matrix[:,1]=copy(transpose(estimatedRatio))
		est_matrixFromScores[:,1]=copy(transpose(estimatedRatio))		
	else
		est_matrix=Array{Float64}(0,0)
		est_matrixFromScores=deepcopy(est_matrix)
		MatrixOfLeafNumbers=Array{Int}(0,0)
	end
	
	scores=zeros(Int,obs)	
	estFromScores=zeros(obs)
	#vectorOfRulePathsToLeavesArrays=Array{Array{Array{Rulepath,1},1}}(iterations+1)
	#vectorOfRulePathsToLeavesArrays[1]=Array{Array{Rulepath,1}}(0) #the first entry is not defined	
	vectorOfLeafArrays=Array{Array{Leaf,1}}(iterations+1)
	vectorOfLeafArrays[1]=Array{Leaf}(0) #the first entry is not defined		
	
    p = Progress(iterations, 2, "Progress of Boosting Model:") # minimum update interval: x seconds (2)
    ((adaptiveLearningRate>=1.0)||(adaptiveLearningRate<=0.0)) ? (BoolAdaptiveLearningRate=false) : (BoolAdaptiveLearningRate=true)
	#currently not supported
		BoolAdaptiveLearningRate&&error("DTM: this is not yet updated to reflect numerator/denominator objectives!!")
	#this header needs to be in line with the output of the function createTrnValStatsForThisIteration
	statsPerIteration=global_statsperiter_header	
	currentRelativity=Array{Float64}(obs)
	scoreBandLabels=createScorebandsUTF8List(sett.scorebandsstartingpoints,sett.nscores,addtotal=true)
	#this row is only here such that the variables are initialized outside the for loop
	local fristRowOfThisTable,obsPerScore,vectorWeightPerScore,numPerScore,denomPerScore,scores,scoresVAL,cumulativeStatsPerScoreBand,maxRawRelativityPerScoreSorted,rawObservedRatioPerScore,MAPPINGSmoothedEstimatePerScore

	#the following two lines are for the initialization in case any subsampling (boostrapping) of data is done
	#the idea is to preallocate many of the arrays which are used for sampling (and to save allocation/memory space)
		sampleSizeCanBeNEGATIVE=convert(Int,round(sett.subsampling_prop*obstrn))
		if abs(sampleSizeCanBeNEGATIVE)<1
			sampleSizeCanBeNEGATIVE=convert(Int,sign(sett.subsampling_prop))
		end
		abssampleSize,sampleVector=initBootstrapSample(sampleSizeCanBeNEGATIVE) #,trnidx,validx)	
	for iter=1:iterations
        @timeConditional(showTimeUsedByEachIteration, begin
		#Build Iteration iter
			#showTimeUsedByEachIteration&&println("Calculating iteration $(iter). Time: $(now())")
			current_mdf=moderationvector[min(iter,size(moderationvector,1))] #mdf remains constant if a vector of size 1 was provided (instead of a vector of size iterations)
			estimatedNumerator=estimatedRatio.*denominator #these are for instance the estimated losses for a LR model
			#@show estimatedNumerator[1:10]
			#@show mean(estimatedNumerator),sum(estimatedNumerator),mean(estimatedRatio),sum(estimatedRatio)
			if showTimeUsedByEachIteration
			#todo tbd: this can be done more efficiently: We can avoid the tmpTree object (also we do not need to call "sometreesettings" every time!, improve this
				@time 	tmpTree=sample_data_and_build_tree!(trnidx,validx,indicatedRelativityForApplyTree_reused,candMatWOMaxValues,mappings,deepcopy(sett),actualNumerator,estimatedNumerator,weight,features,sampleSizeCanBeNEGATIVE,abssampleSize,sampleVector)
			else
						tmpTree=sample_data_and_build_tree!(trnidx,validx,indicatedRelativityForApplyTree_reused,candMatWOMaxValues,mappings,deepcopy(sett),actualNumerator,estimatedNumerator,weight,features,sampleSizeCanBeNEGATIVE,abssampleSize,sampleVector)
			end
			#indicatedRelativity
			if (!sett.bINTERNALignoreNegRelsBoosting)&&(minimum(indicatedRelativityForApplyTree_reused)<0.0)
				#todo/tbd check this: if there are negative indicated relativites the boosting approach is not meaningful.
				warn("The boosting model encountered negative indicated realtivities for some leaves.")
				warn("This should not happen if you have only positive values for numerator and denominator.")
				error("Abort.")
			end
			#todo/tbd check this: what if a fitted value in a leaf is 0? what will the indicatedRelativityForApplyTree_reused be? will the boosting work as expected?
			push!(intVarsUsed,tmpTree.intVarsUsed)
			push!(inds_considered,tmpTree.inds_considered)
			res[iter]=tmpTree.rootnode #we do not need the full tree (with all its settings and so on) here, todo/tbd (someone else should) check if this is ok and creates no issues
			vectorOfLeafArrays[iter+1]=create_leaves_array(res[iter])
			#vectorOfRulePathsToLeavesArrays[iter+1]=[x.rule_path for x in vectorOfLeafArrays[iter+1]]
			#Apply this iteration to the validation data set						
			apply_tree_by_leaf!(indicatedRelativityForApplyTree_reused,reused_fitted_leafnr_vector,validx,res[iter],features)
			#update_part_of_this_vector!(validx,indicatedRelativityForApplyTree_reused,indicatedRelativity)						
            if boolProduceEstAndLeafMatrices
            #todo/tbd this could possibly be optimized. We are calling apply_tree_by_leaf! twice (if boolProduceEstAndLeafMatrices==true), once for val and once for trn
                    apply_tree_by_leaf!(indicatedRelativityForApplyTree_reused,reused_fitted_leafnr_vector,trnidx,res[iter],features)
                    MatrixOfLeafNumbers[:,iter+1]=reused_fitted_leafnr_vector #leaf_numbers(vectorOfLeafArrays[1+iter],obs) #TODO / TBD adjust this for subsampling. this may not work properly, as each tree will use a different amount of data AND DIFFERENT OBSERVATIONS! thus we cannot rely on the *.idx field of the leaves!			
            end			
			#Moderate the estimate, this is done for BOTH trn and val!
			_moderate!(estimatedRatio,indicatedRelativityForApplyTree_reused,current_mdf)
			if boolProduceEstAndLeafMatrices
				write_column!(est_matrix,iter+1,estimatedRatio)
			end
			update_current_rels!(currentRelativity,estimatedRatio,trn_meanobservedvalue)
		#Derive Scores for this iteration #NOTE (tbd/todo: keep in mind) For large datasets (>5m rows) the sorting in construct scores may become a dominant part of the algorithm
			maxRawRelativityPerScoreSorted,MAPPINGSmoothedEstimatePerScore,vectorWeightPerScore,obsPerScore,rawObservedRatioPerScore,numPerScore,denomPerScore,nscoresPotentiallyReduced=constructANDderiveScores!(trnidx,validx,sortvec_reused_trn_only,estimatedRatio,currentRelativity,actualNumerator,denominator,weight,trn_meanobservedvalue,iter,sett)
			update_and_derive_scores!(scores,maxRawRelativityPerScoreSorted,currentRelativity)			
		#Define estimated Numerator		
			sett.smoothEstimates=="0" ? referenceForNumEstimates=rawObservedRatioPerScore : referenceForNumEstimates=MAPPINGSmoothedEstimatePerScore
			update_mapped_estFromScores!(estFromScores,referenceForNumEstimates,scores)
			#estFromScores=map(x->referenceForNumEstimates[x],scores)	
			if boolProduceEstAndLeafMatrices
				write_column!(est_matrixFromScores,iter+1,estFromScores)
			end
		#Derive Cumulative Statistics
			 #todo tbd there is potential for optimization here			 
			statsThisIteration,singleRowWithKeyMetrics,columnOfRelativityTrn=createTrnValStatsForThisIteration(scoreBandLabels,iter,sett.scorebandsstartingpoints,view(actualNumerator,trnidx),view(denominator,trnidx),view(weight,trnidx),view(estFromScores,trnidx),view(scores,trnidx),view(actualNumerator,validx),view(denominator,validx),view(weight,validx),view(estFromScores,validx),view(scores,validx),sett)			
			statsPerIteration=vcat(statsPerIteration,singleRowWithKeyMetrics)
			if iter==1
				fristRowOfThisTable=sett.niter+2+empty_rows_after_iteration_stats
				cumulativeStatsPerScoreBand=copy(statsThisIteration)
			else
				fristRowOfThisTable+=size(statsThisIteration,1)
				cumulativeStatsPerScoreBand=vcat(cumulativeStatsPerScoreBand,statsThisIteration)
			end
			#Add the respective chart
				loc=string(excelLetter(size(cumulativeStatsPerScoreBand,2)+2),fristRowOfThisTable)				
				thischart=defineRelativityChart(nameOfModelStatisticsSheet,nameOfModelStatisticsSheet,loc,size(statsThisIteration,1)-2,columnOfRelativityTrn,fristRowOfThisTable,xscale=2.0)
				push!(xlData.charts,deepcopy(thischart))
		end	)
		if showProgressBar_time #TimeUsedByEachIteration
			next!(p)
		end
    end
	gc_enable(true)	
	#todo/tbd: note that we have at least three different estimates should we add all of them to the output?
	
	estimateUnsmoothed=map_these_values(rawObservedRatioPerScore,scores)
	estimateSmoothed=map_these_values(MAPPINGSmoothedEstimatePerScore,scores)	
	
	estimateFromRelativities = currentRelativity.*trn_meanobservedvalue
	finalEstimateForCharts = sett.smoothEstimates=="0" ? estimateSmoothed : estimateUnsmoothed
	
	#add empty lines
	statsPerIteration=vcat(statsPerIteration,repmat([""],empty_rows_after_iteration_stats,size(statsPerIteration,2)))
	#Add Charts which show stats for all iterations
		add_iteration_charts!(xlData,sett,1)
		size(statsPerIteration,2)>size(cumulativeStatsPerScoreBand,2) ? cumulativeStatsPerScoreBand=hcat(cumulativeStatsPerScoreBand,repmat([""],size(cumulativeStatsPerScoreBand,1),size(statsPerIteration,2)-size(cumulativeStatsPerScoreBand,2))) : statsPerIteration=hcat(statsPerIteration,repmat([""],size(statsPerIteration,1),-size(statsPerIteration,2)+size(cumulativeStatsPerScoreBand,2)))
		stats=vcat(statsPerIteration,cumulativeStatsPerScoreBand)
#Create Score List
	scoreheader=["Score" "Row Count" "Weight" "Numerator" "Denominator" "Maximal Raw Relativity" "Observed Ratio" "Smoothed Fit"]	
	scoreMatrix=vcat(scoreheader,hcat(collect(1:size(maxRawRelativityPerScoreSorted,1)),obsPerScore,vectorWeightPerScore,numPerScore,denomPerScore,maxRawRelativityPerScoreSorted,rawObservedRatioPerScore,MAPPINGSmoothedEstimatePerScore))
	scoreMatrix=hcat(scoreMatrix,repmat([""],size(scoreMatrix,1),size(stats,2)-size(scoreMatrix,2)))
#Add Score Chart
	thischart=defineScoreChart(nameOfScoresSheet,nameOfScoresSheet,"B6",size(scoreMatrix,1)-1,1,3,7,8)
	push!(xlData.charts,deepcopy(thischart))
#Create univariate graphs
		predictorsData,predictorCharts=createPredictorData(trnidx,nameOfpredictorsSheet,mappings,candMatWOMaxValues,sett,scores,actualNumerator,denominator,weight,finalEstimateForCharts,features)
	#push predictorCharts into Excelobject
		push!(xlData.charts,predictorCharts...)
		predictorsSheet=ExcelSheet(nameOfpredictorsSheet,convert(DataFrame,predictorsData))
	#Populate ExcelSheet object with data
		statsSheet=ExcelSheet(nameOfModelStatisticsSheet,convert(DataFrame,stats))
		overviewSheet=ExcelSheet(nameOfOverviewSheet,DataFrame())
		scoreMatrixDF=convert(DataFrame,scoreMatrix)
		scoresSheet=ExcelSheet(nameOfScoresSheet,scoreMatrixDF)
		modelsettingsSheet=ExcelSheet(nameOfSettingsSheet,convert(DataFrame,writeAllFieldsToArray(sett)))
#Create two way validation charts
	twoWayValidation,validationCharts=createTwoWayValidationCharts(trnidx,validx,nameOfValidationSheet,scoreBandLabels,mappings,candMatWOMaxValues,sett,scores,actualNumerator,denominator,weight,finalEstimateForCharts,features)
	validationSheet=ExcelSheet(nameOfValidationSheet,convert(DataFrame,twoWayValidation))
	push!(xlData.charts,validationCharts...)
#Define Excel
	xlData.sheets=[modelsettingsSheet,overviewSheet,statsSheet,scoresSheet,predictorsSheet,validationSheet]
#add boolTariffEstStats if requested
	if sett.boolTariffEstStats
		warn("tariff est stats are not updated for DTM2: this needs review")
			errors_num_estimates=addTariffEstimationStatsAndGraphs!(xlData,trnidx,validx,actualNumerator,estimateUnsmoothed,estimateSmoothed,estimateFromRelativities,est_matrix,est_matrixFromScores)
		#attach errors per iteration to modelStatistics== statsSheet (Excel file)
			statsSheetINDEX=findin([x.name for x in xlData.sheets],[nameOfModelStatisticsSheet])
			@assert length(statsSheetINDEX)==1 "Error no statsheet was found in Excel data (when trying to attach Numerator Errors to Modelstatistics)"
			statsSheetINDEX2=statsSheetINDEX[1]
			errors_num_estimates=vcat(errors_num_estimates,repmat([""],size(stats,1)-size(errors_num_estimates,1),size(errors_num_estimates,2)))
			stats=hcat(stats,errors_num_estimates)
			xlData.sheets[statsSheetINDEX2]=ExcelSheet(nameOfModelStatisticsSheet,convert(DataFrame,stats))
	end
	z=deepcopy(stats[1:sett.niter+1,:])
	modelstats=DataFrame(convert(Array{Float64,2},z[2:size(z,1),:]))
	names!(modelstats,Symbol[Symbol(x) for x in view(z,1,:)])
#resulting BT
	#create trnidx
	@assert issorted(trnidx)
	trnidx_one_zero_full_length=map(x->UInt8(length(searchsorted(trnidx,x))),1:length(scores))	
	fp=get_feature_pools(features)
	resultingBT=BoostedTree(res,sett,intVarsUsed,candMatWOMaxValues,mappings,inds_considered,actual_moderationvector,scores,currentRelativity,maxRawRelativityPerScoreSorted,trn_meanobservedvalue,BoolStartAtMean,MAPPINGSmoothedEstimatePerScore,rawObservedRatioPerScore,est_matrix,modelstats,xlData,trnidx_one_zero_full_length,fp)
	return xlData,estimatedRatio,MatrixOfLeafNumbers,vectorOfLeafArrays,rawObservedRatioPerScore,est_matrixFromScores,stats,estimateUnsmoothed,estimateSmoothed,estimateFromRelativities,resultingBT
end

