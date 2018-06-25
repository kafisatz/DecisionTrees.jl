#Bagging Algorithm
function bagged_tree(dtmtable::DTMTable,sett::ModelSettings)		
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
	actual_moderationvector=Array{Float64}(undef,iterations)
	fill!(actual_moderationvector,deepcopy(moderationvector[1]))
else
	actual_moderationvector=deepcopy(moderationvector)
end
if BoolStartAtMean
	estimatedRatio=ones(length(weight)).*trn_meanobservedvalue
else
	error("DTM: This is currently not working: please set BoolStartAtMean=true)")
end

#obs,obsVAL,trn_meanobservedvalue,empty_rows_after_iteration_stats,showTimeUsedByEachIteration,chosen_apply_tree_fn,BoolStartAtMean,adaptiveLearningRate,moderationvector,iterations,nameOflistofpredictorsSheet,nameOfpredictorsSheet,nameOfModelStatisticsSheet,nameOfScoresSheet,nameOfOverviewSheet,nameOfSettingsSheet,nameOfValidationSheet,xlData,showProgressBar_time=initSettingsWhichAreTheSameForBoostingAndBagging(actualNumerator,actualNumeratorVAL,denominator,sett)
iterationsPerCore=calcIterationsPerCore(sett.niter,Distributed.Distributed.nprocs())
sampleSizeForEachTreeCanBeNEG=convert(Int,round(sett.subsampling_prop*length(trnidx)))
if abs(sampleSizeForEachTreeCanBeNEG)<1
	sampleSizeForEachTreeCanBeNEG=convert(Int,sign(sett.subsampling_prop))
end

vecTreesWErrs = Distributed.@distributed (vcat) for ii=1:Distributed.Distributed.nprocs()
	#create_bagged_trees(iterationsPerCore[ii],sampleSizeForEachTreeCanBeNEG,trn_meanobservedvalue,candMatWOMaxValues,mappings,sett,actualNumerator,denominator,weight,numfeatures,charfeatures)		
	#sample_data_and_build_tree!(trnidx,validx,indicatedRelativityForApplyTree_reused,candMatWOMaxValues,mappings,deepcopy(sett),actualNumerator,estimatedNumerator,weight,features,sampleSizeCanBeNEGATIVE,abssampleSize,sampleVector,T_Uint8_or_UInt16)
	create_bagged_trees(iterationsPerCore[ii],trnidx,sampleSizeForEachTreeCanBeNEG,trn_meanobservedvalue,candMatWOMaxValues,mappings,sett,actualNumerator,denominator,weight,features)	
end
		
	#Derive cumulative estimates for each iteration for trn and val data
		#treesOnly=[x.tree for x in vecTreesWErrs]
		#warn("think about the choice of goodness of fit criterium here!")
		weightPerTree=defineWeightPerTree(vecTreesWErrs,sett.baggingWeightTreesError)		
		#this header needs to be in line with the output of the function createTrnValStatsForThisIteration			
			statsPerIteration=global_statsperiter_header
			scoreBandLabels=createScorebandsUTF8List(sett.scorebandsstartingpoints,sett.nscores,addtotal=true)
		#this row is only here such that the variables are initialized outside the for loop
		local	obsPerScore,vectorWeightPerScore,numPerScore,denomPerScore,scores,scoresVAL,cumulativeStatsPerScoreBand,relativities,relativitiesVAL,maxRawRelativityPerScoreSorted,rawObservedRatioPerScore,MAPPINGSmoothedEstimatePerScore
	
		estimatedRatio=ones(obs).*trn_meanobservedvalue	
		currentEstimateOfIndividualTree=zeros(obs)
		leafNumbersThisTree=zeros(Int,obs)
		currentRelativity=ones(obs)
		sortvec_reused_trn_only=zeros(Int,length(trnidx))
		scores=zeros(Int,obs)
		estFromScores=zeros(obs)		
				
		leaves_of_tree=Array{Array{Leaf{T},1}}(undef,iterations)
		
		if boolProduceEstAndLeafMatrices
			est_matrix=Array{Float64}(undef,obs,iterations+1)
			est_matrixFromScores=copy(est_matrix)
			MatrixOfLeafNumbers=Array{Int}(undef,obs,iterations+1)
			MatrixOfLeafNumbers[:,1]=0
			est_matrix[:,1]=copy(transpose(estimatedRatio))
			est_matrixFromScores[:,1]=copy(transpose(estimatedRatio))		
		else
			est_matrix=Array{Float64}(undef,0,0)
			est_matrixFromScores=copy(est_matrix)
			MatrixOfLeafNumbers=Array{Int}(undef,0,0)
		end
		
    vectorOfLeafArrays=Array{Array{Leaf{T},1}}(undef,iterations+1)
		vectorOfLeafArrays[1]=Array{Leaf{T}}(undef,0) #the first entry is not defined
		p = ProgressMeter.Progress(iterations, 5, "Progress of Bagging Model:") # minimum update interval: 5 second    
		cumulativeWeight=0.0
		local fristRowOfThisTable
	for iter=1:iterations
	error("BK: TODO, need to port the functionality regarding sett.fitForStatsAndCharts to bagging")
		#fit current tree			
			vectorOfLeafArrays[iter+1]=create_leaves_array(vecTreesWErrs[iter].tree)	
			#note the individual trees were constructed on a part of the training data, here we apply them to the full training data set
			apply_tree_by_leaf!(currentEstimateOfIndividualTree,leafNumbersThisTree,trnidx,vecTreesWErrs[iter].tree.rootnode,features)		
			apply_tree_by_leaf!(currentEstimateOfIndividualTree,leafNumbersThisTree,validx,vecTreesWErrs[iter].tree.rootnode,features)		
			#currentEstimateOfIndividualTree[:],leafNumbersThisTree[:]=Core.eval(parse(sett.chosen_apply_tree_fn))(vectorOfRulePathsToLeavesArrays[iter+1],vecTreesWErrs[iter].tree.rootnode,numfeatures,charfeatures)
			#vectorOfLeafNumbersTrn[:,iter+1]=leafNumbersThisTree			
		#update weights
			wi=weightPerTree[iter]
			cumulativeWeight += wi
		#derive current cumulative estimate		
			estimatedRatio.= ((estimatedRatio.*cumulativeWeight).+(currentEstimateOfIndividualTree.*wi))./(cumulativeWeight+wi)			
			currentRelativity.=estimatedRatio./trn_meanobservedvalue

			if boolProduceEstAndLeafMatrices
				est_matrix[:,iter+1]=estimatedRatio
				MatrixOfLeafNumbers[:,iter+1]=leafNumbersThisTree #leaf_numbers(vectorOfLeafArrays[1+iter],obs) #TODO / TBD adjust this for subsampling. this may not work properly, as each tree will use a different amount of data AND DIFFERENT OBSERVATIONS! thus we cannot rely on the *.idx field of the leaves!			
			end	
		#Derive Scores for this iteration		
		maxRawRelativityPerScoreSorted,MAPPINGSmoothedEstimatePerScore,vectorWeightPerScore,obsPerScore,rawObservedRatioPerScore,numPerScore,denomPerScore,nscoresPotentiallyReduced				=constructANDderiveScores!(trnidx,validx,sortvec_reused_trn_only,estimatedRatio,currentRelativity,actualNumerator,denominator,weight,trn_meanobservedvalue,iter,sett)
		update_and_derive_scores!(scores,maxRawRelativityPerScoreSorted,currentRelativity)			
	#Define estimated Numerator
		sett.smoothEstimates=="0" ? referenceForNumEstimates=rawObservedRatioPerScore : referenceForNumEstimates=MAPPINGSmoothedEstimatePerScore
		update_mapped_estFromScores!(estFromScores,referenceForNumEstimates,scores)
		#estFromScores=map(x->referenceForNumEstimates[x],scores)	
		if boolProduceEstAndLeafMatrices
			write_column!(est_matrixFromScores,iter+1,estFromScores)
		end
		#	maxRawRelativityPerScoreSorted,MAPPINGSmoothedEstimatePerScore,vectorWeightPerScore,obsPerScore,rawObservedRatioPerScore,numPerScore,denomPerScore,scores,uniqueRelativitiesSorted,nscoresPotentiallyReduced=constructANDderiveScores(estimatedRatio,currentRelativity,sett.nscores,actualNumerator,denominator,weight,trn_meanobservedvalue,iter,sett.smoothEstimates,sett.print_details)
		#Apply this iteration to the validation data set
			#currentEstimateOfIndividualTreeVAL[:],leafNumbersThisTreeVAL[:]=Core.eval(parse(sett.chosen_apply_tree_fn))(vectorOfLeafArrays[iter+1],vecTreesWErrs.tree[iter.rootnode,numfeaturesVAL,charfeaturesVAL)
		#derive current cumulative estimate
			#estimatedRatioVAL[:]=(estimatedRatioVAL.*cumulativeWeight+currentEstimateOfIndividualTreeVAL.*wi)./(cumulativeWeight+wi)
			#currentRelativityVAL[:]=estimatedRatioVAL./trn_meanobservedvalue
			#est_matrixVAL[:,iter+1]=estimatedRatioVAL
			#scoresVAL=derive_scores(maxRawRelativityPerScoreSorted,currentRelativityVAL) 		
		#Define estimated Numerator	
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
		
			#sett.smoothEstimates=="0" ? referenceForNumEstimates=rawObservedRatioPerScore : referenceForNumEstimates=MAPPINGSmoothedEstimatePerScore		
			#estimatedNumeratorDerivedFromScores=[referenceForNumEstimates[x] for x in scores]
			#estimatedNumeratorDerivedFromScoresVAL=[referenceForNumEstimates[x] for x in scoresVAL]
			#est_matrixFromScores[:,iter+1]=estimatedNumeratorDerivedFromScores[:]
			#est_matrixFromScoresVAL[:,iter+1]=estimatedNumeratorDerivedFromScoresVAL[:]
			if showProgressBar_time #TimeUsedByEachIteration
				ProgressMeter.next!(p)
			end
	end
	#todo/tbd: a large part of this code is redundant to boosting, this should be calucalted in a separate function somewhere..... (such that we do not need to maintain two functions with the same code!)

	estimateUnsmoothed=map_these_values(rawObservedRatioPerScore,scores)
	estimateSmoothed=map_these_values(MAPPINGSmoothedEstimatePerScore,scores)	
	estimateFromRelativities = currentRelativity.*trn_meanobservedvalue
	finalEstimateForCharts = sett.smoothEstimates=="0" ? estimateSmoothed : estimateUnsmoothed	
	#add empty lines
	statsPerIteration=vcat(statsPerIteration,repeat([""],empty_rows_after_iteration_stats,size(statsPerIteration,2)))
	#Add Charts which show stats for all iterations
		add_iteration_charts!(xlData,sett,1)
		size(statsPerIteration,2)>size(cumulativeStatsPerScoreBand,2) ? cumulativeStatsPerScoreBand=hcat(cumulativeStatsPerScoreBand,repeat([""],size(cumulativeStatsPerScoreBand,1),size(statsPerIteration,2)-size(cumulativeStatsPerScoreBand,2))) : statsPerIteration=hcat(statsPerIteration,repeat([""],size(statsPerIteration,1),-size(statsPerIteration,2)+size(cumulativeStatsPerScoreBand,2)))
		stats=vcat(statsPerIteration,cumulativeStatsPerScoreBand)
#Create Score List
	scoreheader=["Score" "Row Count" "Weight" "Numerator" "Denominator" "Maximal Raw Relativity" "Observed Ratio" "Smoothed Fit"]	
	scoreMatrix=vcat(scoreheader,hcat(collect(1:size(maxRawRelativityPerScoreSorted,1)),obsPerScore,vectorWeightPerScore,numPerScore,denomPerScore,maxRawRelativityPerScoreSorted,rawObservedRatioPerScore,MAPPINGSmoothedEstimatePerScore))
	scoreMatrix=hcat(scoreMatrix,repeat([""],size(scoreMatrix,1),size(stats,2)-size(scoreMatrix,2)))
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

	z=deepcopy(stats[1:sett.niter+1,:])
	modelstats=DataFrame(convert(Array{Float64,2},z[2:size(z,1),:]))
	DataFrames.names!(modelstats,Symbol[Symbol(x) for x in view(z,1,:)])
#resulting Bagged Tree, construct Bagged Tree
	intVarsUsed=[x.tree.intVarsUsed for x in vecTreesWErrs]
	inds_considered=[x.tree.inds_considered for x in vecTreesWErrs]
	nodesOnly=[x.tree.rootnode for x in vecTreesWErrs]
	errs=[x.err for x in vecTreesWErrs]

	#create trnidx
	@assert issorted(trnidx)
	trnidx_one_zero_full_length=map(x->UInt8(length(searchsorted(trnidx,x))),1:length(scores))		
	fp=get_feature_pools(features)
	resultingBaggedTree=BaggedTree(nodesOnly,weightPerTree,errs,sett,intVarsUsed,candMatWOMaxValues,mappings,inds_considered,scores,currentRelativity,maxRawRelativityPerScoreSorted,trn_meanobservedvalue,BoolStartAtMean,MAPPINGSmoothedEstimatePerScore,rawObservedRatioPerScorexx,est_matrix,modelstats,trnidx_one_zero_full_length,fp)
		
#resulting Bagged Tree
return xlData,estimatedRatio,MatrixOfLeafNumbers,vectorOfLeafArrays,rawObservedRatioPerScore,est_matrixFromScores,stats,estimateUnsmoothed,estimateSmoothed,estimateFromRelativities,resultingBaggedTree
end

function create_bagged_trees(itr::Int,trnidx_which_should_only_be_used_once_here::Vector{Int},sampleSizeCanBeNEGATIVE::Int,trn_meanobservedvalue::Float64,candMatWOMaxValues::Array{Array{Float64,1},1},mappings::Array{Array{String,1},1},sett::ModelSettings,numerator::Array{Float64,1},denominator::Array{Float64,1},weight::Array{Float64,1},features::DataFrame)
#function create_bagged_trees(iterationsPerCore[ii],sampleSizeForEachTreeCanBeNEG,trn_meanobservedvalue,candMatWOMaxValues,mappings,sett,actualNumerator,denominator,weight,features)	
	res=Array{TreeWithErrorStats}(itr)
	resFitted=Array{Array{Float64,1}}(1)
		
	abssampleSize,sampleVector=initBootstrapSample(sampleSizeCanBeNEGATIVE) 
	#abssampleSize,sampleVector,num,denom,w,numf,charf,ooBagnum,ooBagdenom,ooBagw,ooBagnumf,ooBagcharf,ooBagsize=initBoostrapSample(sampleSizeCanBeNEGATIVE,numerator,denominator,weight,charfeatures,numfeatures)
			
	fitted_values_all_data_this_vector_is_modified_by_build_tree=zeros(length(numerator))
	reused_fitted_leafnr_vector=zeros(Int,length(numerator))
	estimatedNumFromRelativities=zeros(length(numerator))
	sortvec_reused_trn_only=zeros(Int,length(sampleVector))
	#modelling loop
	for i=1:itr
		#sample data			
		unusedSamplePart=sampleData!(trnidx_which_should_only_be_used_once_here,sampleSizeCanBeNEGATIVE,sampleVector)							
	#build tree	on subsample of data
			thistree=build_tree!(sampleVector,unusedSamplePart,candMatWOMaxValues,mappings,deepcopy(sett),numerator,denominator,weight,features,fitted_values_all_data_this_vector_is_modified_by_build_tree,T_Uint8_or_UInt16)
	#apply tree to ooBagSample
			leaves_of_tree=create_leaves_array(thistree.rootnode)
			apply_tree_by_leaf!(fitted_values_all_data_this_vector_is_modified_by_build_tree,reused_fitted_leafnr_vector,sampleVector,thistree.rootnode,features)
			if length(unusedSamplePart)>0
				apply_tree_by_leaf!(fitted_values_all_data_this_vector_is_modified_by_build_tree,reused_fitted_leafnr_vector,unusedSamplePart,thistree.rootnode,features)
			end
		#Determine Goodness of Fit
			#I think it is best to consider the goodness of fit per leaf here
			#out of bag performance should be measured here, however this is not coded yet (see also bagging!)
				#ooBagEstNumerator=ooBagEstimates.*ooBagdenom
				#fittedPerLeaf=Float64[x.fitted for x in leaves_of_tree]
				#ooBagcntperLeaf,ooBagsumnumeratorEST,ooBagsumnumerator,ooBagsumdenominator,ooBagsumweight=aggregateByLeafNr(leafNrooBag,ooBagEstNumerator,ooBagnum,ooBagdenom,ooBagw)
				relativities=fitted_values_all_data_this_vector_is_modified_by_build_tree./trn_meanobservedvalue
				estimatedRatio=fitted_values_all_data_this_vector_is_modified_by_build_tree
				resFitted[1]=fitted_values_all_data_this_vector_is_modified_by_build_tree			
			#derive scores				
			maxRawRelativityPerScoreSorted,MAPPINGSmoothedEstimatePerScore,vectorWeightPerScore,obsPerScore,rawObservedRatioPerScore,numPerScore,denomPerScore,nscoresPotentiallyReduced				=constructANDderiveScores!(sampleVector,
			unusedSamplePart,sortvec_reused_trn_only,estimatedRatio,relativities,numerator,denominator,weight,trn_meanobservedvalue,itr,sett)
			#old call constructScores(estimatedRatio,relativities,uniqueRelativitiesSorted,num,denom,w,trn_meanobservedvalue,ns,sett.nscores,sett.print_details)
			#uniqueRelativitiesSorted=unique(relativities)
			#sort!(uniqueRelativitiesSorted,alg=QuickSort)
			#ns=min(sett.nscores,size(uniqueRelativitiesSorted,1))
			#ns<nscores&&warn("There are only $(ns) distinct relativities. Number of scores will be limited to $(ns) instead of $(nscores).")		
			#aggregatedModelledRatioPerScore,maxRawRelativityPerScoreSorted,MAPPINGSmoothedEstimatePerScore,vectorWeightPerScore,obsPerScore,rawObservedRatioPerScore,numPerScore,denomPerScore=	constructScores(estimatedRatio,relativities,uniqueRelativitiesSorted,num,denom,w,trn_meanobservedvalue,ns,sett.nscores,sett.print_details)

		#apply tree to ooBagSample
			if size(unusedSamplePart,1)>0
				#subsampling is performed
				#warn("todo: PICK a good function to determine the goodness of FIT / determine out of bag performance here! which will determine the weight")
				#@show sampleSizeCanBeNEGATIVE
				#if sampleSizeCanBeNEGATIVE<0
				#	resize!(leafNrooBag,size(unusedSamplePart,1))
				#	resize!(ooBagEstimates,size(unusedSamplePart,1))
				#end
				#leaves_of_tree=create_leaves_array(thistree.rootnode)	
				#if sett.boolProduceEstAndLeafMatrices
				#	error("this has not yet been updated for bagging")
					#todo/tbd this could possibly be optimized. We are calling apply_tree_by_leaf! twice (if boolProduceEstAndLeafMatrices==true), once for val and once for trn
					#apply_tree_by_leaf!(fitted_values_all_data_this_vector_is_modified_by_build_tree,reused_fitted_leafnr_vector,sampleVector,thistree.rootnode,features)
					#MatrixOfLeafNumbers[:,iter+1]=copy(reused_fitted_leafnr_vector) #leaf_numbers(vectorOfLeafArrays[1+iter],obs) #TODO / TBD adjust this for subsampling. this may not work properly, as each tree will use a different amount of data AND DIFFERENT OBSERVATIONS! thus we cannot rely on the *.idx field of the leaves!			
				#end	
                #Determine Goodness of Fit
					#I think it is best to consider the goodness of fit per leaf here
					#ooBagEstNumerator=ooBagEstimates.*ooBagdenom
					
				#relativitiesVAL=ooBagEstimates./trn_meanobservedvalue
				#scoresVAL=derive_scores(maxRawRelativityPerScoreSorted,relativitiesVAL) 
				#todo tbd decide which estimate we should use here.....!
				estimatedNumFromRelativities.=indicatedRelativityForApplyTree_reused.*denominator
				#sett.smoothEstimates=="0" ? referenceForNumEstimates=rawObservedRatioPerScore : referenceForNumEstimates=MAPPINGSmoothedEstimatePerScore			
				#estimatedNumeratorDerivedFromScoresVAL=[referenceForNumEstimates[x] for x in scoresVAL]
				# would this make sense? aggregate_per_score(estimatedNumFromRelativitiesVal,ooBagnum,ooBagw,ooBagdenom)
				thisError=goodnessOfFit(unusedSamplePart,indicatedRelativityForApplyTree_reused,numerator,weight,denominator)
			else
			#no subsampling is performed
				thisError=ErrorStats() #Empty/Default Error Stats
			end
			res[i]=TreeWithErrorStats(deepcopy(thistree),deepcopy(thisError))
	end
	#return resFitted,res
	return res #todo tbd check this: currently we do not return fitted, but I think it is probably not needed anyway
end
