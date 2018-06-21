#Bagging of multiple boosted trees
function bagged_boosted_tree(mappings::Array{Array{String,1},1},candMatWOMaxValues::Array{Array{Float64,1},1},sett::ModelSettings,actualNumerator::Array{Float64,1},denominator::Array{Float64,1},weight::Array{Float64,1}, numfeatures::Array{pdaMod,1},charfeatures::Array{pdaMod,1},actualNumeratorVAL::Array{Float64,1},denominatorVAL::Array{Float64,1},weightVAL::Array{Float64,1}, numfeaturesVAL::Array{pdaMod,1},charfeaturesVAL::Array{pdaMod,1})
error("this will not work as it has not been updated for DecisionTrees.")
	@assert sett.boolRandomizeOnlySplitAtTopNode=true #this must be true because of the selection of features in the tree construction ->    num_inds,char_inds=randomFeatureSelection(settings.number_of_num_features,settings.number_of_char_features,settings.subsampling_features_prop)
	obs,obsVAL,meanobservedvalue,empty_rows_after_iteration_stats,showTimeUsedByEachIteration,chosen_apply_tree_fn,BoolStartAtMean,adaptiveLearningRate,moderationvector,iterations,nameOflistofpredictorsSheet,nameOfpredictorsSheet,nameOfModelStatisticsSheet,nameOfScoresSheet,nameOfOverviewSheet,nameOfSettingsSheet,nameOfValidationSheet,xlData,showProgressBar_time=initSettingsWhichAreTheSameForBoostingAndBagging(actualNumerator,actualNumeratorVAL,denominator,sett)
	cBB_niterBoosting=sett.cBB_niterBoosting
	cBB_niterBagging=sett.cBB_niterBagging

	iterationsPerCore=calcIterationsPerCore(cBB_niterBagging,nprocs())
	sampleSizeForEachTreeCanBeNEG=convert(Int,round(sett.subsampling_prop*length(denominator)))
	
	vecBoostedTrees = Distributed.@parallel (vcat) for ii=1:nprocs()		
	#	ii=1;vecTreesWErrs=
		create_Subsamples_and_BoostedTrees(iterationsPerCore[ii],sampleSizeForEachTreeCanBeNEG,meanobservedvalue,candMatWOMaxValues,mappings,sett,actualNumerator,denominator,weight,numfeatures,charfeatures)
    end
	
end

function function create_Subsamples_and_BoostedTrees(itr::Int,sampleSizeCanBeNEGATIVE::Int,meanobservedvalue::Float64,candMatWOMaxValues::Array{Array{Float64,1},1},mappings::Array{Array{String,1},1},originalSett::ModelSettings,numerator::Array{Float64,1},denominator::Array{Float64,1},weight::Array{Float64,1},numfeatures::Array{pdaMod,1},charfeatures::Array{pdaMod,1})
#NOTE: we want to avoid that the individual trees of each boosted tree select variables
#Sett sampling settings for individual tree
	sett=deepcopy(originalSett)
	fixednuminds,fixedcharinds=randomFeatureSelection(length(numfeatures),length(charfeatures),originalSett.subsampling_features_prop)
	sett.fixedcharinds=deepcopy(fixedcharinds)
	sett.fixednuminds=deepcopy(fixednuminds)
	sett.subsampling_features_prop=1.0
#init sampledata
	abssampleSize,sampleVector,num,denom,w,numf,charf,ooBagnum,ooBagdenom,ooBagw,ooBagnumf,ooBagcharf,ooBagsize=initBoostrapSample(sampleSizeCanBeNEGATIVE,numerator,denominator,weight,charfeatures,numfeatures)
	leafNrooBag=Array{Int}(ooBagsize)
	ooBagEstimates=Array{Float64}(ooBagsize)		
	
	for i=1:itr
	#sample data			
    error("this will not work as it has not been updated for DecisionTrees.")
		unusedSamplePart=sampleData!(sampleSizeCanBeNEGATIVE,numerator,denominator,weight,numfeatures,charfeatures,sampleVector,num,denom,w,numf,charf,ooBagnum,ooBagdenom,ooBagw,ooBagnumf,ooBagcharf)
	#build boosted tree			
		xlData,vectorOfLeafNumbersTrn,vectorOfLeafArrays,vectorOfRulePathsToLeavesArrays,rawObservedRatioPerScore,est_matrixFromScores,est_matrixFromScoresVAL,stats,scoresval,est_matrix_val,estimateUnsmoothedVal,estimateSmoothedVal,estimateFromRelativitiesVal,estimateUnsmoothedTrn,estimateSmoothedTrn,estimateFromRelativitiesTrn,resultEnsemble=boosted_tree(mappings,candMatWOMaxValues,sett,num,denom,w,numf,charf,ooBagnum,ooBagdenom,ooBagw,ooBagnumf,ooBagcharf)
	end
error(".... in the works....")
	return true
end

#@show "this is the end of the file bagged_boosted_tree.jl"