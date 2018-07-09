#susing Revise
#using Profile
#using JLD2
using CSV
#using DataFrames
using DecisionTrees  

#include(joinpath("..","test\\runtests.jl"))

datafile=string("data\\freMTPL2\\freMTPL2.csv")
@assert isfile(datafile);
@time fullData=CSV.read(datafile,rows_for_type_detect=100000,allowmissing=:none,categorical=false);

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

dtmtable,sett,dfprepped=prepare_dataframe_for_dtm!(fullData,trnvalcol="trnTest",numcol="ClaimNb",denomcol="Exposure",weightcol="Exposure",independent_vars=selected_explanatory_vars);

updateSettingsMod!(sett,minWeight=-0.03,model_type="boosted_tree",iterations=40,learningRate=0.025,subsampling_features_prop=.7,calculatePoissonError=true)

resultingFiles,resM=dtm(dtmtable,sett) #31s on notebook (x260 with profiling disabled but Revise on!)






##
#consider a boosting model
##
import DecisionTrees: calculateSplitValue,increasing_subsets,build_listOfMeanResponse,_split_feature,some_tree_settings,build_tree!,global_statsperiter_header,createScorebandsUTF8List,initBootstrapSample,_split,sample_data_and_build_tree!,find_max_type,emptyNode,ExcelData,Tree,Leaf,get_feature_pools,build_tree_iteration!,ExcelSheet,Chart

       trnidx=dtmtable.trnidx
    validx=dtmtable.validx
    features=dtmtable.features
    weight=dtmtable.weight
    actualNumerator=dtmtable.numerator
    denominatorDTM=dtmtable.denominator
    mappings=dtmtable.mappings
    candMatWOMaxValues=dtmtable.candMatWOMaxValues
    obstrn,obsval,trn_meanobservedvalue,val_meanobservedvalue,trn_numtot,val_numtot,trn_denomtot,val_denomtot,empty_rows_after_iteration_stats,showTimeUsedByEachIteration,chosen_apply_tree_fn,startAtMean,adaptiveLearningRate,moderationvector,iterations,nameOflistofpredictorsSheet,nameOfpredictorsSheet,nameOfModelStatisticsSheet,nameOfScoresSheet,nameOfOverviewSheet,nameOfSettingsSheet,nameOfValidationSheet,xlData,showProgressBar_time,prroduceEstAndLeafMatrices=DecisionTrees.initSettingsWhichAreTheSameForBoostingAndBagging(trnidx,validx,actualNumerator,denominatorDTM,sett)
    #trn_meanobservedvalue is the mean observed RATIO, i.e. sum(numerator)/sum(denominatorDTM) for the training data
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
    
    if startAtMean
        estimatedRatio=ones(length(weight)).*trn_meanobservedvalue
    else
        error("DTM: This is currently not working: please set startAtMean=true)")
    end
    
    estimatedNumerator=estimatedRatio.*denominatorDTM #these are for instance the estimated losses for a LR model
    
    indicatedRelativityForApplyTree_reused=zeros(obs)
    reused_fitted_leafnr_vector=zeros(Int,obs)
    sortvec_reused_trn_only=zeros(Int,length(trnidx))
    intVarsUsed=Array{Array{Array{Int,1},1}}(undef,0);sizehint!(intVarsUsed,iterations)
    inds_considered=Array{Array{Int,1}}(undef,0);sizehint!(inds_considered,iterations)
        
    res=Array{Union{Leaf,Node{UInt8},Node{UInt16}}}(undef,iterations)
    if prroduceEstAndLeafMatrices
        est_matrix=Array{Float64}(undef,obs,iterations+1)
        est_matrixFromScores=deepcopy(est_matrix)
        MatrixOfLeafNumbers=Array{Int}(undef,obs,iterations+1)
        MatrixOfLeafNumbers[:,1].=0
        est_matrix[:,1]=copy(transpose(estimatedRatio))
        est_matrixFromScores[:,1]=copy(transpose(estimatedRatio))        
    else
        est_matrix=Array{Float64}(undef,0,0)
        est_matrixFromScores=deepcopy(est_matrix)
        MatrixOfLeafNumbers=Array{Int}(undef,0,0)
    end
    
    scores=zeros(Int,obs)    
    estFromScores=zeros(obs)
    #vectorOfRulePathsToLeavesArrays=Array{Array{Array{Rulepath,1},1}}(iterations+1)
    #vectorOfRulePathsToLeavesArrays[1]=Array{Array{Rulepath,1}}(undef,0) #the first entry is not defined    
    vectorOfLeafArrays=Array{Array{Leaf,1}}(undef,iterations+1)
    vectorOfLeafArrays[1]=Array{Leaf}(undef,0) #the first entry is not defined        
    
    p = ProgressMeter.Progress(iterations, 2, "Progress of Boosting Model:") # minimum update interval: x seconds (2)
    ((adaptiveLearningRate>=1.0)||(adaptiveLearningRate<=0.0)) ? (BoolAdaptiveLearningRate=false) : (BoolAdaptiveLearningRate=true)
    #currently not supported
        BoolAdaptiveLearningRate&&error("DTM: this is not yet updated to reflect numerator/denominatorDTM objectives!!")
    #this header needs to be in line with the output of the function createTrnValStatsForThisIteration
    statsPerIteration=global_statsperiter_header    
    currentRelativity=Array{Float64}(undef,obs)
    scoreBandLabels=createScorebandsUTF8List(sett.scorebandsstartingpoints,sett.nScores,addtotal=true)
    #this row is only here such that the variables are initialized outside the for loop
    local fristRowOfThisTable,obsPerScore,vectorWeightPerScore,numPerScore,denomPerScore,scores,scoresVAL,cumulativeStatsPerScoreBand,maxRawRelativityPerScoreSorted,rawObservedRatioPerScore,MAPPINGSmoothedEstimatePerScore

    #the following two lines are for the initialization in case any subsampling (boostrapping) of data is done
    #the idea is to preallocate many of the arrays which are used for sampling (and to save allocation/memory space)
        sampleSizeCanBeNEGATIVE=convert(Int,round(sett.subsampling_prop*obstrn))
        if abs(sampleSizeCanBeNEGATIVE)<1
            sampleSizeCanBeNEGATIVE=convert(Int,sign(sett.subsampling_prop))
        end
        abssampleSize,sampleVector=initBootstrapSample(sampleSizeCanBeNEGATIVE) #,trnidx,validx)    
    
    
    #for iter=1:iterations
        #@timeConditional(showTimeUsedByEachIteration, begin
        #Build Iteration iter
            #showTimeUsedByEachIteration&&println("Calculating iteration $(iter). Time: $(Dates.now())")
            iter=1
            
            current_mdf=moderationvector[min(iter,size(moderationvector,1))] #mdf remains constant if a vector of size 1 was provided (instead of a vector of size iterations)            
            
           T_Uint8_or_UInt16=find_max_type(features)
           
tmpTree=sample_data_and_build_tree!(trnidx,validx,indicatedRelativityForApplyTree_reused,candMatWOMaxValues,mappings,deepcopy(sett),actualNumerator,estimatedNumerator,weight,features,sampleSizeCanBeNEGATIVE,abssampleSize,sampleVector,T_Uint8_or_UInt16)

using BenchmarkTools

@code_warntype sample_data_and_build_tree!(trnidx,validx,indicatedRelativityForApplyTree_reused,candMatWOMaxValues,mappings,deepcopy(sett),actualNumerator,estimatedNumerator,weight,features,sampleSizeCanBeNEGATIVE,abssampleSize,sampleVector,T_Uint8_or_UInt16)


#if no subsampling is done:
nnnode=build_tree!(trnidx,validx,candMatWOMaxValues,mappings,sett,actualNumerator,estimatedNumerator,weight,features,indicatedRelativityForApplyTree_reused,T_Uint8_or_UInt16)


#if no subsampling is done:
@code_warntype build_tree!(trnidx,validx,candMatWOMaxValues,mappings,sett,actualNumerator,estimatedNumerator,weight,features,indicatedRelativityForApplyTree_reused,T_Uint8_or_UInt16)



settings=sett 
fitted_values=indicatedRelativityForApplyTree_reused

intVarsUsed,inds,minweightcalculated=some_tree_settings(trnidx,validx,settings.fixedinds,candMatWOMaxValues,mappings,settings.minWeight,weight,settings.subsampling_features_prop,size(features,2))

	settings.minWeight=minweightcalculated #update minWeight
	if !(length(inds)>0)
        throw(ErrorException(string("Error: no features were selected length(inds)=",length(inds))))
    end
	empty_xl_data=ExcelData(Array{ExcelSheet}(undef,0),Array{Chart}(undef,0))
	fp=get_feature_pools(features)
	resultingTree=Tree(emptyNode,intVarsUsed,candMatWOMaxValues,mappings,inds,settings,empty_xl_data,fp)
	resultingTree.rootnode=build_tree_iteration!(trnidx,validx,settings,resultingTree,actualNumerator,estimatedNumerator,weight,features,0,settings.randomw,Array{Rulepath}(undef,0),Distributed.myid(),fitted_values,T_Uint8_or_UInt16)
    
    @time resultingTree.rootnode=build_tree_iteration!(trnidx,validx,settings,resultingTree,actualNumerator,estimatedNumerator,weight,features,0,settings.randomw,Array{Rulepath}(undef,0),Distributed.myid(),fitted_values,T_Uint8_or_UInt16)
    
    @code_warntype build_tree_iteration!(trnidx,validx,settings,resultingTree,actualNumerator,estimatedNumerator,weight,features,0,settings.randomw,Array{Rulepath}(undef,0),Distributed.myid(),fitted_values,T_Uint8_or_UInt16)
    
    
    ##
    
    
    thisTree=resultingTree
    
    	#!!!! the current concept forsees that features is always the FULL DataFrame
	boolRandomizeOnlySplitAtTopNode=settings.boolRandomizeOnlySplitAtTopNode
	local inds
	#if fixedinds are provided we never change the selection of features for the whole tree! (bagged boosting model)
	if boolRandomizeOnlySplitAtTopNode||(length(settings.fixedinds)>0)
		#these settings are randomly chosen by some_tree_settings in this case
		inds=thisTree.inds_considered
	else
		#inds need to be randomly picked for each iteration
		inds=randomFeatureSelection(size(features,2),settings.subsampling_features_prop)
	end
	intVarsUsed=thisTree.intVarsUsed

	minweight=settings.minWeight
	crit=settings.crit
	spawnsmaller=settings.spawnsmaller
	catSortByThreshold=settings.catSortByThreshold
	catSortBy=settings.catSortBy

	nobs=size(actualNumerator,1)
	fnames=names(features)
	#@code_warntype _split(settings.number_of_num_features,trnidx,validx,numerator,denominator,weight,fnames,features, minweight, depth,randomweight,crit,parallel_level_threshold,parallel_weight_threshold,inds,catSortByThreshold,catSortBy)
	#T<:Union{UInt8,UInt16}=find_max_type(features)<:Union{UInt8,UInt16}
    #T=find_max_type(features)::DataType #Union{UInt8,UInt16}
    
    depth=0
    randomweight=sett.randomw
    
    best_split = _split(one(T_Uint8_or_UInt16),settings.number_of_num_features,trnidx,validx,actualNumerator,estimatedNumerator,weight,fnames,features, minweight, depth,randomweight,crit,inds,catSortByThreshold,catSortBy)
	
    @code_warntype _split(one(T_Uint8_or_UInt16),settings.number_of_num_features,trnidx,validx,actualNumerator,estimatedNumerator,weight,fnames,features, minweight, depth,randomweight,crit,inds,catSortByThreshold,catSortBy)
    
    
    T=T_Uint8_or_UInt16
    val_of_some_UInt_type=one(T)
    number_of_num_features=sett.number_of_num_features
    
    #This function selects the maximal possible split defined by crit (thus depending on the impurity function, we need to put a minus sign in front of it)
		tmpsz=0
		if sum(view(weight,trnidx))<2*minweight;
		#	return const_default_splitdef;
		end
		tmp_splitlist=Vector{Splitdef{T}}(undef,0)
		#for i in inds
        i=inds[1]
			#ATTENTION: for char variables we pass the variable i with a negative sing!!
			#this allows us to distinguish whether we are working on a char or num variable later on
			modified_i = eltype(features[i])<:AbstractString ? number_of_num_features - i  :  i 
			tmplist=_split_feature(val_of_some_UInt_type,number_of_num_features,trnidx,validx,actualNumerator,estimatedNumerator,weight,fnames[i],features[i],minweight,crit,modified_i,randomweight,catSortByThreshold,catSortBy)
            
            
            
            
            
    
@warn("be careful here, features is overwritten")
features=features[i]
    
crit_type=typeof(crit)
#This function is now for numeric and character variables!
#feature_column_id is negative in case of character variables
best_value=-Inf
best_thresh=best_wl=best_wr=NaN
best_subset=Array{UInt8}(undef,0)
trnfeatures=view(features,trnidx)
elt=T#elt=eltype(trnfeatures.parent.refs)
  labellist_sorted=collect(one(elt):convert(elt,length(trnfeatures.parent.pool))) #this used to be levels(features) #this also contains the val feature levels here! It is considerably faster than levels(view) \factor 100 or so
	# THIS IS CRITICAL
	# THE WAY A PooledArray IS CONSTRUCTED, WE WILL ALWAYS HAVE
	# pda.pool[2] == "some string" -> any pda.refs[x].==0x02 will be "some string" 
	#todo/tbd countsort here might be obsolete: we should check if levels is always sorted by construction
  #also we will later sort the labels in a different order anyway (then again the list probably needs to be sorted in the natural manner such that build_listOfMeanResponse is working properly)
  if size(labellist_sorted,1) <= 1
    if T==UInt8
        return UInt8VECTORemptySplitDef
    else
        return UInt16VECTORemptySplitDef #collect(Vector{Splitdef{T}}(undef,0))::Vector{Splitdef{T}}
    end
  else
	#countsort!(labellist_sorted)
    #here I am (somewhat) misusing multiple dispatch since I was too lazy to parametrize this at the moment (todo/tbd in the future)
	if (crit_type==DifferenceSplit||crit_type==PoissonDevianceSplit||crit_type==MaxValueSplit||crit_type==MaxMinusValueSplit)
		labellist,sumnumerator,sumdenominator,sumweight,countlistfloat=build_listOfMeanResponse(crit,trnidx,validx,actualNumerator,estimatedNumerator,weight,trnfeatures,labellist_sorted,minweight)
        
        # @code_warntype build_listOfMeanResponse(crit,trnidx,validx,actualNumerator,estimatedNumerator,weight,trnfeatures,labellist_sorted,minweight)
        
	elseif (crit_type==NormalDevianceSplit)		
		labellist,sumnumerator,sumdenominator,sumweight,countlistfloat,moments_per_pdaclass=build_listOfMeanResponse(crit,numerator,denominator,weight,trnfeatures,labellist_sorted,minweight)
	#else #we catch this possibility earlier when checking the settings
	#	throw(ErrorException(string("Invalid Splitting criterion $(crit)")))
	end
  #todo/tbd
  #here we can introduce the possiblity to sort the labellist (e.g. by meanobserved or median (to be calculated)).
  #then we could only loop through the "increasing list" of sorted labels (instead of doing the 2^ncategories exhaustive search (bitflip_graycode_subsets))
  if feature_column_id>0 #id>0 -> we are working on a numeric column
		#only consider to split at the candidate split points
		subs=increasing_subsets(labellist)
  else #id<0 -> we are working on a character column
		#distinguish between exhaustive and "increasing" search for split point
	if size(labellist_sorted,1)>catSortByThreshold
		if (crit_type==DifferenceSplit||crit_type==MaxValueSplit||crit_type==MaxMinusValueSplit)
			sortlists!(catSortBy,labellist,sumnumerator,sumdenominator,sumweight,countlistfloat) #catSortBy::SortBy=SORTBYMEAN
		elseif (crit_type==NormalDevianceSplit)
			sortlists!(catSortBy,labellist,sumnumerator,sumdenominator,sumweight,countlistfloat,moments_per_pdaclass)
		end
		subs=increasing_subsets(labellist)
	else
	#perform exhaustive search
		subs=bitflip_graycode_subsetsHALF(labellist)
	end
	end

    
    
    fname=:hopefully_not_important
    
	if randomweight==0.0
	if (crit_type==DifferenceSplit||crit_type==MaxValueSplit||crit_type==MaxMinusValueSplit)
		tmp_result=calculateSplitValue(crit,fname,number_of_num_features,labellist,sumnumerator,sumdenominator,sumweight,countlistfloat,minweight,subs)