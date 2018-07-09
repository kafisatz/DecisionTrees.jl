
import Base: length,start,next,done,eltype,==,hash

export updateSettings!,updateSettingsMod!
export Splitdef,Rulepath,Leaf,Node,Tree,BoostedTree,SplittingCriterion,CVOptions,BaggedTree
export ModelSettings,copySettingsToCurrentType

#struct pdaMod end #this is a legacy type. It is only defined here, because we have some 'old' unused functions in the code which have a singature with this type (and the signature needs an update....)
    
#Splitting Criterions use multiple dispatch
	abstract type SplittingCriterion  end
	struct DifferenceSplit      <: SplittingCriterion end
	struct NormalDevianceSplit      <: SplittingCriterion end
	struct PoissonDevianceSplit      <: SplittingCriterion end
	struct GammaDevianceSplit      <: SplittingCriterion end
	struct MaxValueSplit      <: SplittingCriterion end
	struct MaxAbsValueSplit      <: SplittingCriterion end
	struct MaxMinusValueSplit      <: SplittingCriterion end
	struct MSESplit             <: SplittingCriterion end
	struct RankOptSplit         <: SplittingCriterion end
	struct ROptMinRLostPctSplit	   <: SplittingCriterion end
	struct ROptMinRLostSplit	   <: SplittingCriterion end
    
    PoissonOrGamma = Union{PoissonDevianceSplit,GammaDevianceSplit}
    
#sortby options for categorical splits
	abstract type SortBy end
	struct SortByMean 	<: SortBy end

function canbeint(a)
	try convert(Int,parse(Float64,a))
	   return true
	catch
		return false
	end
end

function canbefloat(a)
	try Float64(a)
	   return true
	catch
		return false
	end
end

mutable struct Chart	
	sheet::String
	chartDict::Dict{AbstractString,Dict{AbstractString,Any}}
	location::String
end

mutable struct ExcelSheet
	name::String #I think this name is currently not used anywhere
	data::DataFrame
end

mutable struct ExcelData
	sheets::Array{ExcelSheet,1}
    charts::Array{Chart,1}
    function ExcelData()
        return new(Array{ExcelSheet}(undef,0),Array{Chart}(undef,0))
    end
    function ExcelData(a,b)
        return  new(a,b)
    end
end

#errorStats
struct ErrorStats
	mse::Float64
	mae::Float64
	mrse::Float64
	mrae::Float64
	pearsoncorrelation::Float64
	function ErrorStats()
		new(NaN,NaN,NaN,NaN,NaN)
	end
	function ErrorStats(a::Float64,b::Float64,c::Float64,d::Float64,e::Float64)
		new(a,b,c,d,e)
	end
end

#subsets which are considered by splits
abstract type DTSubsets end

#custom loop over subsets (Gray Code) ONLY UNTIL the middle - USING SYMMETRY of the list as we consider binary splits
	struct MyGrayCodeSubsetsHALF <: DTSubsets
		size_of_set::Int
		stopat::Int #(1 << it.size_of_set)-2
	end

	eltype(it::MyGrayCodeSubsetsHALF) = Int #Array{eltype(1),1}
	length(it::MyGrayCodeSubsetsHALF) = it.stopat+1 #1 << it.size_of_set
	start(it::MyGrayCodeSubsetsHALF) = zero(1) #state = simply count through the subsets
	done(it::MyGrayCodeSubsetsHALF,state) = state>it.stopat  #if we have equality here, then the full set will be shifted / this is not needed for a binary split, hence we stop beforehand

	function next(it::MyGrayCodeSubsetsHALF, state)
	#returns the element which needs to be flipped
	  bit_to_flip=1+trailing_ones(state) #in the function print_graycode bit_to_flip is defined without the plus 1
	  state+=one(1)
	  bit_to_flip,state
	end
    @inline iterate(x::MyGrayCodeSubsetsHALF) = next(x, start(x))
    @inline iterate(x::MyGrayCodeSubsetsHALF, i) = done(x, i) ? done : next(x, i)

	bitflip_graycode_subsetsHALF(xs) = MyGrayCodeSubsetsHALF(length(xs),(1 << (length(xs)-1))-2)

#normal subset list for numeric splitting points ("increasing" subset)
	struct MyIncreasingSubsets <: DTSubsets
		size_of_set::Int
	end

	eltype(it::MyIncreasingSubsets) = Array{eltype(1),1}
	length(it::MyIncreasingSubsets)=it.size_of_set
	start(it::MyIncreasingSubsets)=zero(1) #state = simply count through the subsets
	done(it::MyIncreasingSubsets,state)=state>=it.size_of_set-1 #we do not need to flip/change the last element (as this would resutl in an empty left child)
    @inline iterate(x::MyIncreasingSubsets) = next(x, start(x))
    @inline iterate(x::MyIncreasingSubsets, i) = done(x, i) ? done : next(x, i)

	function next(it::MyIncreasingSubsets,state)
		state+=one(1)
		return copy(state),state
	end
	increasing_subsets(xs)=MyIncreasingSubsets(length(xs))

#in constrast to types, immutables cannot (or should not?) be changed over time
#however immutables should increase performance (todo: need to check this or read more about it)
mutable struct DTMTable
	key::Vector{String}
	trnidx::Vector{Int}
	validx::Vector{Int}
	numerator::Vector{Float64}
	denominator::Vector{Float64}
	weight::Vector{Float64}
	features::DataFrame
	candMatWOMaxValues::Vector{Vector{Float64}}
	mappings::Vector{Vector{String}}
	function DTMTable(key,trnidx,validx,numerator,denominator,weight,features,cands,maps)
		@assert length(validx)+length(trnidx)==length(key)
		@assert length(key)==length(numerator)
		@assert length(key)==length(denominator)
		@assert length(key)==length(weight)
		@assert length(key)==size(features,1)
		@assert length(cands)+length(maps)==size(features,2)
		@assert length(findall((in)(validx),trnidx))==0 "DTM: trnidx and validx are not mutually exclusive"
		return new(key,trnidx,validx,numerator,denominator,weight,features,cands,maps)
	end
end

Base.size(d::DTMTable) = size(d.features)

function Base.getindex(r::DTMTable,idx,compact_features=false)
    @warn("DTM: getindex(x::DTMTable,idx) was called. This function is highly experimental and untested!")    
    
    key=r.key[idx]
    one_to_n=collect(1:length(r.weight))
    rows_to_keep=one_to_n[idx]        
    new_trnidx=indexin(r.trnidx,rows_to_keep)
    new_validx=indexin(r.validx,rows_to_keep)
    filter!(!iszero,new_trnidx)
    filter!(!iszero,new_validx)
    numerator=r.numerator[idx]    
    denominator=r.denominator[idx]    
    weight=r.weight[idx]    
	
    features=r.features[idx,:]
    if compact_features 
        tmp_number_of_num_features=sum(map(x->eltype(r.features[x].pool),1:size(r.features,2)).!=String)	        
        #"compact" features, it is possible that we could do this in a better/faster way
        for i=tmp_number_of_num_features+1:size(features,2)
            if eltype(features[i])!=UInt8
                features[i]=DecisionTrees.PooledArraysDTM.PooledArray(Array(features[i]))
            end
        end    
	end 
    
    mp=deepcopy(r.mappings)
    @warn("DTM: getindex(x::DTMTable,idx) was called. This function is highly experimental and untested!")
    dt=DTMTable(key,new_trnidx,new_validx,numerator,denominator,weight,features,deepcopy(r.candMatWOMaxValues),mp)
    return dt
end

#todo check the difference between immutable and type!! are we using immutables correctly here?
struct Splitdef{T<:Unsigned} 
	featid::Int
	featid_new_positive::Int
    featurename::Symbol
    subset::Array{T,1}
    splitvalue::Float64 #Depends on the criterion function used (e.g. difference)
    weightl::Float64
    weightr::Float64
end
==(x::Splitdef,y::Splitdef)= (x.featid==y.featid) && (x.featurename==y.featurename) && (x.featid_new_positive==y.featid_new_positive) && (x.splitvalue==y.splitvalue) && (x.subset==y.subset)  && (x.weightl==y.weightl)  && (x.weightr==y.weightr)
#Splitdef(x::Tuple{Int,Float64,Array{UInt8,1},Float64,Float64,Float64})= Splitdef(x[1],x[2],x[3],x[4],x[5],x[6])
hash(x::Splitdef)=hash(x,UInt(9))
function hash(x::Splitdef,h::UInt)
  h=hash(x.featid,h)
  h=hash(x.featid_new_positive,h)
  h=hash(x.featurename,h)
  h=hash(x.subset,h)
  h=hash(x.splitvalue,h)
  h=hash(x.weightl,h)
  h=hash(x.weightr,h)
  return h
end

struct Rulepath{T<:Unsigned}
    featid::Int
    subset::Vector{T}
    isLeftChild::Bool
end

Rulepath{T}() where {T<:Unsigned} = Rulepath(0,T[],false)

==(x::Rulepath,y::Rulepath)= (x.featid==y.featid) && (x.isLeftChild==y.isLeftChild) && (x.subset==y.subset)
hash(x::Rulepath)=hash(x,UInt(9))
function hash(x::Rulepath,h::UInt)
  h=hash(x.featid,h)
  h=hash(x.isLeftChild,h)
  for i in x.subset
    h = hash(i, h)
  end
  return h
end

mutable struct Leaf{T<:Unsigned} #leafs are mutable, since we determine the number of the leaf a posteriori
    #NOTE: the idx array can have an unexpected meaning/size in the case where we perform subsampling!
	rowcount::Int #number of observations of the training data set that fell into thta leaf
    mean::Float64
    fitted::Float64 #this is usually identical to the mean, but can be different for some instances
  	size::Float64 #  this is equal to sum(weight)
    depth::Int
    rule_path::Array{Rulepath{T},1}
	sumnumerator::Float64
	sumdenominator::Float64
    id::Int #a number 1<=id<=#number of leaves in the tree
end
function Leaf{T}() where T
        rp=Vector{Rulepath{T}}(undef,0)
        return Leaf{T}(0,0.0,0.0,0.0,0,rp,0.0,0.0,0)
    end
function Leaf(rowc,me,fi,sz,depth,rp::Vector{Rulepath{T}},sumn,sumd,id) where T
    return Leaf{T}(rowc,me,fi,sz,depth,rp,sumn,sumd,id)
end
==(x::Leaf,y::Leaf)= (x.rowcount==y.rowcount) && (x.mean==y.mean) && (x.fitted==y.fitted)  && (x.size==y.size)  && (x.depth==y.depth)  && (x.rule_path==y.rule_path)  && (x.sumnumerator==y.sumnumerator)  && (x.sumdenominator==y.sumdenominator) && (x.id==y.id)
hash(x::Leaf)=hash(x,UInt(9))
function hash(x::Leaf,h::UInt)
  h=hash(x.rowcount,h)
  h=hash(x.mean,h)
  h=hash(x.fitted,h)
  h=hash(x.size,h)
  h=hash(x.depth,h)
  h=hash(x.rule_path,h)
  h=hash(x.sumnumerator,h)
  h=hash(x.sumdenominator,h)
  h=hash(x.id,h)
  return h
end

struct Node{T<:Unsigned}
    #todo, add depth!
    featid::Int
    featid_new_positive::Int
    #we have redundance in the tree structure, the rulepath contains the same information as the members of the node; todo - remove it!
    subset::Array{T,1}
    left::Union{Leaf{T},Node{T}} #,Node{UInt16}}
    right::Union{Leaf{T},Node{T}} #,Node{UInt16}}
    rule_path::Array{Rulepath{T},1}
end

mutable struct ModelSettings
	model_type::String #1
	minWeight::Float64 #2
	randomw::Float64 #3
	crit::SplittingCriterion # fn version of 4
	maxSplittingPoints::Int #5
	iterations::Int #6
	learningRate::Float64 #7 moderationfactor
	nScores::Int #8
	adaptiveLearningRate::Float64 #9
	prem_buffer::Int #11
	startAtMean::Bool #12
	writeTree::Bool #13	
	number_of_num_features::Int #14
	spawnsmaller::Bool #26
    boolRankOptimization::Bool #35  
	boolRandomizeOnlySplitAtTopNode::Bool #37
	subsampling_prop::Float64 #38 #NOTE Friedman suggests to sample without replacement see Friedman 2002 Stochastic Gradient Boosting,” Computational Statistics and Data Analysis 38(4):367-378
	subsampling_features_prop::Float64 #39
	version::String #40
	preppedJLDFileExists::Bool #41
	catSortByThreshold::Int #42
	catSortBy::SortBy #43
	scorebandsstartingpoints::Array{Int,1} #44
	showTimeUsedByEachIteration::Bool #45
	smoothEstimates::String #46
	deriveFitPerScoreFromObservedRatios::Bool #deriveFitPerScoreFromObservedRatios::Bool (if this is true then the fit per score will be based on the observed ratio per score (instead of the fitted ratio per score)) . Note that this does not influence the structure of the tree, but only the fitted ratio per score (and scoreband)
    roptForcedPremIncr::Bool #47
    writeSasCode::Bool #49
	writeIterationMatrix::Bool #50
	writeResult::Bool #51
	writeStatistics::Bool #52
	boolCreateZipFile::Bool #53
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
	saveJLDFile::Bool #if this is set to false no *.jld2 file is saved #todo we should rename this to ensure it is clear that this refers to the prepped data only!
	saveResultAsJLDFile::Bool #this refers to the model output/result
	print_details::Bool #whether julia should print details to the log or not (note to developers: warnings/errors and certain important infos should always be printed!)
	seed::Int #optional entropy
	graphvizexecutable::String #e.g. C:\Program Files (x86)\Graphviz2.38\ if this
	showProgressBar_time::Bool # bar form iterators package (or some other package)
	prroduceEstAndLeafMatrices::Bool
	write_dot_graph::Bool
	calculateGini::Bool #whether or not to calculate the gini (which requires sorting of the data (which takes time!))
	calculatePoissonError::Bool
    performanceMeasure::String
	fitForStatsAndCharts::String
	ignoreZeroDenominatorValues::Bool

	#the following are treated specially
	df_name_vector::Array{String,1}
    number_of_char_features::Int #this is calculated by Julia from number_of_num_features and ncolsdfIndata (which is defined by the data)

	chosen_apply_tree_fn::String #defined below
	moderationvector::Array{Float64,1} #50 onwards

  function ModelSettings()
	#initialize default settings
	model_type="build_tree" #1
	minWeight=-0.1 #2
	randomw=0.0 #3
	crit=DifferenceSplit()	 #4
	maxSplittingPoints=250 #5 #arbitrary choice. Preliminary testing has shown no performance difference for different values (say 10 or 200) for this setting
	iterations=2 #6
	learningRate=0.1 #7
	nScores=1000 #8
	adaptiveLearningRate=1.0 #9
	prem_buffer=0 #11
	startAtMean=true #12
	writeTree=true #13
	number_of_num_features=-1 #14    
	spawnsmaller=true #26
	boolRankOptimization=false #35
	boolRandomizeOnlySplitAtTopNode=true #37
	subsampling_prop=1.0 #38
	subsampling_features_prop=1.0 #39
	version="not_initialized" #40
	preppedJLDFileExists=false #41
	catSortByThreshold=7 #42
	catSortBy=SortByMean() #43
	scorebandsstartingpoints=[parse(Int,x) for x in split("1,100,200,300,400,500,600,700,800,900",',')] #scorebandsstartingpoints=parse(Int,split("1,100,200,300,400,500,600,700,800,900",',')) #44
	showTimeUsedByEachIteration=false #true #45
	smoothEstimates="1" #46 this is a string for now (as there could be different smoothing methods specified by this string)
	deriveFitPerScoreFromObservedRatios=true
	roptForcedPremIncr=false #47
    writeSasCode=false #49
	writeIterationMatrix=false#50
	writeResult=false #true #51
	writeStatistics=true #52
	boolCreateZipFile=false #true #53
	writeCsharpCode=false #true
	writeVbaCode=false
	nDepthToStartParallelization=-1
	baggingWeightTreesError="mae" # ::String
	#combined bagging&boosting model
	cBB_niterBoosting=0
	cBB_niterBagging=0
	fixedinds=Array{Int}(undef,0)
	boolNumeratorStats=false	
	statsByVariables=Int[]
	statsRandomByVariable=5
	saveJLDFile=true
	saveResultAsJLDFile=false
	print_details=true
	seed=90210
	graphvizexecutable="" #C:\\Program Files (x86)\\Graphviz2.38\\bin\\"
	showProgressBar_time=true
	prroduceEstAndLeafMatrices=false
	write_dot_graph=false
	calculateGini=false
	calculatePoissonError=false
    performanceMeasure="Lift Val"
	fitForStatsAndCharts="rawRelativities" #rawRelativities,unsmoothedPerScore,smoothedPerScore
	ignoreZeroDenominatorValues=false

	#the following are treated specially
	df_name_vector=Array{String}(undef,0) #
    number_of_char_features=-1 #
	chosen_apply_tree_fn="apply_tree_by_leaf" # #apply_tree_by_row does not seem to work for (certain?) boosting models
	moderationvector=[0.1] #

	return new(model_type,minWeight,randomw,crit,maxSplittingPoints,iterations,learningRate,nScores,adaptiveLearningRate,prem_buffer,startAtMean,writeTree,number_of_num_features,spawnsmaller,boolRankOptimization,boolRandomizeOnlySplitAtTopNode,subsampling_prop,subsampling_features_prop,version,preppedJLDFileExists,catSortByThreshold,catSortBy,scorebandsstartingpoints,showTimeUsedByEachIteration,smoothEstimates,deriveFitPerScoreFromObservedRatios,roptForcedPremIncr,writeSasCode,writeIterationMatrix,writeResult,writeStatistics,boolCreateZipFile,writeCsharpCode,writeVbaCode,nDepthToStartParallelization,baggingWeightTreesError,cBB_niterBoosting,cBB_niterBagging,fixedinds,boolNumeratorStats,statsByVariables,statsRandomByVariable,saveJLDFile,saveResultAsJLDFile,print_details,seed,graphvizexecutable,showProgressBar_time,prroduceEstAndLeafMatrices,write_dot_graph,calculateGini,calculatePoissonError,performanceMeasure,fitForStatsAndCharts,ignoreZeroDenominatorValues
	 ,df_name_vector,number_of_char_features,chosen_apply_tree_fn,moderationvector)
  end  # ModelSettings()

end #type definition

function edit_statsByVariables(col,stringValue,s::ModelSettings)
if col=="statsbyvariables"&&length(stringValue)>0
	if in(',',stringValue)||in(';',stringValue)
		println("statsbyvariables=$(stringValue)")
		error("DTM: statsbyvariables must be space separated not comma separated!")
	end
	try
	chosenvars=[lowercase(x) for x in split(stringValue,' ')]
	integerindices=findall(in(chosenvars), [lowercase(x) for x in s.df_name_vector])
	if length(integerindices)!=length(chosenvars)
		@warn("DMT: Some variables of statsbyvariables were not found in the data")
		@show s.df_name_vector
		@show chosenvars
		@show stringValue
		@show integerindices
		@show [lowercase(x) for x in s.df_name_vector]
		error("DTM: Abort")
		else
		integerindices2=[x<=s.number_of_num_features ? x : -1*(x-s.number_of_num_features) for x in integerindices]
		return string(integerindices2)[2:end-1] #remove leading and trailing '[' and ']'
	end
	catch er
		@show er
		error("DTM: Unable to parse statsbyvariables option")
	end
end
return stringValue
end

function copySettingsToCurrentType(oldSetting)
	s=ModelSettings()
	newfields=fieldnames(typeof(s))
	oldfields=fieldnames(typeof(oldSetting))    
	
	for field in fieldnames(typeof(oldSetting))
		if in(field,newfields)
			v=getfield(oldSetting,field)
			setfield!(s,field,v)
		end
    end
	
	#consider new and dropped fields
		newfs=setdiff(newfields,oldfields)
		droppedfs=setdiff(oldfields,newfields)
	@info "DTM: Copied settings to new type."
	if length(newfs)>0
		@info "DTM: New type has the following new fields:"
		for f in newfs        
			v=getfield(s,f)
			@show f,v
		end
	end
	if length(droppedfs)>0
		@info "DTM: Old type had the following fields which were dropped:"
		for f in droppedfs        
			v=getfield(oldSetting,f)
			@show f,v
		end
	end
    return s
end

"""
updateSettings!(s::ModelSettings;args...)
use this function to update the model settings 
try for instance:

updateSettings!(sett,minWeight=-0.3,iterations=30)
"""
function updateSettings!(s::ModelSettings;args...)
	updateSettingsMod!(s;args...)
end

"""
updateSettingsMod!(s::ModelSettings;args...)
use this function to update the model settings 
try for instance:

updateSettingsMod!(sett,minWeight=-0.3,iterations=30)
"""
function updateSettingsMod!(s::ModelSettings;args...)
	seen=[]
	for (smember,v) in args
		if in(smember,seen)
			error("DTM: You provided multiple values for the field $(smember)")
		end
		push!(seen,smember)

		oldValue=getfield(s,smember)		
		if typeof(v)!=typeof(oldValue)
			valConverted=convertFromString(oldValue,string(v)) #tbd maybe string(v) could/should be replaced with v here
			setfield!(s,smember,valConverted)
		else
			setfield!(s,smember,v)
		end		
	end
	res=checkIfSettingsAreValid(s)
	@assert res==nothing "Some settings are invalid! Abort."
end

myParseAsFloat(val::T) where {T<:AbstractString} = parse(Float64,val)
myParseAsFloat(val::T) where {T<:Number} = float(val)
convertFromString(oldvalue::T,val) where {T <: Any}=convert(T,val) #generic method catchall
convertFromString(oldvalue::T,val) where {T <: AbstractString}=convert(T,string(val))

#convertFromString(oldvalue::T,val::U) where {T <: AbstractFloat,U<:AbstractString}=parse(Float64,val)
convertFromString(oldvalue::T,val) where {T <: AbstractFloat}=myParseAsFloat(val)
convertFromString(oldvalue::T,val) where {T <: Integer}=convert(Integer,myParseAsFloat(val))
#convertFromString(oldvalue::T,val::U) where {T <: Integer,U<:AbstractString}=convert(Integer,parse(Float64,val)) 
#convertFromString(oldvalue::T,val::U) where {T <: Unsigned,U<:AbstractString}=uint(parse(Float64,val))
convertFromString(oldvalue::T,val) where {T <: Unsigned}=uint(myParseAsFloat(val))
convertFromString(oldvalue::T,val) where {T <: Bool}=(lowercase(val)=="t"||lowercase(val)=="true"||val=="1")
function convertFromString(oldvalue::T,val) where {T <: SortBy}
	@assert lowercase(val)=="sortbymean"||lowercase(val)=="mean"
	return SortByMean()  #currently only mean is possible here
end
function convertFromString(oldvalue::T,critfnstr) where {T <: SplittingCriterion}
	local allowedinputstrings=[lowercase(x) for x in ["difference" "normal" "poisson" "gamma" "maxvalue" "maxabsvalue" "maxminusvalue"]]
	critfnstr=lowercase(critfnstr)
	if (critfnstr=="difference"||critfnstr=="default")
		crit=DifferenceSplit()
	elseif (critfnstr=="mse"||critfnstr=="normaldeviance" ||critfnstr=="normal" ||critfnstr=="gaussian" || critfnstr=="gaussiandeviance" || critfnstr=="rss")
		crit=NormalDevianceSplit()
	elseif (critfnstr=="poissondeviance"||critfnstr=="poisson")
		crit=PoissonDevianceSplit()
	elseif (critfnstr=="gammadeviance"||critfnstr=="gamma")
		crit=GammaDevianceSplit()
	elseif critfnstr=="maxvalue" #pick the split which leads to the maximum average value (in either the right or left child)
		crit=MaxValueSplit()
	elseif critfnstr=="maxabsvalue" #pick the split which leads to the absolute maximum average value (in either the right or left child)
		crit=MaxAbsValueSplit()
	elseif critfnstr=="maxminusvalue" #pick the split which leads to the maximum negative value (in either the right or left child)
		crit=MaxMinusValueSplit()
	elseif critfnstr=="roptminrlostpct"
		crit=ROptMinRLostPctSplit()
	elseif critfnstr=="roptminrlost"
		crit=ROptMinRLostSplit()
	else
		@info "Invalid splitting criterion: $(critfnstr). Currently only the following strings are allowed as input:"
		@info string(allowedinputstrings)
		error("DTM: Abort due to invalid splitting criterion (see above).")
	end
	return crit
end

function convertFromString(oldvalue::Array{Int,1},val)
local arr
    try
		arr=Int[parse(Int,x) for x in split(val,',')]
	catch er
		@show er
		error("DTM: Unable to parse input string as integer vector: input= $(val)")
	end
return arr
end

function convertFromString(oldvalue::Array{Float64,1},val)
local arr
    try
		arr=Float64[parse(Float64,x) for x in split(val,',')]
	catch er
		@show er
		error("DTM: Unable to parse input string as integer vector: input= $(val)")
	end
return arr
end

function checkForMandatorySettings(headerLowercase)
	@assert in("number_of_num_features",headerLowercase) global_number_of_num_f_warning_mandatory_field
	return nothing
end

function checkIfSettingsAreValid(s::ModelSettings)
	  #@assert in(method,globalConstAllowableMethodsForDefineCandidates)
	  @assert in(s.fitForStatsAndCharts,globalValidfitForStatsAndCharts) 	  
      @assert s.statsRandomByVariable<typemax(UInt8) #more than 255 random groups is really not meaningful; we do work with a UInt8 list later on
	  @assert in(s.model_type,["boosted_tree" "build_tree" "bagged_tree" "bagged_boosted_tree"])
      @assert (s.randomw>=0) & (s.randomw<=1)
	  #Note, we could easily support more splitting Points (i.e. UInt), but it is likely that the user has not intended to exceed typemax(UInt16), thus we catch this here
	  max_splitting_points_num_INTERNAL=Int(typemax(UInt16))-1
      @assert s.maxSplittingPoints<max_splitting_points_num_INTERNAL "Maximum number of splitting points is currently limited to $(max_splitting_points_num_INTERNAL)"
	  @assert s.maxSplittingPoints>0
	  if (s.model_type=="bagged_boosted_tree")
		@assert min(cBB_niterBagging,cBB_niterBoosting)>0
		@assert abs(subsampling_prop)<1.0 #without subsampling bagged boosting probably does not make too much sense
		@assert subsampling_features_prop<1.0 #without sampling features bagged boosting probably does not make too much sense
		@assert s.iterations==0 "Please set iterations to 0 for Combined Bagging&Boosting Models" #I want to avoid confusion where the use thinkgs that this parameter is of relevance.
	  else
		@assert s.iterations>0
	  end
	  @assert in(s.baggingWeightTreesError,["mse" "pearsoncorrelation" "mrae" "mrse" "mae" "uniform"])	  
	  #if length(s.moderationvector)=1
	  if in(s.model_type,["build_tree"])
		s.moderationvector=[s.learningRate]
	  end
	  if length(s.moderationvector)==1
		if s.moderationvector[1]!=s.learningRate && s.model_type=="boosted_tree"
			@info "Replacing moderationvector[1] with the value of learningRate=$(s.learningRate)"
			s.moderationvector[1]=s.learningRate
		end
	  end
	  if s.model_type=="boosted_tree"
		if !((length(s.moderationvector)==1)&&(s.moderationvector[1]==s.learningRate))
			@show s.moderationvector
			@show s.learningRate
			error("DTM: Invalid learningRate or moderationvector settings")
	  	end
		end
      @assert (s.learningRate>=0 && s.learningRate<=1)
	  #end
      @assert (s.adaptiveLearningRate<=1.0 && s.adaptiveLearningRate>=0.0)
      #@assert s.subsampling_prop<=1.0 #positive value => without replacement
	  @assert s.subsampling_prop>=-1.0 #negative value => with replacement, smaller than -1 is not possible wherease larger than 1 is possible
	  if s.subsampling_prop<1.0
		@warn("Subsampling of the data is enabled. Ensure the value of subsampling_prop=$(s.subsampling_prop) has a meaningful value relative to minWeight=$(s.minWeight)")
	  end
	  #if (s.model_type!="bagged_tree")&&(abs(s.subsampling_prop)!=1.0);@warn("Subsampling with Boosting methods is experimental!");end
	  #todo/tbd need to fix this! see goodnessOfFit
		if (s.model_type=="bagged_tree")&&(abs(s.subsampling_prop)!=1.0);@warn("Bagging is currently EXPERIMENTAL!") ;end
      #@assert (abs(s.subsampling_prop)==1.0)||(s.model_type=="bagged_tree") "ERROR: (bk) Sumbsampling is currently only available for Bagging, as it will require the usage of apply tree after each iteration (since not all estimates are available after the construction of the tree)! Thus (for now) you need to set subsampling_prop to 1.0."
      @assert (abs(s.subsampling_features_prop>0) && abs(s.subsampling_features_prop<=1))
	  @assert s.scorebandsstartingpoints[1]==1
	  @assert issorted(s.scorebandsstartingpoints)
	  @assert s.scorebandsstartingpoints[end]<s.nScores
	  @assert s.nScores>8 "Please set nScores to a value greater than 8" #Less than 8 scores makes it hard/impossible for the moving average
	  @assert in(s.smoothEstimates,["0" "1"]) #this is currently a string since there could (in the future) be multiple smoothing methods
	  if !(length(s.graphvizexecutable)<1||isfile(s.graphvizexecutable))
		@show s.graphvizexecutable
		error("DTM: Invalid graphvizexecutable: $(graphvizexecutable) is not a file.")
	  end
	if s.roptForcedPremIncr
		@assert (s.crit==ROptMinRLostPctSplit()||s.crit==ROptMinRLostSplit())
	end
	@assert !s.boolRankOptimization "this is currently not working in this version of the Code"
	#@assert s.ishift==0 #is zero for normal models (historically this was 1 for rkoptmodels (labels_orig and so on....)
	#@show s.df_name_vector
	#@show s.number_of_char_features,s.number_of_num_features
	if length(s.df_name_vector)!=(s.number_of_char_features+s.number_of_num_features)
		@show length(s.df_name_vector),s.number_of_char_features,s.number_of_num_features
		@assert length(s.df_name_vector)==(s.number_of_char_features+s.number_of_num_features)
	end
    @assert in(s.performanceMeasure,global_statsperiter_header)

     return nothing
end

abstract type DTModel end
mutable struct Tree <: DTModel
	rootnode::Union{Leaf,Node{UInt8},Node{UInt16}}
	intVarsUsed::Array{Array{Int,1},1}	
	candMatWOMaxValues::Array{Array{Float64,1},1}
	charMappings::Array{Array{String,1},1}
	inds_considered::Array{Int,1}	
	settings::ModelSettings
	variableImp1Dim::Array{Int,1}
	variableImp2Dim::Array{Int,1}	
	modelstats::DataFrame # #overall model statistics and performance (trn and val), this needs to be a matrix (Any) with a header row NOTE: this is a feature of multiple model types! (tree, boosting, bagging) and we need consistency because of the function run_model_multirow_settings
	exceldata::ExcelData
	featurepools::Array{Union{Array{Float64,1},Array{String,1}},1}
	function Tree(rootnode::Union{Leaf,Node{UInt8},Node{UInt16}},intVarsUsed::Array{Array{Int,1},1},candMatWOMaxValues::Array{Array{Float64,1},1},charMappings::Array{Array{String,1},1},inds_considered::Array{Int,1},settings::ModelSettings,exceldata::ExcelData,fp::Array{Union{Array{Float64,1},Array{String,1}},1})
		nf=settings.number_of_char_features+settings.number_of_num_features
		@assert length(settings.df_name_vector)==nf
		variableImp1Dim=zeros(Int,nf)
		variableImp2Dim=zeros(Int,div(nf*(nf-1),2))
		ms=Array{Any}(undef,0,0)
		return new(rootnode,intVarsUsed,candMatWOMaxValues,charMappings,inds_considered,settings,variableImp1Dim,variableImp2Dim,ms,exceldata,fp)
	end
	function Tree()		
        rootnode=Leaf{UInt8}()
        intVarsUsed=Int[]
        candMatWOMaxValues=[Float64[]]
        charMappings=[String[]]
        inds_considered=Int[]
        settings=ModelSettings()
        variableImp1Dim=Int[]
        variableImp2Dim=Int[]
        ms=DataFrame()
        exceldata=ExcelData()
        fp=[Float64[]]
		return new(rootnode,intVarsUsed,candMatWOMaxValues,charMappings,inds_considered,settings,variableImp1Dim,variableImp2Dim,ms,exceldata,fp)
	end
end

#this type was only created, such that the @distributed results can be created more easilty (with a vcat operation)
#there would have been other workarounds which would have worked too
mutable struct TreeWithErrorStats
	tree::Tree
	err::ErrorStats
end

abstract type Ensemble <: DTModel end

struct EmtpyDTModel <:DTModel
    modelstats::DataFrame #this might not be needed  
    function EmtpyDTModel()
        return new(DataFrame())
    end
end

#This is currently only used if multi run fails (as a return value of an 'empty' model)
struct EmptyEnsemble <: Ensemble      
    modelstats::DataFrame #this might not be needed  
    function EmptyEnsemble()
        return new(DataFrame())
    end
end

struct BoostedTree <: Ensemble
    trees::Vector{Union{Leaf,Node{UInt8},Node{UInt16}}} #todo tbd make this to trees!

	settings::ModelSettings
	intVarsUsed::Array{Array{Array{Int,1},1},1}
	candMatWOMaxValues::Array{Array{Float64,1},1}
	charMappings::Array{Array{String,1},1}
	inds_considered::Array{Array{Int,1},1}
	
    moderationvector::Array{Float64,1}
    scores::Array{Int,1}
    rawrelativities::Array{Float64,1}
    maxRawRelativityPerScoreSorted::Array{Float64,1}
    meanobserved::Float64
    startAtMean::Bool
	ScoreToSmoothedEstimate::Array{Float64,1}
	rawObservedRatioPerScore::Array{Float64,1}
	iterationmatrix::Array{Float64,2}
	modelstats::DataFrame #overall model statistics and performance (trn and val), this needs to be a matrix (Any) with a header row NOTE: this is a feature of multiple model types! (tree, boosting, bagging) and we need consistency because of the function run_model_multirow_settings
	exceldata::ExcelData
	trnidx_one_zero_full::Array{UInt,1}
	featurepools::Array{Union{Array{Float64,1},Array{String,1}},1}
end

struct BaggedTree <: Ensemble
    trees::Vector{Node} #todo tbd make this to trees! (? really)
	weightPerTree::Array{Float64,1}
	oobagErr::Array{ErrorStats,1}

	settings::ModelSettings
	intVarsUsed::Array{Array{Array{Int,1},1},1}
	candMatWOMaxValues::Array{Array{Float64,1},1}
	charMappings::Array{Array{String,1},1}
	inds_considered::Array{Array{Int,1},1}
	
    scores::Array{Int,1}
    rawrelativities::Array{Float64,1}
    maxRawRelativityPerScoreSorted::Array{Float64,1}
    meanobserved::Float64
    startAtMean::Bool
	ScoreToSmoothedEstimate::Array{Float64,1}
	rawObservedRatioPerScore::Array{Float64,1}
	iterationmatrix::Array{Float64,2}
	modelstats::DataFrame
	trnidx_one_zero_full::Array{UInt,1}
	featurepools::Array{Union{Array{Float64,1},Array{String,1}},1}
end


"""
cvo=CVOptions(folds, training_proportion, use_all_data)

folds is an integer. If folds<0 then we consider n disjoint training sets

use_all_data::Bool #if set to false, the CVOptions will never use the 'original' validation rows in the data

training_proportion::Float64 #size of trn data. The value must be >0 and <1

More precisely the evaluation is as follows

if cvo.training_proportion>0

\t   training_proportion > 0 provided. Performing random sampling without replacement 'cvo.folds' samples)

else   

\t if cvo.folds<0

\t DTM: Performing 'DISJOINT training sets' k-fold cross validation 'k=cvo.folds'

\t else

\t\t #we have folds>=0
\t\t Performing k-fold cross validation (disjoint validation sets, 'k=cvo.folds')

\t end     

end
"""
mutable struct CVOptions	
	folds::Int #number of cv models
	#if folds<0 then we consider n disjoint training sets
	training_proportion::Float64 #size of trn data
	use_all_data::Bool #if false, the CVOptions will never use the 'original' validation rows in the data
	function CVOptions()
		x=new(10,0.7,true)		
		@assert check_cvoptions(x)
		return x
	end	
	function CVOptions(a::T,b::N,c::Bool) where {T<:Integer,N<:Number}
		x=new(Int(a),float(b),c)
		@assert check_cvoptions(x)
		return x
	end
end

function check_cvoptions(cvo::CVOptions)
	# folds::UInt16 #number of cv models
	# disjoint::Bool #whether we split the data into 'folds' disjoing subsets, if this is false then random sampling is performed
	# training_proportion::Float64 #size of trn data
	# use_all_data::Bool #if false, the CVOptions will never use the 'original' validation rows in the data
	#@assert cvo.folds>0
	# training_proportion == 1 could be possible but is a somewhat degenerate case
	@assert (cvo.training_proportion>=0 && cvo.training_proportion<1)	
	return true
end

function Base.convert(::Type{Splitdef{UInt16}},x::Splitdef{UInt8})
    return Splitdef{UInt16}(x.featid,x.featid_new_positive,x.featurename,Vector{UInt16}(x.subset),x.splitvalue,x.weightl,x.weightr)
end
