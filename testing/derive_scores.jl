

#update_current_rels!(currentRelativity,estimatedRatio,trn_meanobservedvalue)	
#warn("remove this")	
td=load("C:\\temp\\sc.jld")	;
#		save("C:\\temp\\sc.jld","trnidx",trnidx,"validx",validx,"estimatedRatio",estimatedRatio,"currentRelativity",currentRelativity,"sett",sett,"actualNumerator",actualNumerator,"denominator",denominator,"weight",weight,"trn_meanobservedvalue",trn_meanobservedvalue,"iter",iter)
#Derive Scores for this iteration #NOTE (tbd/todo: keep in mind) For large datasets (>5m rows) the sorting in construct scores may become a dominant part of the algorithm
trnidx=td["trnidx"];
validx=td["validx"];
weight=td["weight"];
denominator=td["denominator"];
actualNumerator=td["actualNumerator"];
currentRelativity=td["currentRelativity"];
estimatedRatio=td["estimatedRatio"];
iter=td["iter"];
sett=td["sett"];
trn_meanobservedvalue=td["trn_meanobservedvalue"];

maxRawRelativityPerScoreSorted,MAPPINGSmoothedEstimatePerScore,vectorWeightPerScore,obsPerScore,rawObservedRatioPerScore,numPerScore,denomPerScore,uniqueRelativitiesSorted,nscoresPotentiallyReduced				=DTM2.constructANDderiveScores(trnidx,validx,estimatedRatio,currentRelativity,actualNumerator,denominator,weight,trn_meanobservedvalue,iter,sett)
@assert isequal(DTM2.constructANDderiveScores_OLD_VERSION(trnidx,validx,estimatedRatio,currentRelativity,sett.nscores,actualNumerator,denominator,weight,trn_meanobservedvalue,iter,sett.smoothEstimates, sett.print_details),DTM2.constructANDderiveScores(trnidx,validx,estimatedRatio,currentRelativity,actualNumerator,denominator,weight,trn_meanobservedvalue,iter,sett))
@benchmark DTM2.constructANDderiveScores(trnidx,validx,estimatedRatio,currentRelativity,actualNumerator,denominator,weight,trn_meanobservedvalue,iter,sett)
#v0 memory estimate:  20.59 MiB  allocs estimate:  164531   179ms 2.4%gc
#v1 

rawObservedRatioPerScore[i]=sum(view(numerator_srt,previdx:nextidx))/sum(view(denominator_srt,previdx:nextidx))

@code_warntype DTM2.constructANDderiveScores(trnidx,validx,estimatedRatio,currentRelativity,actualNumerator,denominator,weight,trn_meanobservedvalue,iter,sett)

relativities=currentRelativity
meanobservedvalue=trn_meanobservedvalue
nscores=sett.nscores
smooth=sett.smoothEstimates
print_details=sett.print_details
boolSmoothingEnabled=!(smooth=="0")
#derive Scores from relativities
uniqueRelativitiesSorted=unique(view(relativities,trnidx))
sort!(uniqueRelativitiesSorted,alg=QuickSort)
ns=min(nscores,size(uniqueRelativitiesSorted,1))
#ns<nscores&&warn("There are only $(ns) distinct relativities. Number of scores will be limited to $(ns) instead of $(nscores).")
aggregatedModelledRatioPerScore,maxRawRelativityPerScoreSorted,MAPPINGSmoothedEstimatePerScore,vectorWeightPerScore,obsPerScore,rawObservedRatioPerScore,numPerScore,denomPerScore=DTM2.constructScores(
    trnidx,validx,estimatedRatio,relativities,uniqueRelativitiesSorted,actualNumerator,denominator,weight,meanobservedvalue,ns,nscores,print_details,boolSmoothingEnabled=boolSmoothingEnabled)	

    @code_warntype DTM2.constructScores(trnidx,validx,estimatedRatio,relativities,uniqueRelativitiesSorted,actualNumerator,denominator,weight,meanobservedvalue,ns,nscores,print_details,boolSmoothingEnabled=boolSmoothingEnabled)	

    estimatedRatioPerRow_full_vec=estimatedRatio
    raw_rel_full_vec=relativities
    numerator_full_vec=actualNumerator
    denominator_full_vec=denominator
    weight_full_vec=weight
    nscoresPotentiallyReduced=ns


    import DTM2.deriveRowCountAndWeight,DTM2.deriveRowCountAndWeight,DTM2.cumulativeToIncremental!, DTM2.aggregate_values_per_score, DTM2.smooth_scores, DTM2.insert_gaps_to_scores!,DTM2.derive_scores_main_aggregation_step!


    function mxi()
        x=rand(999)
        srtt=sortperm(x)
        v= view(x,srtt)
       
        return nothing
    end

    mxi()




    x=rand(999)
    srtt=sortperm(x)
    @btime view(x,srtt)





    
function test1(raw_rel)
	srt::Vector{Int}=sortperm(raw_rel,alg=QuickSort) #we could use RadixSort here, it uses about 40m of memory for 5m obs, AND I am unsure how easy it is to get the permutation (compared to sorting in place). But the algorithm is about 50% faster for random floats than QuickSort.
    raw_rel_srt=view(raw_rel,srt)
return nothing
end

x=rand(999);
@benchmark test1(x)
@benchmark sortperm(x)
@benchmark sortperm(x,alg=QuickSort)







estimatedRatioPerRow_full_vec::Array{Float64,1},raw_rel_full_vec::Array{Float64,1},uniqueRelativitiesSorted::Array{Float64,1},numerator_full_vec::Array{Float64,1},denominator_full_vec::Array{Float64,1},weight_full_vec::Array{Float64,1},meanobservedvalue::Float64,nscoresPotentiallyReduced::Int,nscores::Int,print_details::Bool;boolSmoothingEnabled::Bool=true)
trnidx,validx,estimatedRatio,relativities,uniqueRelativitiesSorted,actualNumerator,denominator,weight,meanobservedvalue,ns,nscores,print_details,boolSmoothingEnabled=boolSmoothingEnabled)	

estimatedRatioPerRow_full_vec=estimatedRatio
raw_rel_full_vec=relativities
numerator_full_vec=actualNumerator
denominator_full_vec=denominator
weight_full_vec=weight
nscoresPotentiallyReduced=ns


@assert nscores>=nscoresPotentiallyReduced
maxItersForSmoothing=0
#todo,tbd this can probably be done more efficiently
raw_rel=view(raw_rel_full_vec,trnidx)
weight=view(weight_full_vec,trnidx)
numerator=view(numerator_full_vec,trnidx)
denominator=view(denominator_full_vec,trnidx)
estimatedRatioPerRow=view(estimatedRatioPerRow_full_vec,trnidx)

srt=sortperm(raw_rel,alg=QuickSort) #we could use RadixSort here, it uses about 40m of memory for 5m obs, AND I am unsure how easy it is to get the permutation (compared to sorting in place). But the algorithm is about 50% faster for random floats than QuickSort.
obs=length(raw_rel)
raw_rel_srt=view(raw_rel,srt)
weight_srt=view(weight,srt)
estimatedRatioPerRow_srt=view(estimatedRatioPerRow,srt)
denominator_srt=view(denominator,srt)
numerator_srt=view(numerator,srt)

#=
@code_warntype deriveRowCountAndWeight(weight_srt,raw_rel_srt,uniqueRelativitiesSorted)
=#
rowCountPerUniqueRelativity,weightPerUniqueRelativity=deriveRowCountAndWeight(weight_srt,raw_rel_srt,uniqueRelativitiesSorted)

numeratorEstimatedPerRow_srt=estimatedRatioPerRow_srt.*denominator_srt
wtot=sum(weight_srt)
wperscore=wtot/nscoresPotentiallyReduced
scoreEndPoints=zeros(Int,nscoresPotentiallyReduced)
scoreEndPoints[end]=size(weight_srt,1)
vectorWeightPerScore=Array{Float64}(nscoresPotentiallyReduced)
cumulativeNumberOfDistinctRawRelativitesPerScore=zeros(Int,nscoresPotentiallyReduced)

#=
    @code_warntype derive_scores_main_aggregation_step!(wperscore,uniqueRelativitiesSorted,weightPerUniqueRelativity,rowCountPerUniqueRelativity,cumulativeNumberOfDistinctRawRelativitesPerScore,scoreEndPoints,vectorWeightPerScore)
=#
nscoresPotentiallyReducedTWOTimes=derive_scores_main_aggregation_step!(wperscore,uniqueRelativitiesSorted,weightPerUniqueRelativity,rowCountPerUniqueRelativity,cumulativeNumberOfDistinctRawRelativitesPerScore,scoreEndPoints,vectorWeightPerScore)

obsPerScore=copy(scoreEndPoints)
cumulativeToIncremental!(obsPerScore)
cumulativeToIncremental!(vectorWeightPerScore)

@code_warntype cumulativeToIncremental!(vectorWeightPerScore)

#aggregate scores
@code_warntype aggregate_values_per_score(nscoresPotentiallyReducedTWOTimes,scoreEndPoints,raw_rel_srt,numerator_srt,denominator_srt,obs,numeratorEstimatedPerRow_srt)
numPerScore,denomPerScore,maxRawRelativityPerScoreSorted,aggregatedModelledRatioPerScore,rawObservedRatioPerScore=aggregate_values_per_score(nscoresPotentiallyReducedTWOTimes,scoreEndPoints,raw_rel_srt,numerator_srt,denominator_srt,obs,numeratorEstimatedPerRow_srt)
#smooth scores
estimatedRatioPerScore=smooth_scores(rawObservedRatioPerScore,print_details,boolSmoothingEnabled)	
@code_warntype smooth_scores(rawObservedRatioPerScore,print_details,boolSmoothingEnabled)	
#insert gaps if necessary
insert_gaps_to_scores!(wtot,nscores,nscoresPotentiallyReducedTWOTimes,aggregatedModelledRatioPerScore,maxRawRelativityPerScoreSorted,estimatedRatioPerScore,vectorWeightPerScore,obsPerScore,rawObservedRatioPerScore,numPerScore,denomPerScore)



nscoresPotentiallyReduced=39
n=1500000

wei=rand(n)*400
wperscore=sum(wei)/nscoresPotentiallyReduced
rel=rand(n)+rand(n).*2.4+1.0
srt=sortperm(rel)
relativitiesSorted=rel[srt]
weight_srt=wei[srt]
scoreEndPoints=zeros(Int,nscoresPotentiallyReduced)
scoreEndPoints[end]=size(weight_srt,1)
vectorWeightPerScore=Array{Float64}(nscoresPotentiallyReduced)

    DTM2.derive_scores_main_aggregation_step!(nscoresPotentiallyReduced,wperscore,relativitiesSorted,weight_srt,scoreEndPoints,vectorWeightPerScore)
    @code_warntype  DTM2.derive_scores_main_aggregation_step!(nscoresPotentiallyReduced,wperscore,relativitiesSorted,weight_srt,scoreEndPoints,vectorWeightPerScore)

    @show vectorWeightPerScore
    @show scoreEndPoints
    xxx=deepcopy(vectorWeightPerScore);DTM2.cumulativeToIncremental!(xxx);xxx


    @show xxx.-weight_srt
    median(xxx)
    maximum(xxx)

    relativitiesSorted[95002]