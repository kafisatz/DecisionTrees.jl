dd=load("c:\\temp\\j.jld");
scores=dd["scores"];
rawObservedRatioPerScore=dd["rawObservedRatioPerScore"];
MAPPINGSmoothedEstimatePerScore=dd["MAPPINGSmoothedEstimatePerScore"];

@btime estimateUnsmoothed=Float64[rawObservedRatioPerScore[x] for x in scores]
@btime estimateUnsmoothed2=map(x->float(rawObservedRatioPerScore[x]),scores)
@btime estimateUnsmoothed3=map(x->rawObservedRatioPerScore[x],scores)
estimateUnsmoothed=Float64[rawObservedRatioPerScore[x] for x in scores]
estimateUnsmoothed2=map(x->float(rawObservedRatioPerScore[x]),scores)
estimateUnsmoothed3=map(x->rawObservedRatioPerScore[x],scores)
all(estimateUnsmoothed.==estimateUnsmoothed2)
all(estimateUnsmoothed.==estimateUnsmoothed3)

estimateSmoothed=Float64[MAPPINGSmoothedEstimatePerScore[x] for x in scores]	
#warn("remove this")
#	save("C:\\temp\\j.jld","MAPPINGSmoothedEstimatePerScore",MAPPINGSmoothedEstimatePerScore,"scores",scores,"rawObservedRatioPerScore",rawObservedRatioPerScore)
