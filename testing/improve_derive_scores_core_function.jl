using StatsBase
using BenchmarkTools


function cumulativeToIncremental!(x::Vector{T}) where {T <: Number}
	for i=length(x):-1:2
		x[i]-=x[i-1]
	end
	return nothing
end

function incrementalToCumulative!(x::Vector{T}) where {T <: Number}
	for i=2:length(x)
		x[i]+=x[i-1]
	end
	return nothing
end

function derive_scores_main_aggregation_step!(nScores,wperscore,relativitiesSorted,weight_srt,scoreEndPoints,vectorWeightPerScore)
#todo/tbd this can be simplified if length(uniqueRelativitiesSorted)<sett.nScores then we do not need to do any aggregation at all.

#NOTE: (todo/tbd improve this) there are different ways to aggregate the relativities to scores. the approach below is a bottom up approach
#issues can arise when there is a mass (of weight) for a certain relativity towards the end of the scores
#because the approach below assigns at least one relativity to each score "one can run out of data" before all scores are filled up (in certain degenerate cases)
#an alternative would be a "binary" split starting in the middle

#note: if the estimates and the weight are degenerated then the vector scoreEndPoints can have the same value multiple times (i.e. it is constant for a few indices)
uniqueRelativitiesSortedCounter=0
current_cumulW=0.0
startPoint=1
# assert size(rowCountPerUniqueRelativity)==size(weightPerUniqueRelativity) #this is to avoid errors due to the inbounds macro below
stillRunning=true
#wi=weightPerUniqueRelativity[1]
wi=0.0
current_score=0
current_relativity=0.0
current_idx=0 #rowCountPerUniqueRelativity[1]

while stillRunning
current_score+=1
  	while current_idx<length(weight_srt) 
		current_idx+=1
		last_relativity=current_relativity
		@inbounds current_relativity=relativitiesSorted[current_idx]
		@inbounds wi::Float64=weight_srt[current_idx]::Float64 
		current_cumulW::Float64 += wi::Float64
		if current_cumulW>wperscore*current_score
			#-> weight per score is reached
			#we now need to ensure that the 'cut' is made at a point where the relativity changes (otherwise it will be impossible to define the scores in relation to the (moderated) relativities)
			while current_relativity==last_relativity
				current_idx+=1
				if current_idx>length(relativitiesSorted)
					break
				end
				last_relativity=current_relativity				
				@inbounds current_relativity=relativitiesSorted[current_idx]
				@inbounds wi=weight_srt[current_idx]::Float64
				current_cumulW += wi::Float64		
			end
			break
		end
 	end
	@inbounds vectorWeightPerScore[current_score]=current_cumulW
	@inbounds scoreEndPoints[current_score]=current_idx
	stillRunning=current_score<length(scoreEndPoints)
end

nscoresPotentiallyReducedTWOTimes::Int=0
if (length(scoreEndPoints)>1)&&(scoreEndPoints[end-1]>=scoreEndPoints[end])
#warn("Failed to uniformly distribute the scores! Scores will be degenerated. \r\nThis usually indicates that: (i) there are too many scores compared to the granularity of the model or (ii) the data is 'degenerated' (i.e. there is a mass of exposure for certain risk details)")
#Reduce the number of Scores again
	thismax=maximum(scoreEndPoints)
	endlocation::Int=searchsortedfirst(scoreEndPoints,convert(eltype(scoreEndPoints),thismax))
	resize!(scoreEndPoints,endlocation) #drop the tail end
	resize!(vectorWeightPerScore,endlocation) #drop the tail end
	nscoresPotentiallyReducedTWOTimes=length(scoreEndPoints)
else
	nscoresPotentiallyReducedTWOTimes=nScores::Int
end

scoreEndPoints[end]=min(length(weight_srt),scoreEndPoints[end])

return nscoresPotentiallyReducedTWOTimes::Int
end

#generate some data
nnn=20
relativitiesSorted=sample(rand(nnn),nnn,replace=true)
sort!(relativitiesSorted)
weight_srt=201*rand(nnn)+130*rand(nnn)

nScores=floor(Int,min(100,nnn/2))
scoreEndPoints=zeros(Int,nScores)
obsPerScore=zeros(Int,nScores)
scoreEndPoints[end]=size(weight_srt,1)
vectorWeightPerScore=zeros(Float64,nScores) #Array{Float64}(undef,nScores)
wtot=sum(weight_srt)
wperscore=wtot/nScores

#@btime derive_scores_main_aggregation_step!(nScores,wperscore,relativitiesSorted,weight_srt,scoreEndPoints,vectorWeightPerScore)

#@code_warntype derive_scores_main_aggregation_step!(nScores,wperscore,relativitiesSorted,weight_srt,scoreEndPoints,vectorWeightPerScore)

#@btime quantile(relativitiesSorted,StatsBase.fweights(weight_srt),range(1/nScores,stop=1-1/nScores,length=nScores))

qtls=quantile(relativitiesSorted,StatsBase.fweights(weight_srt),range(1/nScores,stop=1,length=nScores))
qtls=unique(qtls)
nscoresPotentiallyReducedTWOTimes=length(qtls)
obsPerScore=zeros(Int,nscoresPotentiallyReducedTWOTimes)
scoreEndPoints=zeros(Int,nscoresPotentiallyReducedTWOTimes)
vectorWeightPerScore=zeros(Float64,nscoresPotentiallyReducedTWOTimes)

thisI=1
thisJ=1
while thisI<=length(relativitiesSorted)
    @inbounds reli=relativitiesSorted[thisI]
    @inbounds wi=weight_srt[thisI]
    if reli<=qtls[thisJ]
        obsPerScore[thisJ]+=1
        vectorWeightPerScore[thisJ]+=wi        
    else 
        if thisJ<length(obsPerScore)
            thisJ+=1
        end
        obsPerScore[thisJ]+=1
        vectorWeightPerScore[thisJ]+=wi
    end 
    scoreEndPoints[thisJ]=thisI 
    thisI+=1
end

#check consistency of scoreEndPoints
for i=1:length(scoreEndPoints)
    @show sum(view(weight_srt,1:scoreEndPoints[i]))-sum(vectorWeightPerScore2[1:i])
end


obsPerScore2=deepcopy(obsPerScore)
vectorWeightPerScore2=deepcopy(vectorWeightPerScore)
scoreEndPoints2=deepcopy(scoreEndPoints)

    
obsPerScore=copy(scoreEndPoints)
cumulativeToIncremental!(obsPerScore)
cumulativeToIncremental!(vectorWeightPerScore)

obsPerScore2.-obsPerScore
vectorWeightPerScore2.-vectorWeightPerScore
scoreEndPoints2.-scoreEndPoints

#=
include(joinpath(pwd(),"testing","improve_derive_scores_core_function.jl"))
=#