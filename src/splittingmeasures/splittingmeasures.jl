#include("splitting_measures_mse.jl")
include("difference.jl")
include("normalDeviance.jl")
include("poissonGammaNormal.jl")

#these are not maintained any more
include("maxsplitvalue.jl")
include("maxabssplitvalue.jl")
include("maxminussplitvalue.jl")

#these are not maintained and most most certainly not working .
#include("ROptMinRLost.jl")
#include("ROptMinRLostPct.jl")
#include("rankopt.jl")

#=
#some old code snippets
function _minimize_lost_ranks_maximize_5pct(labels::Array{Float64,1}, labels_orig::Array{Float64,1}, labels_new::Array{Float64,1}, features::Array{Float64,1}, thresh::Float64,minweight::Float64,prem_buffer::Int,moderationfactor::Float64)
    sz=size(labels,1)
    error("this is not working anymore")
	idx=features.<thresh
    labelsl=labels[idx] #this creates a copy of the array which is very inefficient
    labelsr=labels[!idx] #this creates a copy of the array which is very inefficient
    countl=sum(idx)
    countr=sz-countl
return max(quantile(labelsl,0.05),quantile(labelsr,0.05)),Float64(countl),Float64(countr)
end

function _minimize_lost_ranks_maximize_minresidual(labels::Array{Float64,1}, labels_orig::Array{Float64,1}, labels_new::Array{Float64,1}, features::Array{Float64,1}, thresh::Float64,minweight::Float64,prem_buffer::Int,moderationfactor::Float64)
    minl=minr=Inf
    countl=countr=0
	f,l=0.0
    for i=1:size(features,1)
    l=labels[i]
    f=features[i]
      if f<thresh
        if l<minl;minl=l;end
        countl+=1
      else
        if l<minr;minr=l;end
        countr+=1
      end
    end
return max(minr,minl),Float64(countl),Float64(countr)
end

function _minimize_lost_ranks_maximize_minresidual(labels::Array{Float64,1}, labels_orig::Array{Float64,1}, labels_new::Array{Float64,1}, features,subset::Array{UInt8,1},minweight::Float64,prem_buffer::Int, moderationfactor::Float64)
minl=minr=l=0.0
  countl=countr=0
    for i in 1:length(features.pda)
        l=labels[i]
        if features.pda.refs[i] in subset
            if l<minl;minl=l;end
            countl += 1
        else
            if l<minr;minr=l;end
            countr += 1
        end
    end
return max(minr,minl),Float64(countl),Float64(countr)
end

function _minimize_lost_ranks_difference(labels::Array{Float64,1}, labels_orig::Array{Float64,1}, labels_new::Array{Float64,1}, features::Array{Float64,1}, thresh::Float64,minweight::Float64,prem_buffer::Int,moderationfactor::Float64)
 suml=sumr=l=f=0.0
    countl=countr=0
    for i=1:size(features,1)
    l=labels[i]
    f=features[i]
      if f<thresh
        suml+=l
        countl+=1
      else
        sumr+=l
        countr+=1
      end
    end
return abs(sumr/Float64(countr)-suml/Float64(countl)), Float64(countl),Float64(countr)
end

function _minimize_lost_ranks_difference(labels::Array{Float64,1}, labels_orig::Array{Float64,1}, labels_new::Array{Float64,1}, features,subset::Array{UInt8,1},minweight::Float64,prem_buffer::Int, moderationfactor::Float64)
suml=sumr=0.0
  countl=countr=0
    for i in 1:length(features.pda)
        if features.pda.refs[i] in subset
            suml += labels[i]
            countl += 1
        else
            sumr += labels[i]
            countr += 1
        end
    end
return abs(sumr/Float64(countr)-suml/Float64(countl)), Float64(countl),Float64(countr)
end

function _min_mean(labels::Array{Float64,1}, features::Array{Float64,1}, thresh::Float64,minweight::Float64)
    l=f=suml=sumr=0.0
    countl=countr=0
    for i=1:size(features,1)
    l=labels[i]
    f=features[i]
      if f<thresh
        suml+=l
        countl+=1
      else
        sumr+=l
        countr+=1
      end
    end
return min(sumr/Float64(countr),suml/Float64(countl)), Float64(countl),Float64(countr)
end

function _min_mean(labels::Array{Float64,1}, features,subset::Array{UInt8,1},minweight::Float64)
  suml=sumr=0.0
  countl=countr=0
    for i in 1:length(features.pda)
        if features.pda.refs[i] in subset
            suml += labels[i]
            countl += 1
        else
            sumr += labels[i]
            countr += 1
        end
    end
return min(sumr/Float64(countr),suml/Float64(countl)), Float64(countl),Float64(countr)
end

function _max_mean(labels::Array{Float64,1}, features::Array{Float64,1}, thresh::Float64,minweight::Float64)
    l=f=suml=sumr=0.0
    countl=countr=0
    for i=1:size(features,1)
    l=labels[i]
    f=features[i]
      if f<thresh
        suml+=l
        countl+=1
      else
        sumr+=l
        countr+=1
      end
    end
return max(sumr/Float64(countr),suml/Float64(countl)), Float64(countl),Float64(countr)
end

function _max_mean(labels::Array{Float64,1}, features,subset::Array{UInt8,1},minweight::Float64)
  suml=sumr=0.0
  countl=countr=0
    for i in 1:length(features.pda)
        if features.pda.refs[i] in subset
            suml += labels[i]
            countl += 1
        else
            sumr += labels[i]
            countr += 1
        end
    end
return max(sumr/Float64(countr),suml/Float64(countl)), Float64(countl),Float64(countr)
end



function _min_std(labels::Array{Float64,1}, features::Array{Float64,1}, thresh::Float64,minweight::Float64)
    x=f=detal=detanl=terml=m1l=m2l=detar=detanr=termr=m1r=m2r=0.0
    countl=countr=0
    for i=1:size(features,1)
    x=labels[i]
    f=features[i]
      if f<thresh
        countl+=1;detal=x-m1l;detanl=detal/countl;terml=detal*detanl*(countl-1);m1l+=detanl;m2l+=terml;
      else
        countr+=1;detar=x-m1r;detanr=detar/countr;termr=detar*detanr*(countr-1);m1r+=detanr;m2r+=termr;
      end
    end
    varl=m2l/(countl-1)
    varr=m2r/(countr-1)

  count=countl+countr
  deta=m1l-m1r
  deta2=deta*deta
  m2=m2l+m2r+deta2*countl*countr/count
  vartot=m2/(count-1)
  return min(sqrt(varl),sqrt(varr)), Float64(countl),Float64(countr)
end

function _min_std(labels::Array{Float64,1}, features,subset::Array{UInt8,1},minweight::Float64)
 x=detal=detanl=terml=m1l=m2l=detar=detanr=termr=m1r=m2r=0.0
   countl=countr=0
   for i in 1:length(features.pda)
        x=labels[i]
        if features.pda.refs[i] in subset
        countl+=1;detal=x-m1l;detanl=detal/countl;terml=detal*detanl*(countl-1);m1l+=detanl;m2l+=terml;
      else
        countr+=1;detar=x-m1r;detanr=detar/countr;termr=detar*detanr*(countr-1);m1r+=detanr;m2r+=termr;
      end
    end
   varl=m2l/(countl-1)
    varr=m2r/(countr-1)

  count=countl+countr
  deta=m1l-m1r
  deta2=deta*deta
  m2=m2l+m2r+deta2*countl*countr/count
  vartot=m2/(count-1)
  return min(sqrt(varl),sqrt(varr)), Float64(countl),Float64(countr)
end

function _reduction_in_standarddeviation(labels::Array{Float64,1}, features::Array{Float64,1}, thresh::Float64,minweight::Float64)
    x=f=detal=detanl=terml=m1l=m2l=detar=detanr=termr=m1r=m2r=0.0
    countl=countr=0
    for i=1:size(features,1)
    x=labels[i]
    f=features[i]
      if f<thresh
        countl+=1;detal=x-m1l;detanl=detal/countl;terml=detal*detanl*(countl-1);m1l+=detanl;m2l+=terml;
      else
        countr+=1;detar=x-m1r;detanr=detar/countr;termr=detar*detanr*(countr-1);m1r+=detanr;m2r+=termr;
      end
    end
  varl=m2l/(countl-1)
  varr=m2r/(countr-1)

  count=countl+countr
  deta=m1l-m1r
  deta2=deta*deta
  m2=m2l+m2r+deta2*countl*countr/count
  vartot=m2/(count-1)
  return (sqrt(vartot)-(sqrt(varl)+sqrt(varr))), Float64(countl),Float64(countr)
end

function _reduction_in_standarddeviation(labels::Array{Float64,1}, features,subset::Array{UInt8,1},minweight::Float64)
   x=detal=detanl=terml=m1l=m2l=detar=detanr=termr=m1r=m2r=0.0
   countl=countr=0
   for i in 1:length(features.pda)
        x=labels[i]
        if features.pda.refs[i] in subset
        countl+=1;detal=x-m1l;detanl=detal/countl;terml=detal*detanl*(countl-1);m1l+=detanl;m2l+=terml;
      else
        countr+=1;detar=x-m1r;detanr=detar/countr;termr=detar*detanr*(countr-1);m1r+=detanr;m2r+=termr;
      end
    end
    varl=m2l/(countl-1)
  varr=m2r/(countr-1)

  count=countl+countr
  deta=m1l-m1r
  deta2=deta*deta
  m2=m2l+m2r+deta2*countl*countr/count
  vartot=m2/(count-1)
  return (sqrt(vartot)-(sqrt(varl)+sqrt(varr))), Float64(countl),Float64(countr)
end

function _reduction_in_variance(labels::Array{Float64,1}, features::Array{Float64,1}, thresh::Float64,minweight::Float64)
    x=f=detal=detanl=terml=m1l=m2l=detar=detanr=termr=m1r=m2r=0.0
    countl=countr=0
    for i=1:size(features,1)
    x=labels[i]
    f=features[i]
      if f<thresh
        countl+=1;detal=x-m1l;detanl=detal/countl;terml=detal*detanl*(countl-1);m1l+=detanl;m2l+=terml;
      else
        countr+=1;detar=x-m1r;detanr=detar/countr;termr=detar*detanr*(countr-1);m1r+=detanr;m2r+=termr;
      end
    end
    varl=m2l/(countl-1)
    varr=m2r/(countr-1)

  count=countl+countr
  deta=m1l-m1r
  deta2=deta*deta
  m2=m2l+m2r+deta2*countl*countr/count
  vartot=m2/(count-1)
  return (vartot-(varl+varr)), Float64(countl),Float64(countr)
end

function _reduction_in_variance(labels::Array{Float64,1}, features,subset::Array{UInt8,1},minweight::Float64)
   x=detal=detanl=terml=m1l=m2l=detar=detanr=termr=m1r=m2r=0.0
   countl=countr=0
   for i in 1:length(features.pda)
        x=labels[i]
        if features.pda.refs[i] in subset
        countl+=1;detal=x-m1l;detanl=detal/countl;terml=detal*detanl*(countl-1);m1l+=detanl;m2l+=terml;
      else
        countr+=1;detar=x-m1r;detanr=detar/countr;termr=detar*detanr*(countr-1);m1r+=detanr;m2r+=termr;
      end
    end
     varl=m2l/(countl-1)
    varr=m2r/(countr-1)

  count=countl+countr
  deta=m1l-m1r
  deta2=deta*deta
  m2=m2l+m2r+deta2*countl*countr/count
  vartot=m2/(count-1)
  return (vartot-(varl+varr)), Float64(countl),Float64(countr)
end


function _mean_abs_error(labels::Array{Float64,1}, features::Array{Float64,1}, thresh::Float64,minweight::Float64)
    #todo/tbd this can most likely be done with only one pass over the data!
    meanl=meanr=l=f=res=suml=sumr=0.0
    countl=countr=0
    for i=1:size(features,1)
    l=labels[i]
    f=features[i]
      if f<thresh
        suml+=l
        countl+=1
      else
        sumr+=l
        countr+=1
      end
    end
    meanl=sumr/Float64(countr)
    meanr=suml/Float64(countl)
    for i=1:size(features,1)
      l=labels[i]
      f=features[i]
      if f<thresh
        res+=abs(meanl-l)
      else
        res+=abs(meanr-l)
      end
    end
    return res/Float64(countl+countr), Float64(countl),Float64(countr)
end

function _mean_abs_error(labels::Array{Float64,1}, features,subset::Array{UInt8,1},minweight::Float64)
    #NOTE: This function will yield the same resulting tree as _sum_abs_error
  #todo/tbd this can most likely be done with only one pass over the data!

	meanl=meanr=l=res=suml=sumr=0.0
	countl=countr=0
    for i in 1:length(labels)
		l=labels[i]
		f=features.pda.refs[i]
      if in(f,subset)
        suml+=l
        countl+=1
      else
        sumr+=l
        countr+=1
      end
    end
    meanl=sumr/Float64(countr)
    meanr=suml/Float64(countl)
    for i in 1:length(labels)
      l=labels[i]
	  f=features.pda.refs[i]
      if in(f,subset)
        res+=abs(meanl-l)
      else
        res+=abs(meanr-l)
      end
    end
    return res/Float64(countl+countr), Float64(countl),Float64(countr)
end

function _mean_abs_error_rel(labels::Array{Float64,1}, features::Array{Float64,1}, thresh::Float64,minweight::Float64)
    #NOTE: This function will yield the same resulting tree as _sum_abs_error
  #todo/tbd this can most likely be done with only one pass over the data!
    meanl=meanr=l=f=res=suml=sumr=0.0
    countl=countr=0
    for i=1:size(features,1)
    l=labels[i]
    f=features[i]
      if f<thresh
        suml+=l
        countl+=1
      else
        sumr+=l
        countr+=1
      end
    end
    meanl=sumr/Float64(countr)
    meanr=suml/Float64(countl)
    for i=1:size(features,1)
      l=labels[i]
      f=features[i]
      if f<thresh
        res+=abs(meanl-l)/l
      else
        res+=abs(meanr-l)/l
      end
    end
    return res/Float64(countl+countr), Float64(countl),Float64(countr)
end

function _mean_abs_error_rel(labels::Array{Float64,1},  features,subset::Array{UInt8,1},minweight::Float64)
	#todo/tbd this can most likely be done with only one pass over the data!

  res=meanr=meanl=l=suml=sumr=0.0
	countl=countr=0
    for i in 1:length(labels)
		l=labels[i]
		f=features.pda.refs[i]
      if in(f,subset)
        suml+=l
        countl+=1
      else
        sumr+=l
        countr+=1
      end
    end
    meanl=sumr/Float64(countr)
    meanr=suml/Float64(countl)
    for i in 1:length(labels)
      l=labels[i]
	  f=features.pda.refs[i]
      if in(f,subset)
        res+=abs(meanl-l)/l
      else
        res+=abs(meanr-l)/l
      end
    end
    return res/Float64(countl+countr), Float64(countl),Float64(countr)
end

function _sum_abs_error(labels::Array{Float64,1}, features::Array{Float64,1}, thresh::Float64,minweight::Float64)
  #NOTE: This function will yield the same resulting tree as _mean_abs_error
  #todo/tbd this can most likely be done with only one pass over the data!
    meanl=meanr=l=f=res=suml=sumr=0.0
    countl=countr=0
    for i=1:size(features,1)
    l=labels[i]
    f=features[i]
      if f<thresh
        suml+=l
        countl+=1
      else
        sumr+=l
        countr+=1
      end
    end
    meanl=sumr/Float64(countr)
    meanr=suml/Float64(countl)
    for i=1:size(features,1)
      l=labels[i]
      f=features[i]
      if f<thresh
        res+=abs(meanl-l)
      else
        res+=abs(meanr-l)
      end
    end
    return res, Float64(countl),Float64(countr)
end

function _sum_abs_error(labels::Array{Float64,1},features,subset::Array{UInt8,1},minweight::Float64)
  #NOTE: This function will yield the same resulting tree as _mean_abs_error
  #todo/tbd this can most likely be done with only one pass over the data!

  res=meanr=meanl=l=suml=sumr=0.0
	countl=countr=0
    for i in 1:length(labels)
		l=labels[i]
		f=features.pda.refs[i]
      if in(f,subset)
        suml+=l
        countl+=1
      else
        sumr+=l
        countr+=1
      end
    end
    meanl=sumr/Float64(countr)
    meanr=suml/Float64(countl)
    for i in 1:length(labels)
      l=labels[i]
	  f=features.pda.refs[i]
      if in(f,subset)
        res+=abs(meanl-l)
      else
        res+=abs(meanr-l)
      end
    end
    return res, Float64(countl),Float64(countr)
end

=#