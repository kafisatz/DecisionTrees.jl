
function calculateSplitValue(a::RankOptSplit,labels::Array{Float64,1}, labels_orig::Array{Float64,1}, labels_new::Array{Float64,1}, features::Array{Float64,1}, thresh::Float64,minweight::Float64,prem_buffer::Int64,moderationfactor::Float64)
   res=all_minus=all_minus=all_akt=l_plus=l_minus=l_akt=r_plus=r_minus=r_akt=l_orig=l_new=meanl=meanr=l=f=suml=sumr=l_plus=l_minus=r_plus=r_minus=all_plus=all_minus=0.0
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
    l_orig=labels_orig[i]
    l_new=labels_new[i]
    f=features[i]
    all_akt=max(l_orig,l_new-Float64(prem_buffer))-l
      if f<thresh
          l_akt=max(l_orig,l_new-Float64(prem_buffer))-l+meanl*moderationfactor
          l_plus=l_plus+min(l_new-l_orig,l_akt-l_orig)
          l_minus=l_minus-max(l_akt-l_new,0)
      else
          r_akt=max(l_orig,l_new-Float64(prem_buffer))-l+meanr*moderationfactor
          r_plus=r_plus+min(l_new-l_orig,r_akt-l_orig)
          r_minus=r_minus-max(r_akt-l_new,0)
      end
      all_plus=all_plus+min(l_new-l_orig,all_akt-l_orig)
      all_minus=all_minus-max(all_akt-l_new,0)
    end
    res=l_plus+l_minus+r_plus+r_minus-(all_plus+all_minus)

    return res, Float64(countl),Float64(countr) #sum(weight[this_split]), sum(weight[!this_split]
end

function calculateSplitValue(a::RankOptSplit,labels::Array{Float64,1}, labels_orig::Array{Float64,1}, labels_new::Array{Float64,1}, features::pdaMod,subset::Array{UInt8,1},minweight::Float64,prem_buffer::Int64, moderationfactor::Float64)
   res=all_minus=all_minus=all_akt=l_plus=l_minus=l_akt=r_plus=r_minus=r_akt=l_orig=l_new=meanl=meanr=l=f=suml=sumr=l_plus=l_minus=r_plus=r_minus=all_plus=all_minus=0.0
    countl=countr=0
   for i in 1:length(features.pda)

      if features.pda.refs[i] in subset
        suml+=labels[i]
        countl+=1
      else
        sumr+=labels[i]
        countr+=1
      end
    end
    meanl=sumr/Float64(countr)
    meanr=suml/Float64(countl)
    for i in 1:length(features.pda)
    l=labels[i]
    l_orig=labels_orig[i]
    l_new=labels_new[i]
    all_akt=max(l_orig,l_new-Float64(prem_buffer))-l
      if features.pda.refs[i] in subset
          l_akt=max(l_orig,l_new-Float64(prem_buffer))-l+meanl*moderationfactor
          l_plus=l_plus+min(l_new-l_orig,l_akt-l_orig)
          l_minus=l_minus-max(l_akt-l_new,0)
      else
          r_akt=max(l_orig,l_new-Float64(prem_buffer))-l+meanr*moderationfactor
          r_plus=r_plus+min(l_new-l_orig,r_akt-l_orig)
          r_minus=r_minus-max(r_akt-l_new,0)
      end
      all_plus=all_plus+min(l_new-l_orig,all_akt-l_orig)
      all_minus=all_minus-max(all_akt-l_new,0)
    end
    res=l_plus+l_minus+r_plus+r_minus-(all_plus+all_minus)

    return res,Float64(countl),Float64(countr) #,w1,(sum(weight)-w1)
end