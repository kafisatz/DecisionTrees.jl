using BenchmarkTools 

function me1(elementsInLeftChildBV,numerator,weight,current_meanl,current_meanr)
      dl=0.0
      dr=0.0
      f=weight  #simplified version
      for count in 1:length(f)	
      @inbounds idx = count #f.parent.refs[count] + ooo		
      @inbounds ni = numerator[count]
      @inbounds wi = weight[count]
      @inbounds eli=elementsInLeftChildBV[idx]
            if eli            
                  dl += (-log(ni / (wi*current_meanl)) - (ni - wi*current_meanl)/(wi*current_meanl))
            else            
                  dr += (-log(ni / (wi*current_meanr)) - (ni - wi*current_meanr)/(wi*current_meanr))
            end
      end   
      return dl,dr
end


function me2(elementsInLeftChildBV,numerator,weight,current_meanl,current_meanr)
    dl=0.0
      dr=0.0
      f=weight  #simplified version
    for count in 1:length(f)	
        @inbounds idx = count #f.parent.refs[count] + ooo		
        @inbounds ni = numerator[count]
        @inbounds wi = weight[count]
        @inbounds eli=elementsInLeftChildBV[idx]
        wiTimesMean=ifelse(eli,wi*current_meanl,wi*current_meanr)
        addition = -log(ni / wiTimesMean) - (ni - wiTimesMean)/wiTimesMean
          if eli            
                dl += addition
          else            
                dr += addition
          end
    end   
    return dl,dr
end

n=2000000
eLR=(rand(n).<.3)
nu=rand(n)
weight=rand(n)
mr=rand()
ml=rand()
me1(eLR,nu,weight,ml,mr)
me2(eLR,nu,weight,ml,mr)

@assert me1(eLR,nu,weight,ml,mr)==me2(eLR,nu,weight,ml,mr)
@benchmark me1($eLR,$nu,$weight,$ml,$mr)
@benchmark me2($eLR,$nu,$weight,$ml,$mr)
@btime me1($eLR,$nu,$weight,$ml,$mr)
@btime me2($eLR,$nu,$weight,$ml,$mr)
0
0
@code_warntype me2(eLR,nu,weight,ml,mr)