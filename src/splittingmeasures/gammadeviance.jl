
function get_gamma_deviances(current_meanl::Float64,current_meanr::Float64,lo,ooo,f,numerator::Array{Float64,1},denominator::Array{Float64,1},weight::Array{Float64,1},elementsInLeftChildBV)
    #dr and dl are 'reused' by each iteration
	dr=0.0
	dl=0.0
	#note: inbounds increases efficiency here (about a factor of 2), however if the bounds are violated something nasty might happen (quote: If the subscripts are ever out of bounds, you may suffer crashes or silent corruption.)
	for count in 1:length(f)	
		@inbounds idx = f.parent.refs[count] + ooo		
		@inbounds ni = numerator[count]
        @inbounds wi = weight[count]
        if elementsInLeftChildBV[idx]            
            #poisson: dl += (ni*log(ni / (wi*current_meanl)) - (ni - wi*current_meanl))
            dl += (-log(ni / (wi*current_meanl)) - (ni - wi*current_meanl)/(wi*current_meanl))
		else            
            #poisson: dr += (ni*logy(ni / (wi*current_meanr)) - (ni - wi*current_meanr))
            dl += (-log(ni / (wi*current_meanr)) - (ni - wi*current_meanr)/(wi*current_meanr))
		end
	end
	return dl,dr
end
