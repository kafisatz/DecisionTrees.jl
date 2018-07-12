
mutable struct CustomVariance
    meanratio::Float64
    sn::Float64 #sum nominator
    sd::Float64 #sum denominator    
    n::Int #count
    m2::Float64 # var*(n-1)
    var::Float64 #weigthedVariance
    function CustomVariance(num,denom)
        n=length(num)
        @assert n==length(denom) # length(w)
        sn=sum(num)
        sd=sum(denom)        
        meanratio=sn/sd
        m2=zero(eltype(num))           
        for i=1:n
            @inbounds m2 += denom[i]*(num[i]/denom[i]-meanratio)^2
        end        
        if n>1
            var = m2 / (n-1)
        else 
            var=zero(eltype(num))
        end
        return new(meanratio,sn,sd,n,m2,var)
    end
    function CustomVariance()
        return new(0.0,0.0,0.0,0,0.0,0.0)
    end
end

function removeSinglePoint!(a::CustomVariance,num,denom)
    sn=a.sn-num
    sd=a.sd-denom
    meanratio=sn/sd
    #merge: m2 = a.m2 + b.m2 + 2*(a.sn*(a.meanratio-meanratio)+§*(num/denom-meanratio)) + a.sd*(meanratio^2-a.meanratio^2) + denom*(meanratio^2-num/denom^2)
    #unmerge:    
    a.m2 = a.m2 + 2*(num*(a.meanratio-num/denom) + sn*(a.meanratio-meanratio)) + denom*((num/denom)^2-a.meanratio^2) + sd*(meanratio^2-a.meanratio^2)
    a.sn=sn
    a.sd=sd
    a.n -= 1
    a.meanratio=meanratio
    a.var=a.m2/(a.n-1)
    return nothing 
end

function addSinglePoint!(a::CustomVariance,num,denom)
    sn=a.sn+num
    sd=a.sd+denom   
    a.n=a.n+1
    meanratio=sn/sd
    a.m2 = a.m2 + 2*(a.sn*(a.meanratio-meanratio)+num*(num/denom-meanratio)) + a.sd*(meanratio^2-a.meanratio^2) + denom*(meanratio^2-(num/denom)^2)
    a.sn=sn
    a.sd=sd
    a.meanratio=meanratio
    a.var=a.m2/(a.n-1)
    return nothing 
end

function Base.merge!(a::CustomVariance,b::CustomVariance)
    #merges b into a
    #i.e. the points form b are 'added' to a
    sn=a.sn+b.sn
    sd=a.sd+b.sd   
    n=a.n+b.n
    meanratio=sn/sd
    m2 = a.m2 + b.m2 + 2*(a.sn*(a.meanratio-meanratio)+b.sn*(b.meanratio-meanratio)) + a.sd*(meanratio^2-a.meanratio^2) + b.sd*(meanratio^2-b.meanratio^2)
    var = m2/(n-1)
    a.sn=sn
    a.sd=sd
    a.n=n
    a.meanratio=meanratio
    a.m2=m2
    a.var=var
    return nothing
end

function Base.merge(a::CustomVariance,b::CustomVariance)
    c=deepcopy(a)
    merge!(c,b)
    return c
end 

"""
function unmerge!(a::CustomVariance,b::CustomVariance)
this function updates the mean and variance of a
The assumption is that the underlying data of b is a subset of a.
If this is not the case, the result of this function is most likely meaningless
That is to say "subtract" the series b from the series a and update the mean and variance of a accordingly
"""
function unmerge!(a::CustomVariance,b::CustomVariance)
    sn=a.sn-b.sn
    sd=a.sd-b.sd   
    n=a.n-b.n
    meanratio=sn/sd
    #merge: m2 = a.m2 + b.m2 + 2*(a.sn*(a.meanratio-meanratio)+b.sn*(b.meanratio-meanratio)) + a.sd*(meanratio^2-a.meanratio^2) + b.sd*(meanratio^2-b.meanratio^2)
    #unmerge:    
    m2 = a.m2 - b.m2 + 2*(b.sn*(a.meanratio-b.meanratio) + sn*(a.meanratio-meanratio)) + b.sd*(b.meanratio^2-a.meanratio^2) + sd*(meanratio^2-a.meanratio^2)
    var = m2/(n-1)
    a.sn=sn
    a.sd=sd
    a.n=n
    a.meanratio=meanratio
    a.m2=m2
    a.var=var    
    return nothing
end