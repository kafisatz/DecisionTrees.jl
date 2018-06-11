
x=[2,3]
y=convert(Vector{Union{Missing,Int64}},x)
typeof(y)

"""
tries to convert x::Vector{Union{Missing,T}} to type Vector{T}
"""
function tryToRemoveMissingFromType(x::AbstractVector)
    @show elt=eltype(x)
    if typeof(elt)!=Union
        return nothing
    end
    a=elt.a
    b=elt.b
    tragetT = ifelse(a==Missing,b,a)
    try
        res=convert(Vector{tragetT},x) 
        return res
    catch
       @warn("oops")
    end
    return nothing
end


z=tryToRemoveMissingFromType(y)
@show eltype(z),eltype(x),eltype(y)
