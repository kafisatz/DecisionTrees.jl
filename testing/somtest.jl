struct Mi{T<:Unsigned}
    ii::Int
    m::Int
end

Mi{UInt8}(3,2)


h=Mi(3,Mi(2,3))