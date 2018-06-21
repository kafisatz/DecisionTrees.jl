
#function lrIndices(trnidx::Vector{Int},f::PooledArray{String,T,1},subset::Array) where T<:Unsigned

trn=[1,2,9]
f=PooledArray(["A","B","A","AB","A","B","A","AB","ASX"])
sbset=UInt8[1,3]