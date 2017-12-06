using DataArrays
using PooledArrays
n=33
x=ceil.(4*rand(n))
pa=PooledArray(x)

pda=compact(PooledDataArray(x))

pda.pool==pa.pool
pda.refs==pa.refs
@assert isequal(pda.pool,pa.pool)
@assert isequal(pda.refs,pa.refs)
@assert isequal(pda,pa)
@assert pda===pa
pda===pa


a=Array{PooledArray{String,UInt8}}(2)
