import StatsBase
n = 300000

mydata = rand(n)

intIdx = StatsBase.sample(1:n, fld(n, 2), replace=false)
sort!(intIdx)
complIdx = Vector{eltype(intIdx)}(undef, 0)
for i = 1:n
    found = searchsorted(intIdx, i)
    if length(found) < 1
        push!(complIdx, i)
    end
end

function sumlr(dat, l, r)
    sl = sr = zero(eltype(dat))

    for i in l
        sl += dat[i]
    end
    for i in r
        sr += dat[i]
    end
    sl, sr
end

u32a = convert(Vector{UInt32}, intIdx)
u32b = convert(Vector{UInt32}, complIdx)
@benchmark sumlr(mydata, u32a, u32b)
@benchmark sumlr(mydata, intIdx, complIdx)


