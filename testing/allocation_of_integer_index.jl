using BenchmarkTools 

nn = 600000

function newIdx(n)
    return zeros(Float64, nn)
end

@btime newIdx(nn)