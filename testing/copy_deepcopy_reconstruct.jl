using BenchmarkTools 

function reconstruct(a)
    b = similar(a)
    for i in eachindex(a)
        @inbounds b[i] = a[i]
    end
    return b
end


nn = 200000
trn = rand(Int, nn)

@btime deepcopy($trn)
@btime copy($trn)
@btime reconstruct($trn)
@assert copy(trn) == deepcopy(trn) == reconstruct(trn)


nn = 200000
xrand = rand()
trn = collect(ceil(Int, xrand / 2) * nn:ceil(Int, (1 - xrand) / 20) * nn)

@btime deepcopy($trn)
@btime copy($trn)
@btime reconstruct($trn)
@assert copy(trn) == deepcopy(trn) == reconstruct(trn)