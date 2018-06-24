@testset "Helper functions Testing" begin

@test isContiguous(collect(1:23))
@test isContiguous(collect(22:23))
@test isContiguous(collect(23:23))
@test !isContiguous([1,9,3])
@test !isContiguous([1,3,2])

end