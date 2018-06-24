@testset "Helper functions Testing" begin

@test DecisionTrees.isContiguous(collect(1:23))
@test DecisionTrees.isContiguous(collect(22:23))
@test DecisionTrees.isContiguous(collect(23:23))
@test !DecisionTrees.isContiguous([1,9,3])
@test !DecisionTrees.isContiguous([1,3,2])

end