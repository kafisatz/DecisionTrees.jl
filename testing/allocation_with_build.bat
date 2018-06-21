cd %HOMEPATH%
cd Documents\async\home\code\julia\DecisionTrees.jl\
julia --track-allocation=user -e "using Pkg;Pkg.build(\"DecisionTrees\"); Pkg.test(\"DecisionTrees\"); include(\"testing\\allocation.jl\");"