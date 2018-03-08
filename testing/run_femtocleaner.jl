

Pkg.clone("https://github.com/Keno/AbstractTrees.jl")
Pkg.clone("https://github.com/JuliaComputing/Deprecations.jl")
Pkg.clone("https://github.com/JuliaComputing/FemtoCleaner.jl")

using FemtoCleaner
FemtoCleaner.cleanrepo("C:\\temp\\ori\\DecisionTrees.jl"; show_diff = true, delete_local = false)
