

Pkg.clone("https://github.com/Keno/AbstractTrees.jl")
Pkg.clone("https://github.com/JuliaComputing/Deprecations.jl")
Pkg.clone("https://github.com/JuliaComputing/FemtoCleaner.jl")

Pkg.update()

using FemtoCleaner
FemtoCleaner.cleanrepo("C:\\Users\\bernhard.konig\\Documents\\ASync\\home\\Code\\Julia\\DecisionTrees.jl"; show_diff=true, delete_local=false)
