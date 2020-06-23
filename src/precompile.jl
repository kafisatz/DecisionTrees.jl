# the precompile_*.jl files were generated with SnoopCompile

precompFiles = readdir(joinpath("src", "precompile"))

# the runtests.jl uses other packages such as CSV
# therefore some of the 'precompile' files do not run (one needs to modify them)
# for now we only include a few of the files 

# notably it is unclear if the files which do not belong to this package have any impact at all
only_run_these = ["precompile_DecisionTrees.jl"] # ,"precompile_Base.jl","precompile_StatsBase.jl","precompile_PyCall.jl","precompile_Random.jl","precompile_Core.jl"]

# @show size(precompFiles)
ijxt = 0
for thisf in precompFiles
    if in(thisf, only_run_these)
        # @show thisf;@show ijxt;ijxt+=1
        include(joinpath("precompile", thisf))                
        _precompile_()
    end
end
