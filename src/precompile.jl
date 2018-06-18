#the files were generated with SnoopCompile

precompFiles=readdir("precompile")
for thisf in precompFiles
    @show thisf
    include(joinpath("precompile",thisf))
    _precompile_()
end