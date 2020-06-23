@info("DTM: running build.jl")

try
    import Pkg 
    println("")
    println("Building PyCall ...")
    ENV["PYTHON"]="" 
    Pkg.build("Conda") 
    Pkg.build("PyCall")
    using PyCall 

    global const pyModPandas = PyNULL()
    global const pyModxlsxwriter = PyNULL()
    
    copy!(pyModPandas, pyimport_conda("pandas","pandas"))
    copy!(pyModxlsxwriter, pyimport_conda("xlsxwriter","xlsxwriter"));
    writer=pyModPandas.ExcelWriter("blabla_deletme.xlsx", engine = "xlsxwriter")
    @show typeof(writer)
    println("If this worked, then pandas and xlsxwriter should be installed propoerly");
    println("Done building PyCall.")
catch e
    @warn("DTM: Build script failed.")
    println("DTM: Build script failed.")
    @show e
end
exit()

@info("DTM: build script finished.")