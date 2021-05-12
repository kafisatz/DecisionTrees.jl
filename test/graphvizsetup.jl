    global graphvizsetup_is_working = false
    @info("Testing Graphviz Setup")
    @testset verbose = true "GraphViz Setup" begin 
        if haskey(ENV,"graphvizdot")
            #if key exists this takes precedence
            global graphvizexe = ENV["graphvizdot"]
            #=
                ENV["graphvizdot"]="c:\\program files\\Graphviz\\bin\\dot.exe"
            =#
        else 
            global graphvizexe = "dot" #linux & appleOS assume it is in in path
        end
        
        gcmd = `$(graphvizexe) "-V"`
        try 
            gres = run(gcmd)
            @test gres.exitcode == 0
            @show gres.exitcode
            global graphvizsetup_is_working = true
            @test true #graphviz setup works
        catch e
            if Sys.iswindows()
                try
                    global graphvizexe = "c:\\program files\\Graphviz\\bin\\dot.exe" 
                    gcmd = `$(graphvizexe) "-V"`
                    gres = run(gcmd)
                    @test gres.exitcode == 0
                    @show gres.exitcode
                    global graphvizsetup_is_working = true
                    @test true #graphviz setup works
                catch eiwn
                    @warn("DTM: unable to invoke gcmd")
                    @show gcmd
                    print(ewin)
                    @test false        
                end
            else 
                @warn("DTM: unable to invoke gcmd")
                @show gcmd
                print(e)
                @test false
            end
        end
    end