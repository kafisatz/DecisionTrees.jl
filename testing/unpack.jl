function myunzip(file,directory)
    @assert splitext(file)[2]==".zip"
    if Sys.isunix() && Sys.KERNEL != :FreeBSD
        return (`unzip -x $file -d $directory`)
    end

    if Sys.KERNEL == :FreeBSD
        tar_args = ["--no-same-owner", "--no-same-permissions"]
        return pipeline(
            `/bin/mkdir -p $dir`,
            `/usr/bin/tar -xf $file -C $dir $tar_args`)
    end
    
    if Sys.iswindows()
        #exe7z = joinpath(Sys.BINDIR, "7z.exe")
        #return (`$exe7z x $file -y -o$directory`)
        cmd = `powershell.exe -nologo -noprofile -command "& { Add-Type -A 'System.IO.Compression.FileSystem'; [IO.Compression.ZipFile]::ExtractToDirectory('$(file)', '$(directory)'); }"`
        return cmd        
    end
end #end function


#=
    fi="H:\\Code\\Julia\\DecisionTrees.jl\\data\\freMTPL2\\freMTPL2.zip"
    cmd=myunzip(fi,splitdir(fi)[1])
    run(cmd)
=#