using Coverage
import DataFrames
import DelimitedFiles

thisFolder = Pkg.dir("DecisionTrees")
thisFolder = joinpath(thisFolder, "..", "memHist", "2018 06 24c")
isdir(thisFolder)
mr = analyze_malloc(thisFolder)
mr2 = map(x->[x.bytes,x.filename,x.linenumber], mr);
mr3 = deepcopy(mr2[1])
for i = 2:length(mr2)
	mr3 = hcat(mr3, mr2[i])
end
mr4 = permutedims(mr3, [2,1]);
# mr4=readcsv("C:\\temp\\alloc.csv")


function list_of_files(;fldr=pwd(),ext=".jl")
    li = Vector{Vector{String}}()
    for (root, dirs, files) in walkdir(fldr)
        # println("Directories in $root")
        # for dir in dirs
        #   println(joinpath(root, dir)) # path to directories
        # end
        # println("Files in $root")
        for file in files
            thisf = joinpath(root, file)
            if !occursin(".git", thisf) && (lowercase(thisf[end - length(ext) + 1:end]) == lowercase(ext))
                # println(thisf) # path to files
                fn_only = splitdir(thisf)[2]
                wo_fldr = replace(thisf, fldr => "")
                if wo_fldr[1] == '\\' || wo_fldr[1] == '/'
                    wo_fldr = wo_fldr[2:end]
                end
                push!(li, [thisf,wo_fldr,fn_only])
            end
        end
    end
    return li
end

# add new column
newcolumn = similar(mr4[:,1])
newcolumn2 = similar(mr4[:,1])
mr5 = hcat(mr4, newcolumn, newcolumn2)
mr5[:,4] .= "file not found"
mr5[:,5] .= "code not found"

# lookup filenames of *.mem
filelist = unique(mr4[:,2])
filelist2 = deepcopy(filelist)
for i = 1:length(filelist)
    fi = filelist[i][1:end - 4]
    fi2 = splitdir(fi)[2]
    filelist2[i] = fi2
    # @show isfile(fi)
end

# generate list of all *.jl files
file_list_of_julia_code = list_of_files(fldr=Pkg.dir("DecisionTrees"))
file_list_of_julia_code2 = map(x->x[3], file_list_of_julia_code)

filelocations = Dict{String,String}()
for ii = 1:size(filelist2, 1)
    fi = filelist2[ii]
    fiMem = filelist[ii]
    idx = indexin([fi], file_list_of_julia_code2)
    if idx[1] == nothing
        @warn("Oops: this one was not found")
        @show fi
    else 
        # @show file_list_of_julia_code[idx[1]][1]
        # @show fi
        filelocations[fiMem] = file_list_of_julia_code[idx[1]][1]
    end
end

# add location
for i = 1:size(mr5, 1)
    if haskey(filelocations, mr5[i,2])
        mr5[i,4] = filelocations[mr5[i,2]]
    end 
end

codeDict = Dict{String,Vector{String}}()
for fi in values(filelocations)
    if isfile(fi)
        codeDict[fi] = readlines(fi)
    end
end


for memline = 1:size(mr5, 1)
    # memline=2118
    fi = mr5[memline,4]
    lineNr = mr5[memline,3]
    # read the file 
    if haskey(codeDict, fi)
        rl = codeDict[fi]
        if length(rl) >= lineNr
            codeOfThatLine = rl[lineNr]
            mr5[memline,5] = codeOfThatLine
        end
    end
end

fnOnly = map(x->x[2], splitdir.(mr5[:,4]))
mr6 = hcat(mr5, fnOnly)

# sort data 
srtp = sortperm(convert(Vector{Int}, mr6[:,1]), rev=true)

# define header
hdr = ["Allocations" "MemFile" "LineNumber" "CodeFile" "CodeInThatLine" "FileName Only"]
mr7 = vcat(hdr, mr6)

# change columns 
tmp = mr7[:,4]
mr7[:,4] = mr7[:,3]
mr7[:,3] = tmp

DelimitedFiles.writedlm("C:\\temp\\allocations.csv",mr7,',');

#= 
include(joinpath(Pkg.dir("DecisionTrees"),"testing","allocation_coverage_analysis.jl")) =#