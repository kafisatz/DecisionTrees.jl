using Coverage
import DelimitedFiles

thisFolder="C:\\Users\\bernhard.konig.ROOT_MILLIMAN\\Documents\\async\\home\\code\\julia\\memHistory"
mr=analyze_malloc(thisFolder);
mr2=map(x->[x.bytes,x.filename,x.linenumber],mr);
mr3=deepcopy(mr2[1])
for i=2:length(mr2)
	mr3=hcat(mr3,mr2[i])
end
mr4=permutedims(mr3,[2,1]);

srtp=sortperm(convert(Vector{Int},mr4[:,1]),rev=true)

#txt1=readlines(fi)
#filelist=(x->x[1:end-4] for x in unique(mr4[:,2]))
filelist=unique(mr4[:,2])
filelist2=deepcopy(filelist)
for i=1:length(filelist)
    fi=filelist[i][1:end-4]
    fi2=splitdir(fi)[2]
    filelist2[i]=fi2
    #@show isfile(fi)
end

DelimitedFiles.writedlm("C:\\temp\\alloc.csv",mr4,',');