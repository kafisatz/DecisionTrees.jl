using Coverage

thisFolder="C:\\Users\\bernhard.konig.ROOT_MILLIMAN\\Documents\\async\\home\\code\\julia\\memHistory"
mr=analyze_malloc(thisFolder);
mr2=map(x->[x.bytes,x.filename,x.linenumber],mr);
mr3=deepcopy(mr2[1])
for i=2:length(mr2)
	mr3=hcat(mr3,mr2[i])
end
#mr3=hcat(mr2[1],mr2[2],mr2[3],mr2[4]);
mr4=permutedims(mr3,[2,1]);

srtp=sortperm(convert(Vector{Int},mr4[:,1]),rev=true)

import DelimitedFiles
DelimitedFiles.writedlm("C:\\temp\\alloc.csv",mr4,',');