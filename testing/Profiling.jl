using Revise
using Profile
using JLD2
using CSV,DataFrames
using DecisionTrees  

include(joinpath("..","test\\runtests.jl"))

datafile=string("data\\freMTPL2\\freMTPL2.csv")
@assert isfile(datafile);
@time fullData=CSV.read(datafile,rows_for_type_detect=100000,allowmissing=:none,categorical=false);

areasSorted=sort(unique(fullData[:Area]))
AreaInteger=map(x->findall((in)([x]),areasSorted)[1],fullData[:Area])
fullData[:AreaInteger]=AreaInteger

#correct for unreasonable observations
for i=1:size(fullData,1)
    fullData[:ClaimNb][i]=min(4,fullData[:ClaimNb][i])
    fullData[:Exposure][i]=min(1.0,fullData[:Exposure][i])
end
    
#set independent variables
selected_explanatory_vars=["Area","AreaInteger","VehPower","VehAge","DrivAge","BonusMalus","VehBrand","VehGas","Density","Region"]

dtmtable,sett,dfprepped=prepare_dataframe_for_dtm!(fullData,trnvalcol="trnTest",numcol="ClaimNb",denomcol="Exposure",weightcol="Exposure",independent_vars=selected_explanatory_vars);

updateSettingsMod!(sett,minw=-0.03,model_type="boosted_tree",niter=40,mf=0.025,subsampling_features_prop=.7,boolCalculatePoissonError=true)

resultingFiles,resM=dtm(dtmtable,sett)

@profile dtm(dtmtable,sett)


li, lidict = Profile.retrieve()

profile_loc="R:\\temp\\profile.jlprof"
dd=splitdir(profile_loc)[1]
if isdir(dd)
    @warn("This is currently failing with JLD2. I do no want to install JLD as it causes some compatibility issues (for now)")
#    JLD2.@save "mi.jlddd" li lidict
else
    @warn("Folder $(dd) not found!")
end

#consider
Profile.print(noisefloor=2.0,sortedby=:count)

#open("C:\\temp\\logs\\profile.txt", "w") do s
#    Profile.print(IOContext(s, :displaysize => (24, 500)))
#end

#=
#to open and read the data
using HDF5, JLD2, ProfileView
JLD2.@load "/tmp/profdata.jld"
ProfileView.view(li, lidict=lidict)
=#