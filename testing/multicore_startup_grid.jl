number_of_threads=8 #needs to be equal to the number of cores (or processes which will "work")
addprocs(number_of_threads);
@everywhere this_dir=string(dictGlobalSettings["Julia Code Folder"],"\\src");
@everywhere push!(LOAD_PATH,this_dir)
using DTM

#=
	chosenj="boosting_tiny";
	ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".jld") string("out_",chosenj)];
	resbool=DTM.run_model(ARGS)
=#

minw_list=[70 75 80 90]
randomw_list=[0.0]
subsampling_features_prop_list=[0.2 0.5]
subsampling_prop_list=[0.5 1.0]
niter_list=[700]
mf_list=[0.14 .15 .16]
nscores_list=[2000]
boolRandomizeOnlySplitAtTopNode_list=[0]
smoothEstimates_list=[0]

loc=string(dictGlobalSettings["Data Folder"],"\\temp\\azsui_vK\\");
dta=string(loc,"data.jld")
thesesett=string(loc,"sett1.csv");
ARGS=[thesesett dta "outres"];

rxx=DTM.grid_search(ARGS,minw_list,randomw_list,subsampling_features_prop_list,smoothEstimates_list,niter_list,mf_list,nscores_list,boolRandomizeOnlySplitAtTopNode_list,subsampling_prop_list)
