
using DecisionTrees
using Base.Test


#peform CV for boosting
niter=50 #5 #125
mf=0.02
randomw=0.0
minw=-0.08
subsampling_features_prop=0.8
subsampling_prop=1.0
smoothEstimates=1

updateSettingsMod!(this_sett,minw=minw,write_dot_graph=true,write_csharp_code=true,write_vba_code=true,randomw=randomw,mf=mf,niter=niter,subsampling_features_prop=subsampling_features_prop,subsampling_prop=subsampling_prop,smoothEstimates=smoothEstimates,boolRandomizeOnlySplitAtTopNode="true",model_type="build_tree",scorebandsstartingpoints="1,100,200,300,400,500,600,700,800,900",showTimeUsedByEachIteration="false",moderationvector="$(mf)",write_result=false,print_details="true",boolSaveResultAsJLDFile="false",bool_write_tree="true",write_iteration_matrix="false",write_sas_code="true",graphvizexecutable="c:\\Program Files (x86)\\Graphviz2.38\\bin\\dot.exe");
datafolder="C:\\temp\\"
a,b=run_model_actual(dtmtable,this_sett,string(datafolder,"\\tree_",minw,".txt"));

st2=deepcopy(this_sett);
st2.model_type="build_tree"
st2.minw=300
this_sett.minw=800
a,b=run_model_actual(dtmtable,[this_sett,st2 ],string(datafolder,"\\tree_",minw,".txt"));

#=

#outputfldr="C:\\temp\\joutput\\"

dta_fld=dictGlobalSettings["Data Folder"]
run_legacy("tree_normal",dta_fld)

mdls=["tree_normal"  "tree_single_leaf" "tree_small" "boosting_small" "boosting_tiny" "boosting_tiny_subs" "boosting_tiny"][:]
#"bagging_small"
for m in mdls
    @test eltype(run_legacy(m,dta_fld))==String
end
#dictGlobalSettings["Julia Code Folder"]="R:\\TariffWatch New Repo\\algorithms\\Julia\\Code\\dev\\"

chosenj="tree_small" #this runs a small version to compile all functions
ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)]
DTM.run_model(ARGS)
#using DTM # reload("DTM")

warn("todo, find out why this does not work......!")
chosenj="boosting_minimalexample";ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];resbool=DTM.run_model(ARGS)

#multirow settings
chosenj="boosting_minimalexample_multirow_settings";ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];resbool=DTM.run_model(ARGS)
chosenj="tree_multirow";ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];resbool=DTM.run_model(ARGS)
#tiny tree
chosenj="tree_minimalexample";ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];resbool=DTM.run_model(ARGS)
#tree which results in a single leaf
chosenj="tree_single_leaf";ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];resbool=DTM.run_model(ARGS)
#boosting which results in single leafs
chosenj="boosting_single_leafs";ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];resbool=DTM.run_model(ARGS)
chosenj="boostmed_probably_generates_single_leaf_iteration";ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];resbool=DTM.run_model(ARGS)
#other trees
chosenj="tree_normal";ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];resbool=DTM.run_model(ARGS)
chosenj="boosting_minimalexample";ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];resbool=DTM.run_model(ARGS)
chosenj="bagging_minimalexample";ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];resbool=DTM.run_model(ARGS)
chosenj="tree_tiny";ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];resbool=DTM.run_model(ARGS)
#tiny boosting

chosenj="boosting_tiny";ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];resbool=DTM.run_model(ARGS)
#tiny bagging
chosenj="bagging_tiny";ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];resbool=DTM.run_model(ARGS)
#boosting variations
chosenj="boost_tiny_subs";ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];@time resbool=DTM.run_model(ARGS)
chosenj="boost_tiny_subs_withrepl";ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];@time resbool=DTM.run_model(ARGS)
chosenj="boosting_minimalexample";ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];@time resbool=DTM.run_model(ARGS)
chosenj="boostmed_1iter";ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];@time resbool=DTM.run_model(ARGS)
chosenj="boostmed";ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];@time resbool=DTM.run_model(ARGS)
chosenj="boostmed2";ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];resbool=DTM.run_model(ARGS)
chosenj="boostmed_smoothed";ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];@time resbool=DTM.run_model(ARGS)
chosenj="boost_az_mob_test";ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];@time resbool=DTM.run_model(ARGS)
#subsampling
subsamplingtests=["boost_tiny_subs_withrepl" "boost_tiny_subs_withrepl_onlytopn" "boost_tiny_subs_worepl_allnodes" "boost_tiny_subs_worepl_onlytopn"]
for chosenj in subsamplingtests
	ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];@time resbool=DTM.run_model(ARGS)
end

baggingtests=["bagging_minimalexample" "bagging_minimalexample2" "bagging_minimalexample3"]
for chosenj in baggingtests
	ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];@time resbool=DTM.run_model(ARGS)
end

chosenj="bag_with_repl_200k";ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];@time resbool=DTM.run_model(ARGS)
chosenj="bag200k";ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];@time resbool=DTM.run_model(ARGS)

ARGS=["e:\\temp\\boominw75iters200_TKsettings.csv" "e:\\temp\\boominw75iters200_TK.csv" "out_boominw75iters200_TK"]
@time resbool=DTM.run_model(ARGS)



#boosting
chosenj="boosting_small"
ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)]
@time resbool=DTM.run_model(ARGS) #with julia startup 41s, time for rist execution of run_model: 32s, time for next execution of the same model 5s (as code is already compiled)
# @time resbool=DTM.run_model(ARGS) #time 5.4s

chosenj="julia200k";ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];@time resbool=DTM.run_model(ARGS)
chosenj="julia200k_200s";ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];@time resbool=DTM.run_model(ARGS)
chosenj="julia200k_10s";ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];@time resbool=DTM.run_model(ARGS)
#chosenj="julia200k_100s";ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];@time resbool=DTM.run_model(ARGS)
chosenj="julia200k_boosting_subsampling";ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];@time resbool=DTM.run_model(ARGS)

chosenj="boosting_testmodel_small";ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];@time resbool=DTM.run_model(ARGS)

chosenj="boosting_minimalexample";ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];resbool=DTM.run_model(ARGS)


#bagging
chosenj="bagging_small"
info(string("\nrunning ",chosenj,"\n\n"))
ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)]
resbool=DTM.run_model(ARGS)

#boosting deterministics
chosenj="boosting_small_deterministic"
info(string("\nrunning ",chosenj,"\n\n"))
ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)]
resbool=DTM.run_model(ARGS) # run once to trigger compilation



chosenj="julia1m_deterministic"
info(string("\nrunning ",chosenj,"\n\n"))
ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)]
@time resbool=DTM.run_model(ARGS) # run once to trigger compilation

chosenj="uniqa5m"
ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings2.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".jld") string("out_",chosenj)]
ARGS=[string("T:\\temp\\",chosenj,".settings2.csv") string("T:\\temp\\",chosenj,".jld") string("out_",chosenj)]
@time resbool=DTM.run_model(ARGS) #with julia startup 41s, time for rist execution of run_model: 32s, time for next execution of the same model 5s (as code is already compiled)
#the modelling for the 5m uniqa data set uses about 6GM of RAM.


chosenj="uniqa10m"
ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings2.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".jld") string("out_",chosenj)]
ARGS=[string("T:\\temp\\",chosenj,".settings.csv") string("T:\\temp\\",chosenj,".jld") string("out_",chosenj)]
ARGS=[string("T:\\temp\\",chosenj,".settings.csv") string("T:\\temp\\",chosenj,".csv") string("out_",chosenj)]
@time resbool=DTM.run_model(ARGS) #with julia startup 41s, time for rist execution of run_model: 32s, time for next execution of the same model 5s (as code is already compiled)
#the modelling for the 5m uniqa data set uses about 6GM of RAM.

ARGS=["a:\\temp\\random.settings.csv" "a:\\temp\\random.jld" "outmxe"];
resbool=DTM.run_model(ARGS);

ARGS=["a:\\temp\\boost1settings.csv" "a:\\temp\\boost1.csv" "outmxe"];
resbool=DTM.run_model(ARGS);
=#