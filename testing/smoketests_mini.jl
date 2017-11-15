#=
this_dir="C:\\Users\\bernhard.konig\\Documents\\ASync\\irobot\\TariffWatch New Repo\\algorithms\\Julia\\Code\\dev\\src";cd(this_dir);push!(LOAD_PATH,this_dir)
using DTM
this_dir="C:\\Users\\bernhard.konig\\Documents\\ASync\\irobot\\TariffWatch New Repo\\algorithms\\Julia\\Code\\dev\\src";cd(this_dir);push!(LOAD_PATH,this_dir)
 dictGlobalSettings["Data Folder"]="C:\\Users\\bernhard.konig\\Documents\\ASync\\irobot\\data_tmp\\"
 chosenj="boosting_small"
 ARGS=[string(dictGlobalSettings["Data Folder"],chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],chosenj,".csv") string("out_",chosenj)]
 @time resbool=DTM.run_model(ARGS) #with julia startup 41s, time for rist execution of run_model: 32s, time for next execution of the same model 5s (as code is already compiled)
 # @time resbool=DTM.run_model(ARGS) #time 5.4s

 chosenj="tree_normal"
 ARGS=[string(dictGlobalSettings["Data Folder"],chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],chosenj,".csv") string("out_",chosenj)]
 @time resbool=DTM.run_model(ARGS) #with julia startup 41s, time for rist execution of run_model: 32s, time for next execution of the same model 5s (as code is already compiled)
 # @time resbool=DTM.run_model(ARGS) #time 5.4s

 =#

this_dir=string(dictGlobalSettings["Julia Code Folder"],"\\src");cd(this_dir);unshift!(LOAD_PATH,this_dir)
using DTM
using Base.Test
#outputfldr="C:\\temp\\joutput\\"

dta_fld=dictGlobalSettings["Data Folder"]
run_legacy("tree_normal",dta_fld)

mdls=["tree_normal"  "tree_single_leaf" "tree_small" "boosting_small" "boosting_tiny" "boosting_tiny_subs" "boosting_tiny"][:]
#"bagging_small"
for m in mdls
    @test eltype(run_legacy(m,dta_fld))==String
end
@show 23
#dictGlobalSettings["Julia Code Folder"]="R:\\TariffWatch New Repo\\algorithms\\Julia\\Code\\dev\\"
