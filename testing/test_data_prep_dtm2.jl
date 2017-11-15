this_dir="C:\\Users\\bernhard.konig\\Documents\\ASync\\irobot\\TariffWatch New Repo\\algorithms\\Julia\\Code\\";
#this_dir="R:\\TariffWatch New Repo\\algorithms\\Julia\\Code\\dev\\src"
unshift!(LOAD_PATH,string(this_dir,"\\DTM2_dev\\src"))
unshift!(LOAD_PATH,string(this_dir,"\\dev\\src"))
using Compat
using BenchmarkTools
using Revise    
using DTM2
using DTM
using DataFrames,JLD,HDF5
datafolder="C:\\temp\\";outfilename="C:\\temp\\DTM2Model";dataFilename="c:\\temp\\DTM2Model.jld";outfileStringOnly="DTM2Model";fn=string(datafolder,outfileStringOnly)

#=
include("C:\\Users\\bernhard.konig\\Documents\\ASync\\irobot\\TariffWatch New Repo\\algorithms\\Julia\\Code\\DTM2_dev\\testing\\define_df_prepped.jl");
datafolder="C:\\temp\\";outfilename="C:\\temp\\DTM2Model";dataFilename="c:\\temp\\DTM2Model.jld";outfileStringOnly="DTM2Model";fn=string(datafolder,outfileStringOnly)
@time dtmtable=DTM2.prep_data_from_df(df_prepped,this_dt2msett,fn);
save("C:\\Users\\bernhard.konig\\Documents\\ASync\\irobot\\data_tmp\\prepped_dtm2table.jld","df_prepped",df_prepped,"dtmtable",dtmtable,"this_sett",this_sett,"this_dt2msett",this_dt2msett);
=#
di9=load("C:\\Users\\bernhard.konig\\Documents\\ASync\\irobot\\data_tmp\\prepped_dtm2table.jld");dtmtable=di9["dtmtable"];this_dt2msett=di9["this_dt2msett"];df_prepped=di9["df_prepped"];this_sett=di9["this_sett"];

settstr="""model_type="boosted_tree",write_statistics="true",minw=-.051,randomw=0.0,mf=.05,subsampling_features_prop=1,niter=5,subsampling_prop=1.00,write_iteration_matrix="false",graphvizexecutable="C:\\\\Program Files (x86)\\\\Graphviz2.38\\\\bin\\\\dot.exe",write_result=false"""
eval(parse(string("DTM2.updateSettingsMod!(this_dt2msett,",settstr,")")));DTM2.updateSettingsMod!(this_dt2msett,boolProduceEstAndLeafMatrices=false)
@time st2,dtm2model=DTM2.run_model_actual(dtmtable,this_dt2msett,fn);

settstr="""model_type="bagged_tree",write_statistics="true",minw=-.08,randomw=0.0,mf=.051,subsampling_features_prop=.4,niter=40,subsampling_prop=.5,write_iteration_matrix="false",graphvizexecutable="C:\\\\Program Files (x86)\\\\Graphviz2.38\\\\bin\\\\dot.exe",write_result=false"""
eval(parse(string("DTM2.updateSettingsMod!(this_dt2msett,",settstr,")")));DTM2.updateSettingsMod!(this_dt2msett,boolProduceEstAndLeafMatrices=false)
@time st2,dtm2model=DTM2.run_model_actual(dtmtable,this_dt2msett,fn);

settstr="""model_type="boosted_tree",niter=50,print_details="true",write_sas_code="false",bool_write_tree="false",write_statistics="true",write_iteration_matrix="false",minw=-.1,randomw=0.0,mf=.1,subsampling_features_prop=1,subsampling_prop=1.00,graphvizexecutable="C:\\\\Program Files (x86)\\\\Graphviz2.38\\\\bin\\\\dot.exe",write_result=false"""
eval(parse(string("DTM2.updateSettingsMod!(this_dt2msett,",settstr,")"))); DTM2.updateSettingsMod!(this_dt2msett,boolProduceEstAndLeafMatrices="false");
@time st2,dtm2model=DTM2.run_model_actual(dtmtable,this_dt2msett,fn); 

settstr="""model_type="build_tree",moderationvector="-0.1",niter=10,write_sas_code="false",bool_write_tree="false",write_statistics="true",write_iteration_matrix="false",minw=-.01,randomw=0.0,mf=.1,subsampling_features_prop=1,subsampling_prop=1.00,graphvizexecutable="C:\\\\Program Files (x86)\\\\Graphviz2.38\\\\bin\\\\dot.exe",write_result=false"""
eval(parse(string("DTM2.updateSettingsMod!(this_dt2msett,",settstr,")"))); DTM2.updateSettingsMod!(this_dt2msett,boolProduceEstAndLeafMatrices="false");
@time st2,dtm2model=DTM2.run_model_actual(dtmtable,this_dt2msett,fn);



settstr="""model_type="boosted_tree",niter=30,write_sas_code="false",bool_write_tree="false",write_statistics="false",write_iteration_matrix="false",minw=-.05,randomw=0.0,mf=.1,subsampling_features_prop=1,subsampling_prop=1.00,graphvizexecutable="C:\\\\Program Files (x86)\\\\Graphviz2.38\\\\bin\\\\dot.exe",write_result=false"""
eval(parse(string("DTM2.updateSettingsMod!(this_dt2msett,",settstr,")"))); DTM2.updateSettingsMod!(this_dt2msett,boolProduceEstAndLeafMatrices="false");
Profile.clear()          # in case we have any previous profiling data
Profile.init(10^7, 0.01) #n,delay #Profiling settings
@profile  st2,dtm2model=DTM2.run_model_actual(dtmtable,this_dt2msett,fn);
# profile() #this shows profile graphically in Juno
using ProfileView
val=ProfileView.view()




settstr="""model_type="boosted_tree",niter=10,minw=-.01,bool_write_tree="false",write_statistics="false",write_iteration_matrix="false",randomw=0.0,mf=.1,subsampling_features_prop=1,subsampling_prop=1.00,graphvizexecutable="C:\\\\Program Files (x86)\\\\Graphviz2.38\\\\bin\\\\dot.exe",write_result=false"""
eval(parse(string("DTM2.updateSettingsMod!(this_dt2msett,",settstr,")")))
@time st2,dtm2model=DTM2.run_model_actual(dtmtable,this_dt2msett,fn); #2.1s roughly

dt2model_save=deepcopy(dtm2model);

dtmtable.features[2]
s=dtm2model.rootnode.subset
dtmtable.features[2].pool[s]
dtmtable.mappings[2][s]
dtm2model.rootnode.subset
dtm2model.rootnode.featid

#save("C:\\Users\\bernhard.konig\\Documents\\ASync\\irobot\\data_tmp\\prepped_dtm2.jld","df_prepped",df_prepped,"this_dt2msett",this_dt2msett);
dddict=load("C:\\Users\\bernhard.konig\\Documents\\ASync\\irobot\\data_tmp\\prepped_dtm2.jld");
df_prepped=dddict["df_prepped"];this_dt2msett=dddict["this_dt2msett"];

df_prepped_mini=df_prepped[1:20,:];
@time dtmtable=DTM2.prep_data_from_df(df_prepped,this_dt2msett,fn);
minw=-0.1;mf=0.02;niter=2;subsampling_features_prop=1.0;subsampling_prop=1.0;randomw=0.0;smoothEstimates=1
DTM2.updateSettingsMod!(this_dt2msett,write_statistics=false,minw=minw,randomw=randomw,mf=mf,niter=niter,subsampling_features_prop=subsampling_features_prop,subsampling_prop=subsampling_prop,smoothEstimates=smoothEstimates,boolRandomizeOnlySplitAtTopNode="true",model_type="boosted_tree",scorebandsstartingpoints="1,100,200,300,400,500,600,700,800,900",showTimeUsedByEachIteration="false",moderationvector="$(mf)",write_result=false,print_details="false",boolSaveJLDFile="false",bool_write_tree="false",write_iteration_matrix="false",write_sas_code="false");
#save("C:\\Users\\bernhard.konig\\Documents\\ASync\\irobot\\data_tmp\\prepped_dtm2table.jld","dtmtable",dtmtable,"this_dt2msett",this_dt2msett);


this_dir="C:\\Users\\bernhard.konig\\Documents\\ASync\\irobot\\TariffWatch New Repo\\algorithms\\Julia\\Code\\DTM2_dev\\src";
#this_dir="R:\\TariffWatch New Repo\\algorithms\\Julia\\Code\\dev\\src"
unshift!(LOAD_PATH,this_dir);using Revise;using Compat;using BenchmarkTools,JLD,HDF5
using DTM2;datafolder="C:\\temp\\";outfilename="C:\\temp\\tree_13pct_sample_x";dataFilename="c:\\temp\\mira.jld";outfileStringOnly="tree_13pct_sample_x";fn=string(datafolder,outfileStringOnly)
@time dtmtable=DTM2.prep_data_from_df(df_prepped,this_dt2msett,fn);
#C:\\Users\\bernhard.konig\\Documents\\ASync\\publicnl\\Projects\\Active Projects\\\\China Guanjun Lead\\Working Files\\c01_agg_data\\"
di9=load("C:\\Users\\bernhard.konig\\Documents\\ASync\\irobot\\data_tmp\\prepped_dtm2table.jld");dtmtable=di9["dtmtable"];this_dt2msett=di9["this_dt2msett"];

@time resulting_model=DTM2.run_model_actual(dtmtable,this_dt2msett,fn);

#=
    reload("DTM2")
=#

minw=-0.1;mf=0.02;niter=3;subsampling_features_prop=1.0;subsampling_prop=1.0;randomw=0.00;smoothEstimates=1
DTM2.updateSettingsMod!(this_dt2msett,model_type="build_tree",write_statistics=true,minw=minw,randomw=randomw,mf=mf,niter=niter,subsampling_features_prop=subsampling_features_prop,subsampling_prop=subsampling_prop,smoothEstimates=smoothEstimates,boolRandomizeOnlySplitAtTopNode="true",model_type="boosted_tree",scorebandsstartingpoints="1,100,200,300,400,500,600,700,800,900",showTimeUsedByEachIteration="false",moderationvector="$(mf)",write_result=false,print_details="false",boolSaveJLDFile="false",bool_write_tree="false",write_iteration_matrix="false",write_sas_code="false");
@time resulting_model=DTM2.run_model_actual(dtmtable,this_dt2msett,fn);

minw=-0.1;mf=0.02;niter=2;subsampling_features_prop=.6;subsampling_prop=1.0;randomw=0.05;smoothEstimates=1
DTM2.updateSettingsMod!(this_dt2msett,write_statistics=true,minw=minw,randomw=randomw,mf=mf,niter=niter,subsampling_features_prop=subsampling_features_prop,subsampling_prop=subsampling_prop,smoothEstimates=smoothEstimates,boolRandomizeOnlySplitAtTopNode="true",model_type="boosted_tree",scorebandsstartingpoints="1,100,200,300,400,500,600,700,800,900",showTimeUsedByEachIteration="false",moderationvector="$(mf)",write_result=false,print_details="false",boolSaveJLDFile="false",bool_write_tree="false",write_iteration_matrix="false",write_sas_code="false");
@time resulting_model=DTM2.run_model_actual(dtmtable,this_dt2msett,fn);





for niter in niters 
for mf in mfs #[.03 .05] #unclear what about .08?
for randomw in randomws #[0.0 .02]  #preference 0
for minw in minws #[-.075] #preference -.75
for subsampling_features_prop in subsampling_features_props #[1.0] # .7 .5]
for subsampling_prop in subsampling_props #[.6666 1.0] # yet unclear which is the best paramter here.... #preference .6666
for smoothEstimates in smoothEstimatess #[1 0] #we have not checke this yet!		
        try
            outfileStringOnly=string("tree_rand",randomw,"_minw",minw);
            updateSettingsMod!(this_dt2msett,minw=minw,randomw=randomw,mf=mf,niter=niter,subsampling_features_prop=subsampling_features_prop,subsampling_prop=subsampling_prop,smoothEstimates=smoothEstimates,boolRandomizeOnlySplitAtTopNode="true",model_type="build_tree",scorebandsstartingpoints="1,50,100,150,200,300,400,500,600,700,800,850,900,950",showTimeUsedByEachIteration=1,moderationvector="$(mf)");
            @time result_of_runmodel,resulting_model=run_model_actual(key,trn,numeratortrn,denominatortrn,weighttrn,trn_numfeatures,trn_charfeatures_PDA,keyval,val,val_numfeatures,val_charfeatures_PDA,numeratorval,denominatorval,weightval,mappings,datafolder,outfilename,outfileStringOnly,DTM.global_const_shift_cols,sett,num_levels,char_levels,all_levels,all_levels_as_string_vector,names_and_levels,candMatWOMaxValues,dataFilename);
            nothing
        catch
            warn("BK: this model did not work: $(my_count)")
        end
    my_count+=1
end
end
end
end
end
end 
end
    







settstr="""model_type="build_tree",write_iteration_matrix="false",write_statistics="true",minw=-.051,randomw=0.0,mf=.1,subsampling_features_prop=1,niter=6,subsampling_prop=1.00,graphvizexecutable="C:\\\\Program Files (x86)\\\\Graphviz2.38\\\\bin\\\\dot.exe",write_result=false"""
eval(parse(string("DTM2.updateSettingsMod!(this_dt2msett,",settstr,")")));DTM2.updateSettingsMod!(this_dt2msett,boolProduceEstAndLeafMatrices=false)
@time st2,dtm2model=DTM2.run_model_actual(dtmtable,this_dt2msett,fn);
