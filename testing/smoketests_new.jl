include("C:\\Users\\bernhard.konig\\Documents\\ASync\\publicnl\\Projects\\Active Projects\\China Guanjun Lead\\Working Files\\Julia Code\\analyze_pw_only_data_init_only_for_CV.jl")
datafolder="C:\\temp\\"
#peform CV for boosting
niter=3
mf=0.015
randomw=0.0
minw=-0.01
subsampling_features_prop=0.66666
subsampling_prop=1.0
smoothEstimates=1
str_core=string("boost_it",niter); #"_mf",mf,"_rand",randomw,"_minw",minw,"_subsp",subsampling_features_prop,"_subsampl",subsampling_prop,"_smooth",smoothEstimates);
DTM.updateSettingsMod!(this_sett,write_statistics=false,minw=minw,randomw=randomw,mf=mf,niter=niter,subsampling_features_prop=subsampling_features_prop,subsampling_prop=subsampling_prop,smoothEstimates=smoothEstimates,boolRandomizeOnlySplitAtTopNode="true",model_type="boosted_tree",scorebandsstartingpoints="1,100,200,300,400,500,600,700,800,900",showTimeUsedByEachIteration="false",moderationvector="$(mf)",write_result=false,print_details="false",boolSaveJLDFile="false",bool_write_tree="false",write_iteration_matrix="false",write_sas_code="false");
outfilename=string(datafolder,"tmp_del");dataFilename=string(datafolder,"prepped.jld");
outfileStringOnly=string("tree_13pct_sample_x");outfilename=string(datafolder,outfileStringOnly);
key,trn,numeratortrn,denominatortrn,weighttrn,trn_numfeatures,trn_charfeatures_PDA,keyval,val,val_numfeatures,val_charfeatures_PDA,numeratorval,denominatorval,weightval,mappings,sett,num_levels,char_levels,all_levels,all_levels_as_string_vector,names_and_levels,candMatWOMaxValues=DTM.prep_data_from_df_modernized(df_prepped,this_sett,datafolder,outfilename,dataFilename,outfileStringOnly);
@time result_of_runmodel,resulting_model=DTM.run_model_actual(key,trn,numeratortrn,denominatortrn,weighttrn,trn_numfeatures,trn_charfeatures_PDA,keyval,val,val_numfeatures,val_charfeatures_PDA,numeratorval,denominatorval,weightval,mappings,datafolder,outfilename,outfileStringOnly,DTM.global_const_shift_cols,sett,num_levels,char_levels,all_levels,all_levels_as_string_vector,names_and_levels,candMatWOMaxValues,dataFilename);
using JLD,HDF5,Base.Test
#save(string(dictGlobalSettings["Data Folder"],"\\expected_china_b.jld"),"exceldata",resulting_model.exceldata)

di=load(string(dictGlobalSettings["Data Folder"],"\\expected_china_b.jld"));
expect=di["exceldata"];
@test isequal(expect.sheets[3].data,resulting_model.exceldata.sheets[3].data)
@test isequal(expect.sheets[4].data,resulting_model.exceldata.sheets[4].data)


unshift!(LOAD_PATH,"C:\\Users\\bernhard.konig\\Documents\\ASync\\irobot\\TariffWatch New Repo\\algorithms\\Julia\\Code\\DTM2_dev\\src")
using DTM2

@time result_of_runmodel,resmdldtm2=DTM2.run_model_actual(key,trn,numeratortrn,denominatortrn,weighttrn,trn_numfeatures,trn_charfeatures_PDA,keyval,val,val_numfeatures,val_charfeatures_PDA,numeratorval,denominatorval,weightval,mappings,datafolder,outfilename,outfileStringOnly,DTM.global_const_shift_cols,sett,num_levels,char_levels,all_levels,all_levels_as_string_vector,names_and_levels,candMatWOMaxValues,dataFilename);

@test isequal(expect.sheets[3].data,resmdldtm2.exceldata.sheets[3].data)
@test isequal(expect.sheets[4].data,resmdldtm2.exceldata.sheets[4].data)