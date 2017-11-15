this_dir="C:\\Users\\bernhard.konig\\Documents\\ASync\\irobot\\TariffWatch New Repo\\algorithms\\Julia\\Code\\dev\\src";cd(this_dir);push!(LOAD_PATH,this_dir)
using DTM

# chosenj="boosting_minimalexample_multirow_settings";ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];resbool=DTM.run_model(ARGS)



outfileStringOnly="boostedtree_2" #400s for minw 10% and 50 iters (->8s per iter)
updateSettingsMod!(this_sett,model_type="boosted_tree",minw=-.2,randomw=0.05,mf=.01,niter=2)
@time result_of_runmodel,resulting_model=run_model_actual(key,trn,numeratortrn,denominatortrn,weighttrn,trn_numfeatures,trn_charfeatures_PDA,keyval,val,val_numfeatures,val_charfeatures_PDA,numeratorval,denominatorval,weightval,mappings,datafolder,outfilename,outfileStringOnly,DTM.global_const_shift_cols,sett,num_levels,char_levels,all_levels,all_levels_as_string_vector,names_and_levels,candMatWOMaxValues,dataFilename);

global const defaulted_modelstats_df=deepcopy(resulting_model.modelstats)


df_modelstats_loc="R:\\TariffWatch New Repo\\algorithms\\Julia\\Code\\defaulted_modelstats_df.jld"
JLD.save(df_modelstats_loc,"defaulted_modelstats_df",defaulted_modelstats_df)


loaded_dict = JLD.load(df_modelstats_loc)
global const defaulted_modelstats_df = loaded_dict["defaulted_modelstats_df"]