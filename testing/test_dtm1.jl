#df_prepped should be defined 'by now' 
@time (df_prepped_do_not_use_me,this_sett)=DTM.prepare_dataframe_for_dtm!(df_tmp,treat_as_categorical_variable=forced_categorical_variables,trnvalcol=chosen_trn_val_col,weightcol=chosen_weight_col,numcol=chosen_num_col,denomcol=chosen_denom_col,independent_vars=expl);
@show isequal(df_prepped_do_not_use_me,df_prepped)

outfileStringOnlyDTM1="DTM1Model"
datafolder="c:\\temp\\"
outfilenameDTM1=string(datafolder,outfileStringOnlyDTM1)
dataFilenameDTM1=string(datafolder,"DTM1Model.jld")
@time key,trn,numeratortrn,denominatortrn,weighttrn,trn_numfeatures,trn_charfeatures_PDA,keyval,val,val_numfeatures,val_charfeatures_PDA,numeratorval,denominatorval,weightval,mappings,sett,num_levels,char_levels,all_levels,all_levels_as_string_vector,names_and_levels,candMatWOMaxValues=DTM.prep_data_from_df_modernized(df_prepped,this_sett,datafolder,outfilenameDTM1,dataFilenameDTM1,outfileStringOnlyDTM1);

eval(parse(string("DTM.updateSettingsMod!(this_sett,",settstr,")")));
@time st1,dtm1model=DTM.run_model_actual(key,trn,numeratortrn,denominatortrn,weighttrn,trn_numfeatures,trn_charfeatures_PDA,keyval,val,val_numfeatures,val_charfeatures_PDA,numeratorval,denominatorval,weightval,mappings,datafolder,outfilename,outfileStringOnlyDTM1,DTM.global_const_shift_cols,sett,num_levels,char_levels,all_levels,all_levels_as_string_vector,names_and_levels,candMatWOMaxValues,dataFilenameDTM1);

fieldnames(dtm1model) 
fieldnames(dtm2model)
dtm1model.exceldata;    
dtm2model.exceldata;
isequal(dtm1model.exceldata,dtm2model.exceldata)