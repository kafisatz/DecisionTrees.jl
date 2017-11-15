#this_dir=string(dictGlobalSettings["Julia Code Folder"],"\\src");cd(this_dir);push!(LOAD_PATH,this_dir)
this_dir="C:\\Users\\bernhard.konig\\Documents\\ASync\\irobot\\TariffWatch New Repo\\algorithms\\Julia\\Code\\dev\\src";cd(this_dir);push!(LOAD_PATH,this_dir)
#include(string(this_dir,"\\DTM.jl"))
using DTM
using DataFrames

#=
reffile="A:\\JuliaTestData\\bagging_small.csv"
reffile="C:\\Users\\bernhard.konig\\Documents\\ASync\\irobot\\data_tmp\\bagging_small.csv"
eltypev=DTM.define_eltypevector(DataFrames.readtable(reffile,nrows=100),DTM.global_const_shift_cols,35)
reference_db=readtable(reffile,eltypes=eltypev,nrows=-1,normalizenames=true)
=#

tmpf="A:\\CH\\03Testing October 2017\\motordata_10k.csv"
tmpf="C:\\Users\\bernhard.konig\\Documents\\ASync\\irobot\\data_tmp\\motordata_200k.csv"
df_tmp=readtable(tmpf);

expl=["AFF_MM" "FLG_BERSANI" "f_anz_cli_Rca" "f_citt" "anno" "f_cod_fraz" "f_flg_rinnovo_new" "f_num_gar" "F_rin_cons" "f_num_pol_att" "anno_pat" "fonte" "FLG_CNT_INT" "f_can_pol" "f_profes" "f_flg_cla_tpgu" "f_sex" "f_eta" "f_anz_cli_Zuco" "f_etavei" "f_alim" "F_stat_civ_gran" "f_pagamento_new" "f_att_na_old" "f_att_pulito" "f_cod_marc" "f_ant_rit" "f_reg" "COND_ABI" "UT_AUTO" "KM_AN" "SG_PR" "numpol_att_motor" "Ass" "Rit" "inc_fur" "tut" "Nat" "Soc" "ktot" "Crist" "Inf" "bluk" "eta_conseg_patente" "anni_di_patente"][:]

#prep df
warn("add option forced_categorical_variables to specify columns which are interpreted as categorical variables even if they are numerical (e.g. [1 2 4 9 1][:] ...)")
@time (df_prepped,this_sett)=prepare_dataframe_for_dtm(df_tmp,numcol="n_claims_artificial",denomcol="exposure",independent_vars=expl);

#remove all rows which have some NA entries
#an alternative is of course to replace the NAs with some specific value
origrows=size(df_prepped,1)
for i=1:size(df_prepped,2)
    df_prepped = df_prepped[.!isna.(df_prepped[i]),:]
end
@show droppedrows=origrows-size(df_prepped,1)
@show size(df_prepped)

@time updateSettingsMod!(this_sett,model_type="boosted_tree",minw=-.15,randomw=.1,mf=.1,subsampling_features_prop=1,niter=50,subsampling_prop=1.00)

outfileStringOnly="boosti"
datafolder="C:\\temp\\"
outfilename=string("c:\\temp\\",outfileStringOnly)
dataFilename="C:\\temp\\where_to_save_jld_file"
key,trn,numeratortrn,denominatortrn,weighttrn,trn_numfeatures,trn_charfeatures_PDA,keyval,val,val_numfeatures,val_charfeatures_PDA,numeratorval,denominatorval,weightval,mappings,sett,num_levels,char_levels,all_levels,all_levels_as_string_vector,names_and_levels,candMatWOMaxValues=prep_data_from_df_modernized(df_prepped,this_sett,datafolder,outfilename,dataFilename,outfileStringOnly)

tic()
result_of_runmodel,resulting_model=run_model_actual(key,trn,numeratortrn,denominatortrn,weighttrn,trn_numfeatures,trn_charfeatures_PDA,keyval,val,val_numfeatures,val_charfeatures_PDA,numeratorval,denominatorval,weightval,mappings,datafolder,outfilename,outfileStringOnly,DTM.global_const_shift_cols,sett,num_levels,char_levels,all_levels,all_levels_as_string_vector,names_and_levels,candMatWOMaxValues,dataFilename)

#a
