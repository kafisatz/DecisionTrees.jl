using CSV
using DataFrames
using Revise
using Compat
using StatsBase,BenchmarkTools
using DecisionTrees
using HDF5,JLD


#folder and dataset
    rootfolder="C:\\Users\\bernhard.konig\\Documents\\ASync\\publicnl\\Projects\\Active Projects\\"
    #rootfolder="M:\\Non-Life\\Projects\\Active Projects\\"
    #tmpf=string(rootfolder,"China Guanjun Lead\\Working Files\\data\\data_aggfull.csv")
    tmpf=string(rootfolder,"China Guanjun Lead\\Working Files\\data\\data_agg20k.csv")    
    datafolder=string(rootfolder,"\\China Guanjun Lead\\Working Files\\c04\\")
    #datafolder="c:\\temp\\"
    @assert isdir(datafolder)

#repsonse variable and weight
    chosen_denom_col="N_EARNED_PRM"
    chosen_num_col="N_DVL_LOSS"
    chosen_weight_col="N_EARNED_EXPO"
    chosen_trn_val_col="TRAINING_VALIDATION_BOOL"

#explanatory variables
#BK: Oct 19, 2017: inc_year included as per Guanjun's suggestion
# AMOUNT is excluded as it varies by coverage; the tree should rather not use it (otherwise it may use it as a proxy for cov_name which we do not want.)
#expl_cat=vcat(expl_cat,["COV_NAME"]) #excluding cov_name as we do not want it as a predictor
    expl_numeric=["GRP_NEW_NCD_NUMERIC" "NEW_NCD_LAST_NUMERIC" "GRP_CAR_AGE_NUMERIC" "INC_YEAR" "CAR_PRICE" "AGE" "TONNAGE" "HIGH_RISK_IND" "AIRBAG_JY" "DISP"][:]
    #warn("todo consider largest BRAND_NAME_JY")
    #excluded glass as per call with guanjun on October 23,2017.
    expl_cat=["PROD_CODE" "BRAND_NAME_JY" "ORG_CLASS" "CUST_TYPE" "CHANNEL_GROUP" "NCD_BFD" "IS_NEW" "MULTI_COV" "TRUCK_TYPE" "SERIES_LEVEL" "GRP_SEX" "GRP_RENEW"][:];
    forced_categorical_variables=["PROD_CODE"] #those that could be interpreted as a float64 by julia
    expl=uppercase.(vcat(expl_cat,expl_numeric))
    all_required_columns=vcat(expl,chosen_denom_col,chosen_num_col,chosen_weight_col,chosen_trn_val_col)
    #warn("BK: how to handle extremely sparse variables....")
    forty_largest_brands=["S66" "S376" "S87" "S24" "S293" "S382" "S95" "S41" "S20" "S70" "S396" "S34" "S04" "S10" "S284" "S259" "S163" "S13" "S96" "S173" "S33" "S23" "S280" "S398" "S183" "S113" "S177" "S238" "S22" "S261" "S18" "S339" "S132" "S445" "S105" "S263" "S441" "S223" "S252" "S446"][:]

@assert isfile(tmpf)
println("Reading data... \n $(tmpf)")
#@time df_tmp0=CSV.read(tmpf);
@time df_tmp=readtable(tmpf); #maybe around 15s for all the data
#enforce capitalized names
names!(df_tmp,Symbol.(uppercase.(string.(names(df_tmp)))));

#replace rare brands by "OTHER"
@time for i=1:size(df_tmp,1)
    di=df_tmp[i,:BRAND_NAME_JY]
    if di!="OTHER"
        if isna(di)
            @show "brand was na"
            df_tmp[i,:BRAND_NAME_JY]="missing"
        else
            if !in(di,forty_largest_brands)
                @show "brand was not in top 40"
                @show df_tmp[i,:BRAND_NAME_JY]
                @show in(di,forty_largest_brands)
                df_tmp[i,:BRAND_NAME_JY]="OTHER"
            end
        end
    end
end

#warn("removing cust type not equal 04")
keep_this_class_only="C01"
#keep_this_class_only="C04"
custtypeis01=df_tmp[Symbol("CUST_TYPE")].==keep_this_class_only
@show this_sum=sum(custtypeis01)
if this_sum>0
    @show this_sum
    @show this_sum./size(df_tmp,1)
    warn("BK: dropping rows which have cust_type different from $(keep_this_class_only) \t drop percentage: $(1-this_sum./size(df_tmp,1))")
    df_tmp=df_tmp[custtypeis01,:]
    @show new_size=size(df_tmp,1)
    @assert all(df_tmp[:CUST_TYPE].==keep_this_class_only)
else
    error("bk: this should not happen, no cust type 01 found.... why?")
end

#delete unused columns
    info("BK: deleting unused columns")
    for thiscol in names(df_tmp)
        if !(in(string(thiscol),all_required_columns))
            delete!(df_tmp,thiscol)
        end
    end

    warn("BK: dropping 'holdout' rows")
    ho_rows=df_tmp[Symbol(chosen_trn_val_col)].>1
    @show sm=sum(ho_rows)
    @show sm/size(df_tmp,1)
    df_tmp=df_tmp[.!ho_rows,:];

#rows with a zero denominator
    zero_denom_rows=df_tmp[Symbol(chosen_denom_col)].==0
    @show sm=sum(zero_denom_rows)
    if sm>0
        @show sm/size(df_tmp,1)
        warn("dropping rows with zero denominator")
        df_tmp=df_tmp[.!zero_denom_rows,:];
        info("we should compare the performance of this versus filter!")
    end

#prep df
@time (df_prepped,this_sett)=DecisionTrees.prepare_dataframe_for_dtm!(df_tmp,treat_as_categorical_variable=forced_categorical_variables,trnvalcol=chosen_trn_val_col,weightcol=chosen_weight_col,numcol=chosen_num_col,denomcol=chosen_denom_col,independent_vars=expl);

#rows with zero exposure/weight
    idx_zero_exposure=df_prepped[:weight].<0
    idx_missing_exposure=isna.(df_prepped[:weight])
    idx_missing_or_zero_exposure=idx_zero_exposure.|idx_missing_exposure
    @show count_zero_or_missing_exposure=sum(idx_missing_or_zero_exposure)
    @show count_zero_or_missing_exposure/size(df_prepped,1)
    if count_zero_or_missing_exposure>0
        warn("BK: dropping missing exposure rows!")
        df_prepped = df_prepped[.!idx_missing_or_zero_exposure,:];
    end

#remove all rows which have some NA entries
#an alternative is of course to replace the NAs with some specific value
tmp_names=names(df_prepped)
origrows=size(df_prepped,1)
for i=1:size(df_prepped,2)
    isnavec = isna.(df_prepped[i])
    na_count=sum(isnavec)
    if na_count>0
        df_prepped = df_prepped[.!isnavec,:];
        @show na_count,tmp_names[i]
    end
end
droppedrows=origrows-size(df_prepped,1)
if droppedrows>0
    warn("BK: some rows were dropped:")
    @show droppedrows
    @show size(df_prepped)
end

dtmtable=prep_data_from_df(df_prepped,this_sett,"C:\\temp\\tmp.txt");




datafolder="C:\\temp\\"# c:\\Users\\bernhard.konig\\Documents\\ASync\\publicnl\\Projects\\Active Projects\\China Guanjun Lead\\Working Files\\c04\\""
#peform CV for boosting
niter=50 #5 #125
mf=0.02
randomw=0.0
minw=-0.15
subsampling_features_prop=0.8
subsampling_prop=1.0
smoothEstimates=1
str_core=string("b_it",niter,"_mf",mf,"_rand",randomw,"_minw",minw,"_subsp",subsampling_features_prop,"_subsampl",subsampling_prop,"_smooth",smoothEstimates);
updateSettingsMod!(this_sett,minw=minw,write_dot_graph=true,write_csharp_code=true,write_vba_code=true,randomw=randomw,mf=mf,niter=niter,subsampling_features_prop=subsampling_features_prop,subsampling_prop=subsampling_prop,smoothEstimates=smoothEstimates,boolRandomizeOnlySplitAtTopNode="true",model_type="boosted_tree",scorebandsstartingpoints="1,100,200,300,400,500,600,700,800,900",showTimeUsedByEachIteration="false",moderationvector="$(mf)",write_result=false,print_details="false",boolSaveResultAsJLDFile="true",bool_write_tree="true",write_iteration_matrix="false",write_sas_code="true");
a,b=run_model_actual(dtmtable,this_sett,string(datafolder,"\\out.txt"));

