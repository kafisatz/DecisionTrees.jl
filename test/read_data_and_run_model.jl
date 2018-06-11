t0=time_ns()
cd("C:\\Users\\bernhard.konig\\Documents\\ASync\\home\\Code\\Julia\\DecisionTrees.jl")    

using Revise
using DataFrames,CSV
using DecisionTrees  

tela = (time_ns()-t0)/1e9
@show tela

@time df_tmp=readtable(string("data\\data1mini.csv")); 
#@time df_tmp=CSV.read(string("data\\data1mini.csv")); 

selected_explanatory_vars=["ART_DES_WOHNEIGENTUM","GEBURTSDATUM","FAMILIENSTAND","NATIONALITAET","PLZ_WOHNORT","GESCHLECHT","FINANZIERUNGSART","STADT","KENNZEICHEN"]

for x in selected_explanatory_vars
    @show x,eltype(df_tmp[Symbol(x)])
end

@time (df_prepped,sett)=prepare_dataframe_for_dtm!(df_tmp,treat_as_categorical_variable=["PLZ_WOHNORT"],weightcol="EXPOSURE",numcol="LOSS20",denomcol="PREMIUM66",independent_vars=selected_explanatory_vars);

res=prep_data_from_df(df_prepped,sett,"C:\\temp\\result.csv")