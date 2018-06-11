t0=time_ns()
cd(string(ENV["HOMEPATH"],"\\Documents\\ASync\\home\\Code\\Julia\\DecisionTrees.jl"))
@warn("You may need to run 'pkg> instantiate' when you first run this. Use ] to enther the package mode.")

using Revise
using JLD2
using DataFrames,CSV
using DecisionTrees  

tela = (time_ns()-t0)/1e9
@show tela #maybe 500-700 seconds on June 11, 2018

elt=[Int64,	Float64,	Float64,	Float64,	Int64,	String,	String,	String,	Int64,	String,	String,	String,	String,	String,	String,	String,	String,	String,	String,	Int64,	Int64,	Int64,	Int64,	String,	Int64,	String,	String,	Int64,	String,	String,	Int64,	String,	String,	Int64,	String,	Int64,	String,	String,	Int64,	Int64,	String,	Int64,	Int64,	String,	Int64,	Int64,	String,	Int64,	String,	String,	String,	Int64,	String,	String,	String,	String,	Int64,	Int64,	Int64,	String,	String,	String,	String,	String,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Float64,	Float64,	Float64,	Float64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64]

#@time df_tmp=readtable(string("data\\data1mini.csv"),eltypes=elt); 
@time df_tmp=CSV.read(string("data\\data1mini.csv"),allowmissing=:none,categorical=false); 

selected_explanatory_vars=["ART_DES_WOHNEIGENTUM","GEBURTSDATUM","FAMILIENSTAND","NATIONALITAET","PLZ_WOHNORT","GESCHLECHT","FINANZIERUNGSART","STADT","KENNZEICHEN"]

for x in selected_explanatory_vars
    @show x,eltype(df_tmp[Symbol(x)])
end

@time (df_prepped,sett)=prepare_dataframe_for_dtm!(df_tmp,treat_as_categorical_variable=["PLZ_WOHNORT"],weightcol="EXPOSURE",numcol="LOSS20",denomcol="PREMIUM66",independent_vars=selected_explanatory_vars);

#temporary disable this (as FileIO does not seem to work properly)
sett.boolSaveJLDFile=false

res=prep_data_from_df(df_prepped,sett,"C:\\temp\\")

rr=prepare_dataframe_for_dtm!(df_tmp,directory="C:\\temp\\mrx\\",treat_as_categorical_variable=["PLZ_WOHNORT"],weightcol="EXPOSURE",numcol="LOSS20",denomcol="PREMIUM66",independent_vars=selected_explanatory_vars);

#0