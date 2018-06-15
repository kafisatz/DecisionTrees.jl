t0=time_ns()
cd(string(ENV["HOMEPATH"],"\\Documents\\ASync\\home\\Code\\Julia\\DecisionTrees.jl"))
@warn("You may need to run 'pkg> instantiate' when you first run this. Use ] to enter the package mode.")

using Revise

#these lines are primarly to test the precompliation of the dependency packages
using FileIO
using DataFrames
using PyCall 
using CSV 
using Compat
using OnlineStats 
using StatsBase
#using Distributions
using ProgressMeter
using JLD2
using FileIO
using StatsFuns 

using DecisionTrees  

tela = (time_ns()-t0)/1e9
@show tela #maybe 500-700 seconds on June 11, 2018

elt=[Int64,	Float64,	Float64,	Float64,	Int64,	String,	String,	String,	Int64,	String,	String,	String,	String,	String,	String,	String,	String,	String,	String,	Int64,	Int64,	Int64,	Int64,	String,	Int64,	String,	String,	Int64,	String,	String,	Int64,	String,	String,	Int64,	String,	Int64,	String,	String,	Int64,	Int64,	String,	Int64,	Int64,	String,	Int64,	Int64,	String,	Int64,	String,	String,	String,	Int64,	String,	String,	String,	String,	Int64,	Int64,	Int64,	String,	String,	String,	String,	String,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Float64,	Int64,	Float64,	Float64,	Float64,	Float64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64]

#@time df_tmp=readtable(string("data\\data1mini.csv"),eltypes=elt); 
@time df_tmp=CSV.read(string("data\\data1mini.csv"),allowmissing=:none,categorical=false); 

selected_explanatory_vars=["ART_DES_WOHNEIGENTUM","GEBURTSDATUM","FAMILIENSTAND","NATIONALITAET","PLZ_WOHNORT","GESCHLECHT","FINANZIERUNGSART","STADT","KENNZEICHEN"]

for x in selected_explanatory_vars
    @show x,eltype(df_tmp[Symbol(x)])
end

#@time (df_prepped,sett)=prepare_dataframe_for_dtm!(df_tmp,treat_as_categorical_variable=["PLZ_WOHNORT"],weightcol="EXPOSURE",numcol="LOSS20",denomcol="PREMIUM66",independent_vars=selected_explanatory_vars);

#temporary disable this (as FileIO does not seem to work properly)
#sett.boolSaveJLDFile=false

#res=prep_data_from_df(df_prepped,sett,"C:\\temp\\")

dtmtable,sett=prepare_dataframe_for_dtm!(df_tmp,directory="C:\\temp\\",treat_as_categorical_variable=["PLZ_WOHNORT"],weightcol="EXPOSURE",numcol="LOSS20",denomcol="PREMIUM66",independent_vars=selected_explanatory_vars);
sett.minw=-.2
#0

@time df_tmp=CSV.read(string("data\\data1medium.csv"),allowmissing=:none,categorical=false,rows_for_type_detect=10000);

selected_explanatory_vars=[    "VORSCHAEDEN_ANZAHL",    "MALLORCA_POLICE",	"SCHUTZBRIEF_INKL",	"FREIE_WERKSTATTWAHL",	"AUTOMOBILCLUB_MITGLIED_SEIT",	"BAHNCARD",	"ZAHLUNGSWEISE",	"JAHRESKARTE_OEPNV",	"MOTORRAD_BESITZER",	"AUTOMOBILCLUB",	"SFKLASSE_VOLLKASKO",	"SFKLASSE_HAFTPFLICHT",	"STELLPLATZ_ABSCHLIESSBAR",	"NAECHTLICHER_STELLPLATZ",	"NUTZUNGSWEISE",	"JAEHRLICHE_FAHRLEISTUNG",	"TSN",	"ERSTZULASSUNG",	"HSN",	"FINANZIERUNGSART",	"ZULASSUNG_AUF_VERSICHERUNGSNEHM",	"STADT",	"KENNZEICHEN",	"PLZ_DES_HALTER",	"SELBSTGENUTZTES_WOHNEIGENTUM",	"ART_DES_WOHNEIGENTUM",	"GEBURTSDATUM",	"FAMILIENSTAND",	"NATIONALITAET",	"GESCHLECHT",	"FUEHRERSCHEIN_ERWORBEN_AM",	"VORSCHAEDEN0_typeKH",	"VORSCHAEDEN0_typetk",	"VORSCHAEDEN0_month",	"VORSCHAEDEN0_year",	"VORSCHAEDEN1_typetk",	"VORSCHAEDEN1_month",	"VORSCHAEDEN1_year",	"VORSCHAEDEN2_typevk",	"VORSCHAEDEN2_month",	"VORSCHAEDEN2_year",	"adacid",	"name",	"marke",	"modell",	"preis",	"getriebeart",	"antriebsart",	"Fahrzeugklasse",	"co2klasse",	"kw",	"ps",	"tueranzahl",	"Motorart",	"Kraftstoffart",	"Motorbauart",	"Schadstoffklasse",	"Karosserie",	"Sitzanzahl",	"typklasseh_num",	"typklassetk_num",	"typklassevk_num",	"hubraum2",	"drehmoment2",	"breite2",	"radstand2",	"laenge2",	"hoehe2",	"leergewicht2",	"gesamtgewicht2",	"zuladung2",	"kofferraumvolumen_num",	"hoechstgeschwindigkeit2",	"verbrauchgesamt2",	"verbrauchausserorts2",	"verbrauchinnerorts2",	"beschleunigung2",	"tank2",	"kfzsteuer2",	"anzahlgaenge2",	"anzahlzylinder2",	"co2_wert",	"modellstart_y"]

for x in selected_explanatory_vars
    @show x,eltype(df_tmp[Symbol(x)])
end

#@time (df_prepped,sett)=prepare_dataframe_for_dtm!(df_tmp,treat_as_categorical_variable=["PLZ_WOHNORT"],weightcol="EXPOSURE",numcol="LOSS20",denomcol="PREMIUM66",independent_vars=selected_explanatory_vars);

#temporary disable this (as FileIO does not seem to work properly)
#sett.boolSaveJLDFile=false

#res=prep_data_from_df(df_prepped,sett,"C:\\temp\\")
features=dtmtable.features
poolsizes=map(x->length(features[x].pool),1:size(features,2))
pldidx=poolsizes.>255
sett.df_name_vector[pldidx]

dtmtable,sett=prepare_dataframe_for_dtm!(df_tmp,directory="C:\\temp\\",treat_as_categorical_variable=["PLZ_WOHNORT"],weightcol="EXPOSURE",numcol="LOSS20HALF",denomcol="PREMIUM66",independent_vars=selected_explanatory_vars);
dtm(dtmtable,sett)
sett.minw=-.2