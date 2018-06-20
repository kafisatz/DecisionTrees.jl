#smoketest_tree.jl

@testset "Smoketests" begin

##################################################
#Read the data
##################################################
elt=[Int64,	Float64,	Float64,	Float64,	Float64,	Float64,	Int64,	String,	String,	String,	Int64,	String,	String,	String,	String,	String,	String,	String,	String,	String,	String,	Int64,	String,	Int64,	String,	String,	Int64,	String,	String,	Int64,	String,	String,	Int64,	String,	String,	Int64,	String,	Int64,	String,	String,	Int64,	Int64,	String,	Int64,	Int64,	String,	Int64,	Int64,	String,	Int64,	String,	String,	String,	Int64,	String,	String,	String,	String,	Int64,	Int64,	Int64,	String,	String,	String,	String,	String,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Float64,	Int64,	Float64,	Float64,	Float64,	Float64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64]

fi="data1Small.csv"
thisfile=joinpath(datadir,fi)
@test isfile(thisfile)
@time df_tmp=CSV.read(thisfile,allowmissing=:none,types=elt,categorical=false,rows_for_type_detect=10000);

selected_explanatory_vars=["PLZ_WOHNORT","ART_DES_WOHNEIGENTUM","GEBURTSDATUM","FAMILIENSTAND","NATIONALITAET","GESCHLECHT","FINANZIERUNGSART","STADT","KENNZEICHEN"]

#keep 10 largest PLZ only
vorig=deepcopy(df_tmp[:PLZ_WOHNORT])
counts,freqs,vals,keep=getCounts(vorig,threshold=10)
vnew=df_tmp[:PLZ_WOHNORT]
mapToOther!(vnew,keep,9999999)

##################################################
#Prepare the data
##################################################

dtmtable,sett,df_prepped=prepare_dataframe_for_dtm!(df_tmp,treat_as_categorical_variable=["PLZ_WOHNORT"],weightcol="EXPOSURE",numcol="LOSS20HALF",denomcol="PREMIUM66",independent_vars=selected_explanatory_vars);
sett.minw=-.2

#selected_explanatory_vars=[    "VORSCHAEDEN_ANZAHL",    "MALLORCA_POLICE",	"SCHUTZBRIEF_INKL",	"FREIE_WERKSTATTWAHL",	"AUTOMOBILCLUB_MITGLIED_SEIT",	"BAHNCARD",	"ZAHLUNGSWEISE",	"JAHRESKARTE_OEPNV",	"MOTORRAD_BESITZER",	"AUTOMOBILCLUB",	"SFKLASSE_VOLLKASKO",	"SFKLASSE_HAFTPFLICHT",	"STELLPLATZ_ABSCHLIESSBAR",	"NAECHTLICHER_STELLPLATZ",	"NUTZUNGSWEISE",	"JAEHRLICHE_FAHRLEISTUNG",	"TSN",	"ERSTZULASSUNG",	"HSN",	"FINANZIERUNGSART",	"ZULASSUNG_AUF_VERSICHERUNGSNEHM",	"STADT",	"KENNZEICHEN",	"PLZ_DES_HALTER",	"SELBSTGENUTZTES_WOHNEIGENTUM",	"ART_DES_WOHNEIGENTUM",	"GEBURTSDATUM",	"FAMILIENSTAND",	"NATIONALITAET",	"GESCHLECHT",	"FUEHRERSCHEIN_ERWORBEN_AM",	"VORSCHAEDEN0_typeKH",	"VORSCHAEDEN0_typetk",	"VORSCHAEDEN0_month",	"VORSCHAEDEN0_year",	"VORSCHAEDEN1_typetk",	"VORSCHAEDEN1_month",	"VORSCHAEDEN1_year",	"VORSCHAEDEN2_typevk",	"VORSCHAEDEN2_month",	"VORSCHAEDEN2_year",	"adacid",	"name",	"marke",	"modell",	"preis",	"getriebeart",	"antriebsart",	"Fahrzeugklasse",	"co2klasse",	"kw",	"ps",	"tueranzahl",	"Motorart",	"Kraftstoffart",	"Motorbauart",	"Schadstoffklasse",	"Karosserie",	"Sitzanzahl",	"typklasseh_num",	"typklassetk_num",	"typklassevk_num",	"hubraum2",	"drehmoment2",	"breite2",	"radstand2",	"laenge2",	"hoehe2",	"leergewicht2",	"gesamtgewicht2",	"zuladung2",	"kofferraumvolumen_num",	"hoechstgeschwindigkeit2",	"verbrauchgesamt2",	"verbrauchausserorts2",	"verbrauchinnerorts2",	"beschleunigung2",	"tank2",	"kfzsteuer2",	"anzahlgaenge2",	"anzahlzylinder2",	"co2_wert",	"modellstart_y"]
##################################################
#run tree
##################################################
strs,resm=dtm(dtmtable,sett)
@test typeof(resm.modelstats)==DataFrame
@test 1==1

##################################################
#run boosting
##################################################

sett.niter=10
sett.model_type="boosted_tree"
strs,resm2=dtm(dtmtable,sett)
@test typeof(resm2.modelstats)==DataFrame
@test 1==1


##################################################
#run bagging
##################################################
#tbd need to update bagging to work under j0.7 version of Decisiontrees
sett.model_type="bagged_tree"
#dtm(dtmtable,sett)
#@test 1==1


##################################################
#check case where we have more than 255 splitting points
##################################################

counts,freqs,vals,keep=getCounts(vorig,threshold=300)
vnew=deepcopy(vorig)
mapToOther!(vnew,keep,9999999)
df_tmp[:PLZ_WOHNORT]=vnew

dtmtable,sett,df_prepped=prepare_dataframe_for_dtm!(df_tmp,treat_as_categorical_variable=["PLZ_WOHNORT"],weightcol="EXPOSURE",numcol="LOSS20HALF",denomcol="PREMIUM66",independent_vars=selected_explanatory_vars);

sett.max_splitting_points_num=300
dtmtable=prep_data_from_df(df_prepped,sett,joinpath(mktempdir(),"dtmsmoketest.csv"))
@test eltype(dtmtable.features[:PLZ_WOHNORT].refs)==UInt16
strs,resm2b=dtm(dtmtable,sett)
@test typeof(resm2b.modelstats)==DataFrame

##################################################
#More Tests "for coverage" purposes
##################################################

##################################################
#Multirun
##################################################
settV=createGridSearchSettings(sett,    
    minw=[-0.1,-0.2]
    ,mf=[0.1,0.05],niter=[2])
a,b,allmodels=dtm(dtmtable,settV)
@test typeof(allmodels[1].modelstats)==DataFrame
    

##################################################
#Cross validation
##################################################
cvsampler=CVOptions(-3,0.0,true)
#this will run 5 models and thus might take some time
statsdf,settsdf,cvModels=dtm(dtmtable,sett,cvsampler)


#set some default settings
updateSettingsMod!(sett,
niter=4,
model_type="boosted_tree",
nscores="1000",
write_statistics="true",
write_sas_code="true",
write_iteration_matrix="true",
write_result="true",
write_csharp_code="true",
write_vba_code="true",
boolSaveJLDFile="true",
boolSaveResultAsJLDFile="false",

showProgressBar_time="true",
boolProduceEstAndLeafMatrices="true",
write_dot_graph="false",
boolCalculateGini="true",

boolCalculatePoissonError="true",
performanceMeasure="Average Poisson Error Val"
)

strs,resm3=dtm(dtmtable,sett)
@test typeof(resm3.modelstats)==DataFrame

#more smoketests
sett.subsampling_prop=0.7
sett.subsampling_features_prop=0.7
strs,resm4=dtm(dtmtable,sett)
@test typeof(resm4.modelstats)==DataFrame


sett.boolRandomizeOnlySplitAtTopNode=false
sett.randomw=0.02
strs,resm5=dtm(dtmtable,sett)
@test typeof(resm5.modelstats)==DataFrame

#minw too big
oldminw=sett.minw
sett.minw=-0.51
strs,resm6=dtm(dtmtable,sett)
@test typeof(resm6.modelstats)==DataFrame



sett.minw=oldminw
#example with weights.==1 and denominator.==1
dtmtable,sett,df_prepped=prepare_dataframe_for_dtm!(df_tmp,treat_as_categorical_variable=["PLZ_WOHNORT"],numcol="LOSS20HALF",independent_vars=selected_explanatory_vars);
    
   end #testset
