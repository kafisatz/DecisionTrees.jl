
using BenchmarkTools
using Revise    
using DecisionTrees
using DataFrames
using JLD2

pkgdir=Pkg.dir("DecisionTrees") #how long will this remain supported?
testdir=joinpath(pkgdir,"test")
datadir=joinpath(pkgdir,"data")

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


di9=load("C:\\Users\\bernhard.konig\\Documents\\ASync\\irobot\\data_tmp\\prepped_dtm2table.jld");dtmtable=di9["dtmtable"];this_dt2msett=di9["this_dt2msett"];

settstr="""model_type="boosted_tree",write_iteration_matrix="false",bool_write_tree="false",write_sas_code="false",write_statistics="false",write_iteration_matrix="false",minw=-.1,randomw=0.0,mf=.1,subsampling_features_prop=1,niter=2,subsampling_prop=1.00,graphvizexecutable="C:\\\\Program Files (x86)\\\\Graphviz2.38\\\\bin\\\\dot.exe",write_result=false"""
eval(Meta.parse(string("DTM2.updateSettingsMod!(this_dt2msett,",settstr,")")))
@time st2,dtm2model=DTM2.run_model_actual(dtmtable,this_dt2msett,fn);

settstr="""model_type="boosted_tree",boolProduceEstAndLeafMatrices=false,niter=50,bool_write_tree="false",write_sas_code="false",write_statistics="false",write_iteration_matrix="false",minw=-.08,randomw=0.0,mf=.1,subsampling_features_prop=1,subsampling_prop=1.00,graphvizexecutable="C:\\\\Program Files (x86)\\\\Graphviz2.38\\\\bin\\\\dot.exe",write_result=false"""
eval(Meta.parse(string("DTM2.updateSettingsMod!(this_dt2msett,",settstr,")")))

Profile.clear_malloc_data()
 st2,dtm2model=DTM2.run_model_actual(dtmtable,this_dt2msett,fn); #2.1s roughly
info("allocation run done.")
quit();



#=
run(`c:\\users\\bernhard.konig\\alloc.bat`)

tmpfilei="C:\\temp\\i.csv";isfile(tmpfilei)&&rm(tmpfilei)
mmfld="C:\\Users\\bernhard.konig\\Documents\\ASync\\irobot\\TariffWatch New Repo\\algorithms\\Julia\\Code\\DTM2_dev\\src\\"
using Coverage
mmm=analyze_malloc(mmfld);idx=map(x->x.bytes>0,mmm);mm2=mmm[idx];
sumKB=sum(map(x->x.bytes/1024, mm2));pct1=map(x->x.bytes/1024/sumKB,mm2);pct2=deepcopy(reverse(pct1));DTM2.incrementalToCumulative!(pct2);reverse!(pct2);
mfiles=map(x->x.filename,mm2);
mfiles_dironly=map(x->splitdir.(x)[1], mfiles);
mfiles_fileonly=map(x->splitdir.(x)[2], mfiles);
kbs=map(x->string(mm2[x].bytes/1024),eachindex(mm2));

#read source code
uq_files=unique(mfiles)
linesdict=Dict{String,Vector{String}}()
ff=uq_files[2]
for i=1:length(uq_files)
    ff=uq_files[i]
    lns=readlines(ff[1:end-4])
    linesdict[splitdir(ff)[2][1:end-4]]=deepcopy(lns)
end

ab=map(x->hcat(string(mm2[x].linenumber),linesdict[mfiles_fileonly[x][1:end-4]][mm2[x].linenumber],mfiles_fileonly[x],kbs[x],string(pct1[x]),string(pct2[x])), eachindex(mm2));writecsv("C:\\temp\\i.csv",ab);
@show "Done. Consider i.csv"
=#