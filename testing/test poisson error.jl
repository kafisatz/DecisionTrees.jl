fi="data1medium.csv"
fi="data1small.csv"
@time df_tmp=CSV.read(joinpath(datadir,fi),allowmissing=:none,types=elt,categorical=false,rows_for_type_detect=10000);

selected_explanatory_vars=["ART_DES_WOHNEIGENTUM","GEBURTSDATUM","FAMILIENSTAND","NATIONALITAET","GESCHLECHT","FINANZIERUNGSART","STADT","KENNZEICHEN"]

#for x in selected_explanatory_vars
#    @show x,eltype(df_tmp[Symbol(x)])
#end

df_tmp_orig=deepcopy(df_tmp)
df_tmp=df_tmp_orig[1:20,:]

dtmtable,sett=prepare_dataframe_for_dtm!(df_tmp,weightcol="EXPOSURE",numcol="LOSS20HALF",denomcol="PREMIUM66",independent_vars=selected_explanatory_vars);
sett.minw=-.2
#0

#selected_explanatory_vars=[    "VORSCHAEDEN_ANZAHL",    "MALLORCA_POLICE",	"SCHUTZBRIEF_INKL",	"FREIE_WERKSTATTWAHL",	"AUTOMOBILCLUB_MITGLIED_SEIT",	"BAHNCARD",	"ZAHLUNGSWEISE",	"JAHRESKARTE_OEPNV",	"MOTORRAD_BESITZER",	"AUTOMOBILCLUB",	"SFKLASSE_VOLLKASKO",	"SFKLASSE_HAFTPFLICHT",	"STELLPLATZ_ABSCHLIESSBAR",	"NAECHTLICHER_STELLPLATZ",	"NUTZUNGSWEISE",	"JAEHRLICHE_FAHRLEISTUNG",	"TSN",	"ERSTZULASSUNG",	"HSN",	"FINANZIERUNGSART",	"ZULASSUNG_AUF_VERSICHERUNGSNEHM",	"STADT",	"KENNZEICHEN",	"PLZ_DES_HALTER",	"SELBSTGENUTZTES_WOHNEIGENTUM",	"ART_DES_WOHNEIGENTUM",	"GEBURTSDATUM",	"FAMILIENSTAND",	"NATIONALITAET",	"GESCHLECHT",	"FUEHRERSCHEIN_ERWORBEN_AM",	"VORSCHAEDEN0_typeKH",	"VORSCHAEDEN0_typetk",	"VORSCHAEDEN0_month",	"VORSCHAEDEN0_year",	"VORSCHAEDEN1_typetk",	"VORSCHAEDEN1_month",	"VORSCHAEDEN1_year",	"VORSCHAEDEN2_typevk",	"VORSCHAEDEN2_month",	"VORSCHAEDEN2_year",	"adacid",	"name",	"marke",	"modell",	"preis",	"getriebeart",	"antriebsart",	"Fahrzeugklasse",	"co2klasse",	"kw",	"ps",	"tueranzahl",	"Motorart",	"Kraftstoffart",	"Motorbauart",	"Schadstoffklasse",	"Karosserie",	"Sitzanzahl",	"typklasseh_num",	"typklassetk_num",	"typklassevk_num",	"hubraum2",	"drehmoment2",	"breite2",	"radstand2",	"laenge2",	"hoehe2",	"leergewicht2",	"gesamtgewicht2",	"zuladung2",	"kofferraumvolumen_num",	"hoechstgeschwindigkeit2",	"verbrauchgesamt2",	"verbrauchausserorts2",	"verbrauchinnerorts2",	"beschleunigung2",	"tank2",	"kfzsteuer2",	"anzahlgaenge2",	"anzahlzylinder2",	"co2_wert",	"modellstart_y"]
#run tree
#strs,resM=dtm(dtmtable,sett)

sett.model_type="boosted_tree"
sett.niter=4
sett.boolCalculatePoissonError=true
strs,resM=dtm(dtmtable,sett)

fitted=resM.meanobserved.*resM.rawrelativities
mean(fitted[dtmtable.trnidx]) #0.8405367103402663
mean(dtmtable.numerator[dtmtable.trnidx]) 
mean(dtmtable.weight[dtmtable.trnidx])
errs=poissonError(dtmtable.numerator,dtmtable.weight,fitted);
trnAvgError=sum(errs[dtmtable.trnidx])/length(dtmtable.trnidx) #10055.619876087196
valAvgError=sum(errs[dtmtable.validx])/length(dtmtable.validx) 