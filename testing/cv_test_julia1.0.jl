using Revise 
st0 = time_ns()

import CSV
import DataFrames
import DataFrames: DataFrame
import Random
import OnlineStats
import StatsBase

using Distributed # only for Tests?
using Statistics
using Pkg # only for Tests?
using Test # only for Tests? 
using DecisionTrees 

# set directories
global pkgdir = ""
global testdir = joinpath(dirname(@__FILE__), "") # joinpath(@__DIR__,"..","test",filename)
# the following if clause is for "direct execution" (e.g. in VScode where @__FILE__ does not seem to work)
if testdir == ""
    try 
        global pkgdir = Pkg.dir("DecisionTrees") # how long will this remain supported?
        global testdir = joinpath(pkgdir, "test")
    catch
    end
else 
    global pkgdir = normpath(joinpath(testdir, ".."))
end
datadir = joinpath(pkgdir, "data")

Random.seed!(1239946)

t0runtests = time_ns()
println("DTM: Running tests...") # Time = $(t0runtests)")
    
global eltypesData1 = [Int,	Float64,	Float64,	Float64,	Float64,	Float64,	Int,	String,	String,	String,	Int,	String,	String,	String,	String,	String,	String,	String,	String,	String,	String,	Int,	String,	Int,	String,	String,	Int,	String,	String,	Int,	String,	String,	Int,	String,	String,	Int,	String,	Int,	String,	String,	Int,	Int,	String,	Int,	Int,	String,	Int,	Int,	String,	Int,	String,	String,	String,	Int,	String,	String,	String,	String,	Int,	Int,	Int,	String,	String,	String,	String,	String,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Float64,	Int,	Float64,	Float64,	Float64,	Float64,	Int,	Int,	Int,	Int,	Int,	Int]



import DataFrames
import Random: rand,randstring

##################################################
# Read the data
##################################################

thisfile = joinpath(datadir, "GermanMotorPremiums", "data1Small.csv")
@test isfile(thisfile)
@time df_tmp = CSV.read(thisfile, allowmissing=:none, types=eltypesData1, categorical=false, rows_for_type_detect=10000);
df_tmp_orig = deepcopy(df_tmp)

selected_explanatory_vars = ["PLZ_WOHNORT","ART_DES_WOHNEIGENTUM","GEBURTSDATUM","FAMILIENSTAND","NATIONALITAET","GESCHLECHT","FINANZIERUNGSART","STADT","KENNZEICHEN"]

# keep 10 largest PLZ only
vorig = deepcopy(df_tmp[:PLZ_WOHNORT])
counts, freqs, vals, keep = getCounts(vorig, threshold=10)
vnew = df_tmp[:PLZ_WOHNORT]
mapToOther!(vnew,keep,9999999)

##################################################
# Run Mini Tree
##################################################
dtmtableMini, sett, df_prepped = prepare_dataframe_for_dtm!(df_tmp[1:100,:], treat_as_categorical_variable=["PLZ_WOHNORT"], weightcol="EXPOSURE", numcol="LOSS20HALF", denomcol="PREMIUM66", independent_vars=selected_explanatory_vars);
sett.minWeight = -.2
# strs,resm=dtm(dtmtableMini,sett)
# boosting
settB = deepcopy(sett)
settB.model_type = "boosted_tree"
strs, resm = dtm(dtmtableMini, settB)

# selected_explanatory_vars=[    "VORSCHAEDEN_ANZAHL",    "MALLORCA_POLICE",	"SCHUTZBRIEF_INKL",	"FREIE_WERKSTATTWAHL",	"AUTOMOBILCLUB_MITGLIED_SEIT",	"BAHNCARD",	"ZAHLUNGSWEISE",	"JAHRESKARTE_OEPNV",	"MOTORRAD_BESITZER",	"AUTOMOBILCLUB",	"SFKLASSE_VOLLKASKO",	"SFKLASSE_HAFTPFLICHT",	"STELLPLATZ_ABSCHLIESSBAR",	"NAECHTLICHER_STELLPLATZ",	"NUTZUNGSWEISE",	"JAEHRLICHE_FAHRLEISTUNG",	"TSN",	"ERSTZULASSUNG",	"HSN",	"FINANZIERUNGSART",	"ZULASSUNG_AUF_VERSICHERUNGSNEHM",	"STADT",	"KENNZEICHEN",	"PLZ_DES_HALTER",	"SELBSTGENUTZTES_WOHNEIGENTUM",	"ART_DES_WOHNEIGENTUM",	"GEBURTSDATUM",	"FAMILIENSTAND",	"NATIONALITAET",	"GESCHLECHT",	"FUEHRERSCHEIN_ERWORBEN_AM",	"VORSCHAEDEN0_typeKH",	"VORSCHAEDEN0_typetk",	"VORSCHAEDEN0_month",	"VORSCHAEDEN0_year",	"VORSCHAEDEN1_typetk",	"VORSCHAEDEN1_month",	"VORSCHAEDEN1_year",	"VORSCHAEDEN2_typevk",	"VORSCHAEDEN2_month",	"VORSCHAEDEN2_year",	"adacid",	"name",	"marke",	"modell",	"preis",	"getriebeart",	"antriebsart",	"Fahrzeugklasse",	"co2klasse",	"kw",	"ps",	"tueranzahl",	"Motorart",	"Kraftstoffart",	"Motorbauart",	"Schadstoffklasse",	"Karosserie",	"Sitzanzahl",	"typklasseh_num",	"typklassetk_num",	"typklassevk_num",	"hubraum2",	"drehmoment2",	"breite2",	"radstand2",	"laenge2",	"hoehe2",	"leergewicht2",	"gesamtgewicht2",	"zuladung2",	"kofferraumvolumen_num",	"hoechstgeschwindigkeit2",	"verbrauchgesamt2",	"verbrauchausserorts2",	"verbrauchinnerorts2",	"beschleunigung2",	"tank2",	"kfzsteuer2",	"anzahlgaenge2",	"anzahlzylinder2",	"co2_wert",	"modellstart_y"]
##################################################
# run tree
##################################################
dtmtable, sett, df_prepped = prepare_dataframe_for_dtm!(df_tmp, treat_as_categorical_variable=["PLZ_WOHNORT"], weightcol="EXPOSURE", numcol="LOSS20HALF", denomcol="PREMIUM66", independent_vars=selected_explanatory_vars);
# strs,resm=dtm(dtmtable,sett)

##################################################
# run boosting
##################################################

sett.iterations = 10
sett.model_type = "boosted_tree"
# strs,resm2=dtm(dtmtable,sett)


# reduce iterations for the remaining tests
sett.iterations = 3

sett.statsByVariables = Int[1,2]
# strs,resm2=dtm(dtmtable,sett)
sett.statsByVariables = Int[]


sett.model_type = "boosted_tree"
cvsampler = CVOptions(-3, 0.0, true)
@info("DTM:testing cvsampler")


try
    statsdf, settsdf, cvModels = dtm(dtmtable, sett, cvsampler)
catch e 
    println("a")
    @show e
    @assert false
end 

cvsampler = CVOptions(3, 0.5, false)
try
    statsdf, settsdf, cvModels = dtm(dtmtableMini, sett, cvsampler)
catch e 
    println("b")
    @show e
    @assert false
end 