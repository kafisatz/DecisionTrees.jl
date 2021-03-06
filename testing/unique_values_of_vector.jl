using CSV
using Test 
using DecisionTrees 

elt = [Int,	Float64,	Float64,	Float64,	Float64,	Float64,	Int,	String,	String,	String,	Int,	String,	String,	String,	String,	String,	String,	String,	String,	String,	String,	Int,	String,	Int,	String,	String,	Int,	String,	String,	Int,	String,	String,	Int,	String,	String,	Int,	String,	Int,	String,	String,	Int,	Int,	String,	Int,	Int,	String,	Int,	Int,	String,	Int,	String,	String,	String,	Int,	String,	String,	String,	String,	Int,	Int,	Int,	String,	String,	String,	String,	String,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Float64,	Int,	Float64,	Float64,	Float64,	Float64,	Int,	Int,	Int,	Int,	Int,	Int]

fi = "data1Medium.csv"
fi = "data1Small.csv"
@time df_tmp = CSV.read(joinpath(datadir, "GermanMotorPremiums", fi), allowmissing=:none, types=elt, categorical=false, rows_for_type_detect=10000);

selected_explanatory_vars = ["PLZ_WOHNORT","ART_DES_WOHNEIGENTUM","GEBURTSDATUM","FAMILIENSTAND","NATIONALITAET","GESCHLECHT","FINANZIERUNGSART","STADT","KENNZEICHEN"]

# for x in selected_explanatory_vars
#    @show x,eltype(df_tmp[Symbol(x)])
# end

dtmtable, sett, df_prepped = prepare_dataframe_for_dtm!(df_tmp, treat_as_categorical_variable=["PLZ_WOHNORT"], weightcol="EXPOSURE", numcol="LOSS20HALF", denomcol="PREMIUM66", independent_vars=selected_explanatory_vars);
sett.minWeight = -.2


using StatsBase
using BenchmarkTools

v = df_prepped[:PLZ_WOHNORT]
# version 1
@benchmark uq1 = StatsBase.countmap(v)


# version 2
function me(y)
    u = unique(y)
    d = Dict([(i, count(x->x == i, y)) for i in u])
    return d
end
@benchmark uq2 = me(v) # this is about 40% faster than countmap. Then again


# sorting / cutoff
function frequencyRel(x)
    d = proportionmap(x)
    srt = sortperm(collect(values(d)))
    return d, srt
end

