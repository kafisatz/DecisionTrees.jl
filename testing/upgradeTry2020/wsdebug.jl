thisfile = joinpath(datadir, "GermanMotorPremiums", "data1Small.csv")
@test isfile(thisfile)
@time df_tmp = CSV.read(thisfile, allowmissing=:none, types=eltypesData1, categorical=false, rows_for_type_detect=10000);
df_tmp_orig = deepcopy(df_tmp)

# example with only char columns
selected_explanatory_vars_char_only = ["ART_DES_WOHNEIGENTUM","FAMILIENSTAND","GESCHLECHT","FINANZIERUNGSART","STADT","KENNZEICHEN"]
dtmtableMini, sett, df_prepped = prepare_dataframe_for_dtm!(df_tmp[1:100,:], weightcol="EXPOSURE", numcol="LOSS20HALF", denomcol="PREMIUM66", independent_vars=selected_explanatory_vars_char_only);
sett.minWeight = -.2
@info("DTM:Testing build tree (dtm mini; only character independent variables)")
strs, resm = dtm(dtmtableMini, sett)

features = dtmtableMini.features

featuresSubArr = view(features, collect(5:8), :)

i = 2

f = featuresSubArr[i]
featuresSubArr[i].parent.pool

fieldnames(typeof(f))

f.indices
f.indices

szTMP = size(f, 1)
scores = rand(szTMP);
numeratorEst = rand(szTMP);
numerator = rand(szTMP);
denominator = rand(szTMP);
weight = rand(szTMP);

