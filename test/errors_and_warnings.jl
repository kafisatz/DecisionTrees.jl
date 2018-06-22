#test errors
#prepare_dataframe_for_dtm

@testset "errorTesting" begin


fi="data1Small.csv"
thisfile=joinpath(datadir,fi)
@test isfile(thisfile)
@time df_tmp=CSV.read(thisfile,allowmissing=:none,types=eltypesData1,categorical=false,rows_for_type_detect=10000);

selected_explanatory_vars=["PLZ_WOHNORT","ART_DES_WOHNEIGENTUM","GEBURTSDATUM","FAMILIENSTAND","NATIONALITAET","GESCHLECHT","FINANZIERUNGSART","STADT","KENNZEICHEN"]

#keep 10 largest PLZ only
vorig=deepcopy(df_tmp[:PLZ_WOHNORT])
counts,freqs,vals,keep=getCounts(vorig,threshold=10)
vnew=df_tmp[:PLZ_WOHNORT]
mapToOther!(vnew,keep,9999999)

##################################################
#Prepare the data
##################################################

#examples of test_throws
@test_throws ErrorException error("some text")
@test_throws ErrorException("some text") error("some text") #this is more exact and also tests the text thrown
@test_throws AssertionError @assert false


dtmtable,sett,df_prepped=prepare_dataframe_for_dtm!(df_tmp,treat_as_categorical_variable=["PLZ_WOHNORT"],weightcol="EXPOSURE",numcol="LOSS20HALF",denomcol="PREMIUM66",independent_vars=selected_explanatory_vars);

@test_throws AssertionError prepare_dataframe_for_dtm!(df_tmp,treat_as_categorical_variable=["PLZ_WOHNORT"],weightcol="EXPOSURE",numcol="",denomcol="PREMIUM66",independent_vars=selected_explanatory_vars);
@test_throws AssertionError prepare_dataframe_for_dtm!(df_tmp,treat_as_categorical_variable=["PLZ_WOHNORT"],weightcol="EXPOSURE",numcol="does_not_exist",denomcol="PREMIUM66",independent_vars=selected_explanatory_vars);

@test_throws AssertionError prepare_dataframe_for_dtm!(df_tmp,treat_as_categorical_variable=["PLZ_WOHNORT"],weightcol="EXPOSURE",numcol="LOSS20HALF",denomcol="does_not_exist",independent_vars=selected_explanatory_vars);

@test_throws AssertionError prepare_dataframe_for_dtm!(df_tmp,treat_as_categorical_variable=["PLZ_WOHNORT"],weightcol="does_not_exist",numcol="LOSS20HALF",independent_vars=selected_explanatory_vars);

@test_throws AssertionError prepare_dataframe_for_dtm!(df_tmp,treat_as_categorical_variable=["PLZ_WOHNORT"],keycol="does_not_exist",numcol="LOSS20HALF",independent_vars=selected_explanatory_vars);

@test_throws ErrorException prepare_dataframe_for_dtm!(df_tmp,treat_as_categorical_variable=["does_not_exist","PLZ_WOHNORT"],numcol="LOSS20HALF",independent_vars=selected_explanatory_vars);


@test_throws Base.UVError prepare_dataframe_for_dtm!(df_tmp,directory="K:\\does_not_existIUDFLKJSDFEAA\\aaEFIEF_EFDXxT\\",treat_as_categorical_variable=["PLZ_WOHNORT"],numcol="LOSS20HALF",independent_vars=selected_explanatory_vars);

@test_throws AssertionError prepare_dataframe_for_dtm!(df_tmp,valpct=-22.1,treat_as_categorical_variable=["PLZ_WOHNORT"],numcol="LOSS20HALF",independent_vars=selected_explanatory_vars);

#@assert length(numcol)>0 "Error: No argument numcol provided."

end #testset