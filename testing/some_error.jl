dtmtableMini, sett, df_prepped = prepare_dataframe_for_dtm!(df_tmp[1:100,:], treat_as_categorical_variable=["PLZ_WOHNORT"], weightcol="EXPOSURE", numcol="LOSS20HALF", denomcol="PREMIUM66", independent_vars=selected_explanatory_vars);



function prepare_dataframe_for_dtm!(dfin::DataFrame;directory::String=mktempdir(),treat_as_categorical_variable::Vector{String}=Vector{String}(),numcol::String="",denomcol::String="",weightcol::String="",trnvalcol::String="",valpct::Float64=0.3,keycol::String="",independent_vars::Vector{String}=Vector{String}(),methodForSplittingPointsSelection::String="basedOnWeightVector")


    dfin = df_tmp[1:100,:]
    treat_as_categorical_variable = ["PLZ_WOHNORT"]
    weightcol = "EXPOSURE"
    numcol = "LOSS20HALF"
    denomcol = "PREMIUM66"
    keycol = "" 
    trnvalcol = ""
    methodForSplittingPointsSelection = "basedOnWeightVector"
    independent_vars = selected_explanatory_vars
    valpct = .3




    function prepareDF!(dfin::DataFrame;treat_as_categorical_variable::Vector{String}=Vector{String}(),numcol::String="",denomcol::String="",weightcol::String="",trnvalcol::String="",valpct::Float64=0.3,keycol::String="",independent_vars::Vector{String}=Vector{String}())


        DecisionTrees.prepareDF!(dfin, treat_as_categorical_variable=treat_as_categorical_variable, numcol=numcol, denomcol=denomcol, weightcol=weightcol, trnvalcol=trnvalcol, valpct=valpct, keycol=keycol, independent_vars=independent_vars);                    