# smoketest_tree.jl

@testset "Smoketests" begin

    import DataFrames
    import Random: rand,randstring

##################################################
# Read the data
##################################################

    thisfile = joinpath(datadir, "GermanMotorPremiums", "data1Small.csv")
    @test isfile(thisfile)
    @time df_tmp = CSV.read(thisfile, DataFrame,strict=true, types=eltypesData1,  pool=true, copycols=true);
    df_tmp_orig = deepcopy(df_tmp)

# example with only char columns
    selected_explanatory_vars_char_only = ["ART_DES_WOHNEIGENTUM","FAMILIENSTAND","GESCHLECHT","FINANZIERUNGSART","STADT","KENNZEICHEN"]
    dtmtableMini, sett, df_prepped = prepare_dataframe_for_dtm!(df_tmp[1:100,:], weightcol="EXPOSURE", numcol="LOSS20HALF", denomcol="PREMIUM66", independent_vars=selected_explanatory_vars_char_only);
    sett.minWeight = -.2
    @info("DTM:Testing build tree (dtm mini; only character independent variables)")
    strs, resm = dtm(dtmtableMini, sett)
    @test true

    settB = deepcopy(sett)
    settB.model_type = "boosted_tree"
    @info("DTM:Testing boosted tree (dtm mini; only character independent variables)")
    strs, resm = dtm(dtmtableMini, settB)
    @test true

# example with only numeric columns
    selected_explanatory_vars_num_only = ["GEBURTSDATUM","VORSCHAEDEN_ANZAHL","AUTOMOBILCLUB_MITGLIED_SEIT"]
    dtmtableMini, sett, df_prepped = prepare_dataframe_for_dtm!(df_tmp[1:100,:], weightcol="EXPOSURE", numcol="LOSS20HALF", denomcol="PREMIUM66", independent_vars=selected_explanatory_vars_num_only);
    sett.minWeight = -.2
    @info("DTM:Testing build tree (dtm mini; only numerical independent variables)")
    strs, resm = dtm(dtmtableMini, sett)
    @test true

    settB = deepcopy(sett)
    settB.model_type = "boosted_tree"
    @info("DTM:Testing boosted tree (dtm mini; only numerical independent variables)")
    strs, resm = dtm(dtmtableMini, settB)


# other tests:
    selected_explanatory_vars = ["PLZ_WOHNORT","ART_DES_WOHNEIGENTUM","GEBURTSDATUM","FAMILIENSTAND","NATIONALITAET","GESCHLECHT","FINANZIERUNGSART","STADT","KENNZEICHEN"]

# keep 10 largest PLZ only
    vorig = deepcopy(df_tmp[!,:PLZ_WOHNORT])
    counts, freqs, vals, keep = getCounts(vorig, threshold=10)
    vnew = df_tmp[!,:PLZ_WOHNORT]
    mapToOther!(vnew, keep, 9999999)

##################################################
# Run Mini Tree
##################################################
    dtmtableMini, sett, df_prepped = prepare_dataframe_for_dtm!(df_tmp[1:100,:], treat_as_categorical_variable=["PLZ_WOHNORT"], weightcol="EXPOSURE", numcol="LOSS20HALF", denomcol="PREMIUM66", independent_vars=selected_explanatory_vars);
    sett.minWeight = -.2
    @info("DTM:Testing build tree (dtm mini)")
    strs, resm = dtm(dtmtableMini, sett)

# boosting
    settB = deepcopy(sett)
    settB.model_type = "boosted_tree"
    @info("DTM:Testing boosted tree (dtm mini)")
    strs, resm = dtm(dtmtableMini, settB)

############################################
# graphviz 
############################################

    local graphvizexe 
    if Sys.iswindows()
        graphvizexe = "dot.exe"
    else 
        graphvizexe = "dot"
    end
    if !isfile(graphvizexe)
        if haskey(ENV,"graphvizdot")
            graphvizexe = ENV["graphvizdot"]
        end 
    end 
    if isfile(graphvizexe)
        #do graphviz tests 

            dtmtableMini, sett, df_prepped = prepare_dataframe_for_dtm!(df_tmp[1:100,:], treat_as_categorical_variable=["PLZ_WOHNORT"], weightcol="EXPOSURE", numcol="LOSS20HALF", denomcol="PREMIUM66", independent_vars=selected_explanatory_vars);
            sett.write_dot_graph = true
            sett.graphvizexecutable = graphvizexe
            sett.minWeight = -.2
            @info("DTM:Testing build tree (dtm mini)")
            strs, resm = dtm(dtmtableMini, sett)
            @test true
            pdf_files_produced = filter(x->endswith(lowercase(x),".pdf"),strs)
            @test size(pdf_files_produced,1)==1
            
        # boosting
            settB = deepcopy(sett)
            settB.model_type = "boosted_tree"
            @info("DTM:Testing boosted tree (dtm mini)")
            strs, resm = dtm(dtmtableMini, settB)
            @test true
            pdf_files_produced = filter(x->endswith(lowercase(x),".pdf"),strs)
            @test size(pdf_files_produced,1)==2
    
    else 
        #GraphViz Exe not found
        @test false 
        @warn("""Graphvizexecutable not found. Please specify the path to dot.exe in ENV["graphvizdot"]""")            
    end


# selected_explanatory_vars=[    "VORSCHAEDEN_ANZAHL",    "MALLORCA_POLICE",	"SCHUTZBRIEF_INKL",	"FREIE_WERKSTATTWAHL",	"AUTOMOBILCLUB_MITGLIED_SEIT",	"BAHNCARD",	"ZAHLUNGSWEISE",	"JAHRESKARTE_OEPNV",	"MOTORRAD_BESITZER",	"AUTOMOBILCLUB",	"SFKLASSE_VOLLKASKO",	"SFKLASSE_HAFTPFLICHT",	"STELLPLATZ_ABSCHLIESSBAR",	"NAECHTLICHER_STELLPLATZ",	"NUTZUNGSWEISE",	"JAEHRLICHE_FAHRLEISTUNG",	"TSN",	"ERSTZULASSUNG",	"HSN",	"FINANZIERUNGSART",	"ZULASSUNG_AUF_VERSICHERUNGSNEHM",	"STADT",	"KENNZEICHEN",	"PLZ_DES_HALTER",	"SELBSTGENUTZTES_WOHNEIGENTUM",	"ART_DES_WOHNEIGENTUM",	"GEBURTSDATUM",	"FAMILIENSTAND",	"NATIONALITAET",	"GESCHLECHT",	"FUEHRERSCHEIN_ERWORBEN_AM",	"VORSCHAEDEN0_typeKH",	"VORSCHAEDEN0_typetk",	"VORSCHAEDEN0_month",	"VORSCHAEDEN0_year",	"VORSCHAEDEN1_typetk",	"VORSCHAEDEN1_month",	"VORSCHAEDEN1_year",	"VORSCHAEDEN2_typevk",	"VORSCHAEDEN2_month",	"VORSCHAEDEN2_year",	"adacid",	"name",	"marke",	"modell",	"preis",	"getriebeart",	"antriebsart",	"Fahrzeugklasse",	"co2klasse",	"kw",	"ps",	"tueranzahl",	"Motorart",	"Kraftstoffart",	"Motorbauart",	"Schadstoffklasse",	"Karosserie",	"Sitzanzahl",	"typklasseh_num",	"typklassetk_num",	"typklassevk_num",	"hubraum2",	"drehmoment2",	"breite2",	"radstand2",	"laenge2",	"hoehe2",	"leergewicht2",	"gesamtgewicht2",	"zuladung2",	"kofferraumvolumen_num",	"hoechstgeschwindigkeit2",	"verbrauchgesamt2",	"verbrauchausserorts2",	"verbrauchinnerorts2",	"beschleunigung2",	"tank2",	"kfzsteuer2",	"anzahlgaenge2",	"anzahlzylinder2",	"co2_wert",	"modellstart_y"]
##################################################
# run tree
##################################################
    dtmtable, sett, df_prepped = prepare_dataframe_for_dtm!(df_tmp, treat_as_categorical_variable=["PLZ_WOHNORT"], weightcol="EXPOSURE", numcol="LOSS20HALF", denomcol="PREMIUM66", independent_vars=selected_explanatory_vars);
    @info("DTM:Testing build tree")
    strs, resm = dtm(dtmtable, sett)
    @test typeof(resm.modelstats) == DataFrame
    @test 1 == 1

# predict 
    DecisionTrees.predict(resm, dtmtable.features)  

##################################################
# run boosting
##################################################

    sett.iterations = 10
    sett.model_type = "boosted_tree"
    @info("DTM:Testing boosted tree. Iterations=$(sett.iterations)")
    strs, resm2 = dtm(dtmtable, sett)
    @test typeof(resm2.modelstats) == DataFrame
    @test 1 == 1

# reduce iterations for the remaining tests
    sett.iterations = 3

    sett.statsByVariables = Int[1,2]
    @info("DTM:Testing statsByVariables")
    strs, resm2 = dtm(dtmtable, sett)
    sett.statsByVariables = Int[]

    @info("DTM:Testing different splitting criterias...")

# predict 
DecisionTrees.predict(resm2, dtmtable.features)  

# try different splitting criteria
    for splitCrit in ["difference","poissondeviance","gammadeviance","mse","sse","msepointwise"], modelTYPE in ["build_tree","boosted_tree"]
        updateSettingsMod!(sett, crit=splitCrit, model_type=modelTYPE)
        try 
            strs, resmT = dtm(dtmtable, sett)
            @test typeof(resmT.modelstats) == DataFrame
        catch thisErr
            @show stacktrace()
            @show thisErr        
            @test "This one failed->" == string("crit=", splitCrit, "; model_type=", modelTYPE)
        end
    end

    # try different splitting criteria with randomw > 0 
    for splitCrit in ["difference","poissondeviance","gammadeviance"], modelTYPE in ["build_tree","boosted_tree"]
        updateSettingsMod!(sett, crit=splitCrit, model_type=modelTYPE)
        updateSettingsMod!(sett, randomw = 0.05)
        try 
            strs, resmT = dtm(dtmtable, sett)        
            @test typeof(resmT.modelstats) == DataFrame
        catch thisErr
            @show stacktrace()
            @show thisErr        
            @test "This one failed->" == string("crit=", splitCrit, "; model_type=", modelTYPE, "; randomw=", rw)
        end
    end

    updateSettingsMod!(sett, crit="difference")
    @info("DTM:Finished testing different splitting criterias.")

##################################################
# run bagging
##################################################
# tbd need to update bagging to work under j0.7 version of Decisiontrees
    sett.model_type = "bagged_tree"
# @info("DTM:testing bagging")
# dtm(dtmtable,sett)
# @test 1==1


##################################################
# check case where we have more than 255 splitting points
##################################################

    counts, freqs, vals, keep = getCounts(vorig, threshold=300)
    vnew = deepcopy(vorig)
    mapToOther!(vnew, keep, 9999999)
    df_tmp[!,:PLZ_WOHNORT] = vnew

    dtmtable, sett, df_prepped = prepare_dataframe_for_dtm!(df_tmp, treat_as_categorical_variable=["PLZ_WOHNORT"], weightcol="EXPOSURE", numcol="LOSS20HALF", denomcol="PREMIUM66", independent_vars=selected_explanatory_vars);

    sett.maxSplittingPoints = 300
    dtmtable = prep_data_from_df(df_prepped, sett, joinpath(mktempdir(), "dtmsmoketest.csv"))
    @test eltype(dtmtable.features[!,:PLZ_WOHNORT].refs) == UInt16
    @info("DTM:testing maxSplittingPoints=$(sett.maxSplittingPoints)")
    strs, resm2b = dtm(dtmtable, sett)
    @test typeof(resm2b.modelstats) == DataFrame

##################################################
# More Tests "for coverage" purposes
##################################################

##################################################
# Multirun Boosting
##################################################
    settV = createGridSearchSettings(sett,    
    minWeight=[-0.1,-0.2]
    ,learningRate=[0.1,0.05],iterations=[2])
    settVBoosting = deepcopy(settV)
    @info("DTM:testing gridsearch for boosting")
    a, b, allmodels = dtm(dtmtable, settV)
    @test typeof(allmodels[1].modelstats) == DataFrame
    

##################################################
# Multirun buildtree
##################################################
    sett.model_type = "build_tree"
    settV = createGridSearchSettings(sett,    
    minWeight=[-0.1,-0.2]
    ,learningRate=[0.1,0.05],iterations=[2])
    @info("DTM:testing gridsearch/multirun for build_tree")
    a, b, allmodels = dtm(dtmtable, settV)
# @test typeof(allmodels[1].modelstats)==DataFrame

##################################################
# Cross validation
##################################################
    sett.model_type = "boosted_tree"
    cvsampler = CVOptions(-3, 0.0, true)
    @info("DTM:testing cvsampler")
    statsdf, settsdf, cvModels = dtm(dtmtable, sett, cvsampler)

    cvsampler = CVOptions(3, 0.65, true)
    @info("DTM:testing cvsampler2")
    statsdf, settsdf, cvModels = dtm(dtmtable, sett, cvsampler)

    @warn("We can possibly enable the following snippet now (it should run fine)")
    @warn("CV sampler tests were breaking the other tests (possibly because dtmtable is modified which should not happen!")
    if true # false
        @info("DTM:testing cvsampler3")
        cvsampler = CVOptions(3, 0.5, false)
        statsdf, settsdf, cvModels = dtm(dtmtable, sett, cvsampler)
    else 
        @test_broken 1 == 2 
    end

    @info("DTM:testing different performance measures")
    cvsampler = CVOptions(-3, 0.0, true)
    for thisPerfMeasure in DecisionTrees.global_statsperiter_header
        updateSettingsMod!(sett, performanceMeasure=thisPerfMeasure)
        statsdf = DataFrame()
        statsdf, settsdf, cvModels = dtm(dtmtableMini, sett, cvsampler)
        @test size(statsdf, 1) > 0
    end

# set some default settings
    updateSettingsMod!(sett,
iterations=3,
model_type="boosted_tree",
nScores="1000",
writeStatistics="true",
writeSasCode="true",
writeIterationMatrix="true",
writeResult="true",
writeCsharpCode="true",
writeVbaCode="true",
saveJLDFile="true",
saveResultAsJLDFile="false",

showProgressBar_time="true",
prroduceEstAndLeafMatrices="true",
write_dot_graph="false",
calculateGini="true",

calculatePoissonError="true",
performanceMeasure="Average Poisson Error Val"
)
    sett.model_type = "build_tree"
    @info("DTM:testing a build_tree model")
    strs, resm3 = dtm(dtmtable, sett)
    sett.model_type = "boosted_tree"
    @info("DTM:testing a boosted_tree")
    strs, resm3 = dtm(dtmtable, sett)
    @test typeof(resm3.modelstats) == DataFrame
# test if files exist 
    for fi in strs
        @test isfile(fi)
    end

# more smoketests
    sett.subsampling_prop = 0.7
    sett.subsampling_features_prop = 0.7
    @info("DTM:testing subsampling_features_prop and subsampling_prop")
    strs, resm4 = dtm(dtmtable, sett)
    @test typeof(resm4.modelstats) == DataFrame


    sett.boolRandomizeOnlySplitAtTopNode = false
    sett.randomw = 0.02
    @info("DTM:testing randomw=$(sett.randomw)")
    strs, resm5 = dtm(dtmtable, sett)
    @test typeof(resm5.modelstats) == DataFrame

    sett.spawnsmaller = false
    sett.randomw = 0.0
    @info("DTM:testing spawnsmaller=$(sett.spawnsmaller) randomw=$(sett.randomw)")
    strs, resm5b = dtm(dtmtable, sett)
    @test typeof(resm5b.modelstats) == DataFrame
    
# minWeight too big
    oldminw = sett.minWeight
    sett.minWeight = -0.51
    @info("DTM:testing minWeight=$(sett.minWeight)")
    strs, resm6 = dtm(dtmtable, sett)
    @test typeof(resm6.modelstats) == DataFrame

# for coverage 
# empty keycol
    @info("DTM:testing errors/warnings of prepare_dataframe_for_dtm")
    prepare_dataframe_for_dtm!(df_tmp, treat_as_categorical_variable=["PLZ_WOHNORT"], keycol="", numcol="LOSS20HALF", independent_vars=selected_explanatory_vars);
    df_tmp2 = deepcopy(df_tmp_orig)
    df_tmp2[!,:PREMIUM66][1:20] .= -20.1
    dtmtable, sett, df_prepped = prepare_dataframe_for_dtm!(df_tmp, treat_as_categorical_variable=["PLZ_WOHNORT"], weightcol="EXPOSURE", numcol="LOSS20HALF", denomcol="PREMIUM66", independent_vars=selected_explanatory_vars);
    df_tmp2[!,:PREMIUM66][1:20] .= 0.0
    dtmtable, sett, df_prepped = prepare_dataframe_for_dtm!(df_tmp, treat_as_categorical_variable=["PLZ_WOHNORT"], weightcol="EXPOSURE", numcol="LOSS20HALF", denomcol="PREMIUM66", independent_vars=selected_explanatory_vars);
    df_tmp2[!,:LOSS20HALF][1:20] .= 0.0
    dtmtable, sett, df_prepped = prepare_dataframe_for_dtm!(df_tmp, treat_as_categorical_variable=["PLZ_WOHNORT"], weightcol="EXPOSURE", numcol="LOSS20HALF", denomcol="PREMIUM66", independent_vars=selected_explanatory_vars);
    df_tmp2[!,:LOSS20HALF][1:20] .= -20
    dtmtable, sett, df_prepped = prepare_dataframe_for_dtm!(df_tmp, treat_as_categorical_variable=["PLZ_WOHNORT"], weightcol="EXPOSURE", numcol="LOSS20HALF", denomcol="PREMIUM66", independent_vars=selected_explanatory_vars);


    sett.minWeight = oldminw
# example with weights.==1 and denominator.==1
    dtmtable, sett, df_prepped = prepare_dataframe_for_dtm!(df_tmp, treat_as_categorical_variable=["PLZ_WOHNORT"], numcol="LOSS20HALF", independent_vars=selected_explanatory_vars);
   

# multithreaded runs
    if false  # currently disabled
        @info("DTM:testing multithreading")
        @test Distributed.nprocs() == 1 # we generally expect that tests are run on only 1 thread
        dtmtable, sett, df_prepped = prepare_dataframe_for_dtm!(df_tmp, treat_as_categorical_variable=["PLZ_WOHNORT"], weightcol="EXPOSURE", numcol="LOSS20HALF", denomcol="PREMIUM66", independent_vars=selected_explanatory_vars);
        sett.model_type = "boosted_tree"
        if Distributed.nprocs() < 2
            Distributed.addprocs(2)
            Distributed.@everywhere using Decisiontrees
        end
        statsdf, settsdf, cvModels = dtm(dtmtable, sett, cvsampler)
    # same for multirun 
        a, b, allmodels = dtm(dtmtable, settV)
        dtm(dtmtable, settVBoosting)
    # remove workers again
        Distributed.rmprocs(workers())
    end

# coverage test of absrel
    DecisionTrees.absrel(rand(3), rand(3))

end # testset
