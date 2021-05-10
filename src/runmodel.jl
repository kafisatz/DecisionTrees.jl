
"""       
    prepare_dataframe_for_dtm!(dfin::DataFrame;directory::String=mktempdir(),treat_as_categorical_variable::Vector{String}=Vector{String}(),numcol::String="",denomcol::String="",weightcol::String="",trnvalcol::String="",valpct::Float64=0.3,keycol::String="",independent_vars::Vector{String}=Vector{String}())

Prepares a DataFrame for the modelling.

returns dtmtable::DTMTable,sett::ModelSettings,df_prepped::DataFrame
"""
function prepare_dataframe_for_dtm!(dfin::DataFrame;directory::String=mktempdir(),treat_as_categorical_variable::Vector{String}=Vector{String}(),numcol::String="",denomcol::String="",weightcol::String="",trnvalcol::String="",valpct::Float64=0.3,keycol::String="",independent_vars::Vector{String}=Vector{String}(),methodForSplittingPointsSelection::String="basedOnWeightVector")
    @assert isdir(directory)
    @assert in(methodForSplittingPointsSelection, globalConstAllowableMethodsForDefineCandidates)
	@time (df_prepped, sett) = prepareDF!(dfin, treat_as_categorical_variable=treat_as_categorical_variable, numcol=numcol, denomcol=denomcol, weightcol=weightcol, trnvalcol=trnvalcol, valpct=valpct, keycol=keycol, independent_vars=independent_vars);
	# set this to false
		sett.saveJLDFile = false
		
	fn = joinpath(directory, "DTMResult.csv")
	dtmtable = prep_data_from_df(df_prepped, sett, fn, methodForSplittingPointsSelection=methodForSplittingPointsSelection)
	return dtmtable, sett, df_prepped
end

function prepareDF!(dfin::DataFrame;treat_as_categorical_variable::Vector{String}=Vector{String}(),numcol::String="",denomcol::String="",weightcol::String="",trnvalcol::String="",valpct::Float64=0.3,keycol::String="",independent_vars::Vector{String}=Vector{String}())
    sz = size(dfin, 1)
    dfnames = propertynames(dfin)

    @assert ((valpct < 1) && (valpct > 0))
    @assert length(numcol) > 0 "Error: No argument numcol provided."
    @assert in(Symbol(numcol), dfnames) "Numerator column not found. We searched for the column $(numcol) in the vector $(string.(dfnames))"
    numDA = dfin[!,Symbol(numcol)]

    if denomcol == ""
        @info "Using constant denominator for each row (.==1)" # of the form ones(size(yourIntputDataFrame,1))"        
        denomDA = ones(sz)        
    else
        @assert in(Symbol(denomcol), dfnames) "Denominator column not found. We searched for the column $(denomcol) in the vector $(string.(dfnames))"
        denomDA = dfin[!,Symbol(denomcol)]
    end

    if weightcol == ""
        @info "Using constant weight for each row (.==1)" # of the form ones(size(yourIntputDataFrame,1))"
        weightDA = ones(sz)
    else
        @assert in(Symbol(weightcol), dfnames) "Weight column not found. We searched for the column $(weightcol) in the vector $(string.(dfnames))"
        weightDA = dfin[!,Symbol(weightcol)]
    end

    if trnvalcol == ""
        @info "Using random choice of Training-Validation column with proportion $(valpct) for validation"
        trnvalDA = map(x->x > valpct ? 1.0 : 0.0, rand(sz))
    else
        @assert in(Symbol(trnvalcol), dfnames) "Training-Validation column not found. We searched for the column $(trnvalcol) in the vector $(string.(dfnames))"
        trnvalDA = dfin[!,Symbol(trnvalcol)]
        # tbd ensure trnvalDA has only 0 and 1 entries
        tmp_found_vals = float.(sort(unique(trnvalDA))[:])
        tmp_expected_vals = Float64[0 1][:]
        if !isequal(tmp_found_vals, tmp_expected_vals)
            @show tmp_found_vals
            @assert false  "Error: Training-Validation column contained entries which are not equal to 1 or 0. Abort."
        end
        @assert isequal(tmp_found_vals, tmp_expected_vals)
    end

    if keycol == ""
        @info "Using canonical choice of identifier column 1:size(yourIntputDataFrame,1)"
        irkeyDA = collect(1:sz)
    else
        @assert in(Symbol(keycol), dfnames) "Identifier column not found. We searched for the column $(keycol) in the vector $(string.(dfnames))"
        irkeyDA = dfin[!,Symbol(keycol)]
    end

# remove union types (especially missing if there are any!)
    removeUnionTypes!(dfin, independent_vars)
    
# construct resulting DF
# requirements on the first 5 columns: irkey (=string identifier), numerator, denominator, weight, training_validation_bool (0 or 1)
    # 'main' columns
    dfres = DataFrame()
    dfres[!,:irkey] = irkeyDA
    dfres[!,:numerator] = numDA
    dfres[!,:denominator] = denomDA
    dfres[!,:weight] = weightDA
    dfres[!,:training_validation_bool] = trnvalDA

# explanatory columns
    if length(independent_vars) == 0
        independent_vars = string.(dfnames)
    end
    # delete main columns if they appear in the indep columns
    for remove_this_col in [numcol denomcol weightcol trnvalcol keycol][:]
        deleteat!(independent_vars, findall((in)([remove_this_col]), independent_vars))
    end
    n_indep = length(independent_vars)

    for x in Symbol.(treat_as_categorical_variable)        
        if !in(x, dfnames)
            @show x
            @show dfnames
            error("DTM: Column $(x) not found in dfnames. You provided that column as part of treat_as_categorical_variable.")            
        end
    end     
    @assert issubset(Symbol.(treat_as_categorical_variable), dfnames)
    # n_char_features=0
    char_features = Vector{Symbol}()
    num_features = Vector{Symbol}()
    for nn = 1:n_indep
        this_is_a_string = false
        this_name = Symbol(independent_vars[nn])
        @assert in(Symbol(this_name), dfnames) "Error: explanatory column $(this_name) not found in data."
        if in(this_name, Symbol.(treat_as_categorical_variable))
            # change variable type in original dataframe
            if !is_categorical_column(dfin[!,this_name], this_name)
                @info "DTM: Converting the column $(this_name) from $(eltype(dfin[!,this_name])) to String."
                dfin[!,this_name] = string.(dfin[!,this_name])
                this_is_a_string = true
            else
                # @info "DTM: Variable $(this_name) is already categorical in the input data. No type conversion performed."
            end
        end
        if is_categorical_column(dfin[!,this_name], this_name)
            push!(char_features, this_name)
        else
            push!(num_features, this_name)
        end

    end

	@assert (length(char_features) + length(num_features) == n_indep) # this should always be the case!
	n_char = length(char_features)
	n_num = length(num_features)

  # @warn("this might benefit from optimization. I am not entirely sure, how the code performs for larger data sets")
    # if n_num>0        
     #   @show num_features
     # !WARN this also sets the ORDER of the COLUMNS, which is imporant (unfortunately....). numeric columns need to come first.
    for x in vcat(num_features, char_features)
        dfres[!, x] = dfin[!, x] 
    end
        # dfres[num_features]=dfin[!,num_features]
    # end
    # if n_char>0        
        # dfres[char_features]=dfin[!,char_features]
    # end
    
    this_sett = ModelSettings()

    this_sett.df_name_vector = string.(propertynames(dfres)[1 + global_const_shift_cols:end])
    updateSettingsMod!(this_sett, number_of_num_features=n_num, number_of_char_features=n_char) 
    @info "Data initialization finished:"
    @show size(dfres)
    @show n_num, n_char
    @show char_features
    @show num_features
    return dfres, this_sett
end


"""
tries to convert all types of Union{Missing,T} to T
"""
function removeUnionTypes!(dfin, independent_vars::Vector{String})
    for v in independent_vars
        vsymb = Symbol(v)
        col = dfin[!,vsymb]
        elt = eltype(col)
        typeOfelt = typeof(elt)
        if typeOfelt != DataType
            @info "DTM: Type of the eltype of $(v) is not a DataType but $(typeOfelt). The type of $(v) is $(elt)"        
            @info "DTM: Trying to convert the type..."        
            try
                a = elt.a
                b = elt.b
                tragetT = ifelse(a == Missing, b, a)
                dfin[vsymb] = convert(Vector{tragetT}, x) 
                return nothing
            catch
                @warn("DTM: Type conversion failed. The model will likely not run as intended.")
            end
        end
    end
    return nothing 
end

function is_categorical_column(x, nm)
    # @warn("BK need to ensure that there are NO union types!, especially we do not want the missing tpiye (for now)!")
    elt = eltype(x)
    typeOfElt = typeof(elt)
    if typeOfElt != DataType  
        # @show x[1:10]
        error("DTM: The type of variable $(nm) is not of type 'DataType' but '$(typeOfElt)'. The type of the variable is $(elt). Currently only DataTypes are supported.")
    end
    return elt <: AbstractString
end

function prep_data_from_df(df_userinput::DataFrame, sett::ModelSettings, fn_with_ext::String;methodForSplittingPointsSelection::String="basedOnWeightVector")
    datafilename, ext = splitext(fn_with_ext) # actually one can provide fn_with_ext="c:\\temp\\my folder\\out" (so no extension is necessary)
    datafolder, outfileStringOnly = splitdir(datafilename)
    if length(outfileStringOnly) == 0
        error("DTM: fn_with_ext should be of the form C:\\some\\folder\\my_output.txt (the extension is irrelevant) \n You provided: $(fn_with_ext) \nAbort.")
    end
    dfIndata = df_userinput
    println("Preparing Data for analysis...")

    @assert check_for_missing_data(dfIndata, global_const_shift_cols) "Error dataset contains missing (NA) values! ABORT."
   
    key = convert(Vector{String}, string.(dfIndata[!,:irkey]))
    trn_val_idx = convert(Vector{UInt8}, UInt8.(dfIndata[!,:training_validation_bool]))
    @assert sort(unique(trn_val_idx)) == [0x00,0x01] "DTM: training validation column contains values not equal to 0 or 1." 
    # construct data
    numerator = convert(Array{Float64,1}, dfIndata[!,:numerator])
    denominator = convert(Array{Float64,1}, dfIndata[!,:denominator])
    weight = convert(Array{Float64,1}, dfIndata[!,:weight])
        
    count_neg_denom = sum(denominator .<= 0)
    if count_neg_denom > 0 
        @warn("DTM: Denominator contains observations with values <= 0.0")
        tot_sz = length(denominator)
        @warn("DTM: Negative rowcount: $(count_neg_denom) of $(tot_sz). I.e. a proportion of $(count_neg_denom / tot_sz)")
    end
    if any(denominator .== 0) 
        @warn("denominator[trn] contains observations with values 0.0")
      # foundobs=findfirst(denominator,0.0)
        foundobs = something(findfirst(isequal(0.0), denominator), 0)
        foundk = key[foundobs]
        println("First observation: key=$(foundk), rownumber_trn=$(foundobs), value=$(denominator[foundobs])\r\nDenominator needs to be nonzero (negative values are allowed)\r\nConsider discarding or changing these values/rows.")
        ignoreZeroDenominatorValues = sett.ignoreZeroDenominatorValues
        if !ignoreZeroDenominatorValues
            error("DTM: Abort.")
        else 
            @warn("DTM: Continuing data prep and modelling with zero denomintor observations. This is EXPERIMENTAL. Algorithm may fail.")
        end
    end
    strNegativeRatioWarning = "Negative Ratios may lead to errors during the modelling process for models which rely on multiplicative residuals! \nIt is suggested that you adjust the data such that no negative values occur!"
    if any(x->x < 0.0, denominator)
        @warn("DTM: There are training obsevations with negative values in the denominator column!")
        @warn(strNegativeRatioWarning)
    end
    # obsratiotrn=numeratortrn./denominatortrn
    if any(numerator .< 0)
        @warn("DTM: There are negative numerator values!")
        @warn(strNegativeRatioWarning)
        if sett.calculatePoissonError
            @warn("DTM: Setting calculatePoissonError to false!")
            sett.calculatePoissonError = false
        end        
        if typeof(sett.crit) == PoissonDevianceSplit
            error("DTM: Negative numerator values are not allowed for Poisson Deviance! Abort. \n Crit = $(sett.critt)")
        end
    end

    df_names = sett.df_name_vector # propertynames(dfIndata)
    # this DataFrame will contain both categorical and numerical features
    features = DataFrame() 
    
    print("Preparing numeric variables...")
    # this function modifies features -> the num features are added as pdas!
    candMatWOMaxValues = add_coded_numdata!(dfIndata, sett, trn_val_idx, sett.maxSplittingPoints, features, weight, methodForSplittingPointsSelection)
    println(" done.")
    # IMPORTANT convention  numerical data 'comes first' (i.e. first x columns are numerical), then categorical
    print("Preparing character variables...")
    for i = 1:sett.number_of_char_features
        if !(eltype(dfIndata[!,global_const_shift_cols + i + sett.number_of_num_features]) <: AbstractString) 
            @show eltype(dfIndata[!,global_const_shift_cols + i + sett.number_of_num_features])
            @show dfIndata[1:min(10, size(dfIndata, 1)),global_const_shift_cols + i + sett.number_of_num_features]
            error("DTM: Character variable $(i) = $(sett.df_name_vector[i + sett.number_of_num_features]) is not <:AbstractString in the dataframe.")
            # @warn("It may be that the variable was exported as character (by SAS) even though it is infact an integer. Please check!")            
        end        
        try
            colnumber = global_const_shift_cols + sett.number_of_num_features + i
            thisname = Symbol(df_names[i + sett.number_of_num_features]) 
            thiscol_as_utf8 = dfIndata[!,thisname] 
            if eltype(thiscol_as_utf8) != String
                thiscol_as_utf8 = convert(Array{String,1}, thiscol_as_utf8)
            end
            features[!,thisname] = DecisionTrees.PooledArraysDTM.PooledArray(thiscol_as_utf8)
            if length(features[!,thisname].pool) > 255
                @show nlevels = length(features[!,thisname].pool)
                @warn("DTM: Variable $(string(thisname)) has more than 255 levels, namely $(nlevels)\nThe tree may take very long to be constructed.")            
            end
        catch this_error
            @show this_error
            error("DTM: Unable to initialize character variable $(i) = $(sett.df_name_vector[i + sett.number_of_num_features]) (see error message above).")
        end      
    end
    mappings = deepcopy(map(i->features[!,i].pool, sett.number_of_num_features + 1:size(features, 2)))
    
    println(" done.")
    wlo, whi = extrema(weight)
    if wlo != whi
        print("Max weight is $(whi), min weight is $(wlo)")
        # @info("Currently the automatic choice of splitting points does not consider the weight distribution!")
    end
    
    pools = map(i->features[!,i].pool, 1:size(features, 2)) 
    pool_lengths = length.(pools)
    for i = 1:length(pool_lengths)
        if pool_lengths[i] == 1
            @info "DTM: Column $(i) = $(sett.df_name_vector[i]) has zero splitting points."
        end
    end
        
    # "Check if string pools are valid utf8 characters. This is specifically important when data is read in a 'wrong' manner due to encoding issues. The PyCall functions which write the Excel files will fail in such cases"
    for i = 1:length(pool_lengths)
    # @show pools        
        if eltype(pools[i]) <: Number
            # nothing to do, numbers are always valid
        else 
            vld = .!(isvalid.(pools[i]))
            # vld=isvalid.(pools[i])
            if any(vld)
                @warn("DTM: The variable $(sett.df_name_vector[i]) contains invalid characters for some values (see below). Please check if your data was read and encoded properly.\r\nDTM will try and remove the invalid characters from the data. You should clean the data prior to invoking DTM.")
                @show pools[i][vld]
                for zz = 1:length(vld)
                    if vld[zz]
                        before = pools[i][zz]
                        after = cleanString(pools[i][zz])
                        println("DTM: $(before) -> $(after)")
                        pools[i][zz] = after
                    end
                end
            end
        end
    end

    trnidx, validx = createIntegerIndices(trn_val_idx)
    dtmtable = DTMTable(key, trnidx, validx, numerator, denominator, weight, features, candMatWOMaxValues, mappings)
    
    if sett.saveJLDFile
        jldfile = string(fileroot(datafilename), ".jld2")
        println("Saving prepped data to file:\n $(jldfile)")
        isfile(jldfile) && rm(jldfile)
        # @time save(jldfile,"dtmtable",dtmtable,"candMatWOMaxValues",candMatWOMaxValues,"features",features,"mappings",mappings,"key",key,"trn_val_idx",trn_val_idx,"numerator",numerator,"denominator",denominator,"weight",weight,"oldsettings",sett)
        @time save(jldfile, "dtmtable", dtmtable, "candMatWOMaxValues", candMatWOMaxValues, "mappings", mappings, "sett", sett)
        println("JLD2 file saved. Starting modelling... \n")
    end

    return dtmtable
end

"""
This is the core modelling function of DecisionTrees
    
    run_model_actual(dtmtable::DTMTable,input_setttings::ModelSettings,fn::String)
"""
function run_model_actual(dtmtable::DTMTable, input_setttings::ModelSettings, fn::String)
# this function is called with prepped julia data, this is the core modelling function (all previous ones are for preparational tasks only)
    t0RunModel0 = time_ns()
    sett = deepcopy(input_setttings)
    @assert sett.chosen_apply_tree_fn == "apply_tree_by_leaf" # currently this is the default choice. Consider the code in boosting.jl and bagging.jl -> the apply tree fn is currently hardcoded
    if sett.calculateGini
        @warn("DTM: Note that Gini and RSS metrics are experimental and may need review. Also it needs to be clarified if the metrics are evaluated only for the numerator or for the ratio=numerator/denominator.")
    end
    x = checkIfSettingsAreValid(input_setttings)
    if !(x == nothing)
        @warn "DTM: You may have provided invalid or inconsistent model settings. The model run may fail."
    end
    if (sett.fitForStatsAndCharts == "smoothedPerScore") && (sett.smoothEstimates != "1")
        @info("DTM: setting sett.smoothEstimates to 1, because sett.fitForStatsAndCharts=$(sett.fitForStatsAndCharts)")    
        sett.smoothEstimates = "1"
    end
    if sett.writeIterationMatrix
        if !sett.prroduceEstAndLeafMatrices
            @info "DTM: writeIterationMatrix=$(sett.writeIterationMatrix). The estimate matrices are required for this. Setting sett.prroduceEstAndLeafMatrices to true."
        end
    end
    if sett.boolNumeratorStats
        if !sett.prroduceEstAndLeafMatrices
            @info "DTM: User has requested NumeratorStats (boolNumeratorStats=$(sett.boolNumeratorStats)). The estimate matrices are required for these statistics. Setting prroduceEstAndLeafMatrices to true."
            sett.prroduceEstAndLeafMatrices = true
        # @show sett.boolNumeratorStats,sett.prroduceEstAndLeafMatrices
        end
    end

    if (typeof(sett.crit) == mseSplit) || (typeof(sett.crit) == sseSplit) || (typeof(sett.crit) == NormalDevianceDifferenceToMeanFitWEIGHTEDSplit)
        if any(iszero, dtmtable.denominator) # must be non zero (ideally even positive)
            error("DTM: Denominator has zero values. This is not allowed for crit=$(sett.crit)")
        end
        if any(dtmtable.denominator .< 0) # must be non zero (ideally even positive)
            @warn("DTM: Denominator has negative values. You have specified crit=$(sett.crit). This may lead to unintended results as the algorithm divides by the denominator (and it is used as weight!)")
        end
    end

    path_and_fn_wo_extension, ext = splitext(fn)
# check for file extension. This does not seem to be required... (then again we should ensure that the user provides a filename stump and not only a folder name)
    if in(lowercase(ext), ["",".csv",".txt",".sas",".jl",".pdf"])
    else
        @show fn
        @show path_and_fn_wo_extension, ext
        error("DTM: Unexpected file name extension $(ext). You provided the $(fn). Please provide a file name extension such as .txt")    
    end
    datafolder, outfileStringOnly = splitdir(path_and_fn_wo_extension)
    @assert isdir(datafolder) "DTM: The specified output folder $(datafolder) does not exist. Abort."

    filelistWithFilesToBeZipped = String[]
    csharp_code_file = string(path_and_fn_wo_extension, ".cs")
    vba_code_file = string(datafolder, "\\", outfileStringOnly, ".vba.txt")
    sas_code_file = string(path_and_fn_wo_extension, ".sas")
    tree_file = string(path_and_fn_wo_extension, ".txt")
    dot_graph_TXT_file = string(path_and_fn_wo_extension, ".dot.txt")
    dot_graph_visualization_file = string(path_and_fn_wo_extension, ".dot.pdf")
    this_outfile = convert(String, string(path_and_fn_wo_extension, ".csv"))
    statsfile = string(path_and_fn_wo_extension, ".stats.csv")
    jldresultsfile = string(path_and_fn_wo_extension, ".jld2")
    statsfileExcel = string(path_and_fn_wo_extension, ".xlsx")

# for now trnindx and validx should always be sorted!
    if !issorted(dtmtable.trnidx) 
        @info "DTM: trnidx was not sorted. Sorting trnindex"
        sort!(dtmtable.trnidx)
    end
    if !issorted(dtmtable.validx) 
        @info "DTM: validx was not sorted. Sorting validx"
        sort!(dtmtable.validx)
    end

    key = dtmtable.key
    trnidx = dtmtable.trnidx
    validx = dtmtable.validx
    features = dtmtable.features
    weight = dtmtable.weight
    numerator = dtmtable.numerator
    denominator = dtmtable.denominator

# trnidx_one_zero_full_length is added to the output file
    @assert issorted(trnidx)
    trnidx_one_zero_full_length = map(x->length(searchsorted(trnidx, x)), 1:length(key))

# define 'version', we store the sha1 hash of the current commit in sett.version
    sett.version = "--" # get_sha1() #need to fix this/ find out how to read package hash in Julia code

# initialize random number generator
    # we may want to make this independent from the data (trn), then again different data should/will lead to a different result.
# hash of DF is not working properly, thus we need to iteration over the columns of the DF (this can be simpplified later on.)
    s00 = hash(trnidx, hash(validx, hash(sett.seed, hash(numerator, hash(denominator, hash(weight))))))
    for x = 1:size(features, 2)
        s00 = hash(features[!,x], s00)
    end
    srandInt = floor(Int, 1 / 3 * hash(931, s00))
    Random.seed!(srandInt)

    if sett.subsampling_prop < 1.0 
        @warn("BK: subsampling_prop <1. This might be slow in the current implementation!")
    end

# there are at least two ways to derive estimates for hold out data: either loop through the leaves of the tree or loop over the rows of the validation data setting
# also there are two types of hold out estimates: either the ones derived from training data (more meaningful case) or the ones derived from validation labels
    sumwtrn = sum(view(weight, trnidx))
    sumwval = sum(view(weight, validx))
    minweighttmp = sett.minWeight < 0.0 ? sumwtrn * (-max(-1.0, sett.minWeight)) : sett.minWeight
    @assert minweighttmp >= 0
    sumwtrn >= 2 * minweighttmp || @warn("DTM: Specified minweight is larger than half of the total training weight. No split is possible!")

    someInteger = sumwtrn / minweighttmp
    if someInteger > 100000
        @warn("DTM: Total training weight is: $(sumwtrn) and you specified minWeight=$(minweighttmp). The tree might become very large")
    end

    prnt = sett.print_details
    general_settings = convert(String, string("Write tree to txt file: $(sett.writeTree), statsByVariables=$(join(sett.statsByVariables, ','))"))
# general_settings=convert(String,string(general_settings,))
    char_levels = map(x->length(features[!,x].pool), sett.number_of_num_features + 1:sett.number_of_num_features + sett.number_of_char_features)
    num_levels = map(x->length(features[!,x].pool), 1:sett.number_of_num_features)
    if length(char_levels) > 0
        max_c0 = maximum(char_levels)
    else
        max_c0 = 0
    end
    if length(num_levels) > 0
        max_n0 = maximum(num_levels)
    else
        max_n0 = 0
    end
    max_levels_uint8 = typemax(UInt8)
    if max(max_c0, max_n0) > max_levels_uint8
        @info("DTM: At least one variable has more than $(max_levels_uint8) potential splitting points. \r\nThis is experimental and not yet fully tested.")        
    end    
    @assert max_c0 < nLevelsThreshold "A variable has too many levels (i.e. more than $(nLevelsThreshold)). You can try and increase this threshold but the algorithm may take quite long for the modelling"
    @assert max_n0 < nLevelsThreshold "A variable has too many levels (i.e. more than $(nLevelsThreshold)). You can try and increase this threshold but the algorithm may take quite long for the modelling"
    if size(sett.moderationvector, 1) > 1;if (sett.adaptiveLearningRate != 1.0 && sett.adaptiveLearningRate != 0.0);println("######################################## \nWarning: (bk) Adaptive Learning has been disabled as a moderation vector has been provided. \n######################################## ");end;sett.adaptiveLearningRate = 1.0;end;
    if prnt
        println("")
        println(stars)
        println("Time:\t\t", Dates.now())
        println(stars)
        println("---General Settings------------------------------------------------------------")
        println(general_settings)
        print_some_settings(sett, ["writeSasCode","writeIterationMatrix","writeResult","writeStatistics","boolCreateZipFile","writeCsharpCode","statsRandomByVariable","saveJLDFile"])
        println("---Data Summary----------------------------------------------------------------")
        println("Number of records - Trn: $(length(trnidx)) Val: $(length(validx)) Val Proportion: $(round(length(validx) / (length(trnidx) + length(validx)), sigdigits=4))")
        println("Weight - Trn: $(round(sumwtrn, sigdigits=6)) Val: $(round(sumwval, sigdigits=6)) Val Proportion: $(round(sumwval / (sumwtrn + sumwval), sigdigits=4))")
        println("Number of numeric/character variables: $(sett.number_of_num_features)/$(sett.number_of_char_features)")
    end
    if size(num_levels, 1) > 0
        prnt && println(string("Number of levels for numeric variables: (Maximal # of splitting points: $(sett.maxSplittingPoints))"))
        num_level_string = join(sort(num_levels, rev=true), ",")
        prnt && println(string(" ", num_level_string))
    else
        prnt && println("The data does not contain any numeric attributes.")
        num_level_string = ""
    end
    if size(char_levels, 1) > 0
        prnt && println(string("Number of levels for character variables: "))
        # println(arr_to_str(sort(char_levels, rev=true)))
        char_level_string = join(sort(char_levels, rev=true), ",")
        prnt && println(string(" ", char_level_string))
    else
        char_level_string = ""
        prnt && println("The data does not contain any character attributes.")
    end
    prnt && println("catSortByThreshold: $(sett.catSortByThreshold), catSortBy: $(replace(string(sett.catSortBy), "DecisionTrees." => "")), smoothEstimates: $(sett.smoothEstimates), deriveFitPerScoreFromObservedRatios: $(sett.deriveFitPerScoreFromObservedRatios)")

  # create a string which contains all model settings (this will be written to output file(s))
    model_setting_string = convert(String, string("Number of records - Trn: $(length(trnidx)) Val: $(sum(validx)) \nNumber of levels for numeric variables: (Maximal # of splitting points: $(sett.maxSplittingPoints)) "))
    model_setting_string = convert(String, string(model_setting_string, "\n", num_level_string, "\nNumber of levels for character variables: \n", char_level_string))
    model_setting_string = convert(String, string(model_setting_string, "\n", general_settings))
    model_setting_string = convert(String, string(model_setting_string, "\n", "Type: $(sett.model_type) \nMinweight: $(sett.minWeight) \nRnd: $(sett.randomw) \nNumMaxSplitPoints: $(sett.maxSplittingPoints) \nCriterion: $(replace(string(sett.crit), "DecisionTrees." => "")) \nIterations: $(sett.iterations) \nModeration Vector: $(join(sett.moderationvector, ",")) \nNScores: $(sett.nScores) \nBoolStartAtMean: $(sett.startAtMean) "))
    model_setting_string = convert(String, string(model_setting_string, "\n", "Datafolder= \n", datafolder, "\nOutfilename= \n", outfileStringOnly, " "))
    model_setting_string = convert(String, string(model_setting_string, "\n\n", "CSV column names (in the order in which Julia processes them):\n", join(sett.df_name_vector, ",")))
      # permutedims is "transpose" in julia 0.5
      # model_setting_string=convert(String,string(model_setting_string, "\n\n", "Names:Levels: \n",join(permutedims(names_and_levels, (2, 1)),":")))
    model_setting_string = convert(String, string(model_setting_string, "\n", "Folder of Julia code (decisiontree.jl): \n", dirname(@__FILE__)))
      # model_setting_string=convert(String,string(model_setting_string, "\n", "Julia started at: \n not_available "))
    model_setting_string = convert(String, string(model_setting_string, "\n", "Julia (Modelling) start time: \n$(Dates.now()) "))

    prnt && println("---Model Settings--------------------------------------------------------------")
    if sett.boolRankOptimization && "build_tree" == sett.model_type
        error("Error: (bk) I do not think that _RankOpt is currently working with a single tree....")
    end

    local resulting_model # ,resultEnsemble
    # resulting_model=EmptyEnsemble()
    # local xlData,estimatedRatio,vectorOfLeafNumbers,vectorOfLeafArrays,rawObservedRatioPerScore,est_matrixFromScores,stats,estimateUnsmoothed,estimateSmoothed,estimateFromRelativities,estimatesPerScoreForCode

    if "build_tree" == sett.model_type
        if prnt
            println("Type: $(sett.model_type), Minweight: $(sett.minWeight), Rnd: $(sett.randomw)")
            println("NumMaxSplitPoints: $(sett.maxSplittingPoints), Criterion: $(replace(string(sett.crit), "DecisionTrees." => ""))")
            println(stars)
            println(string(""))
            println(string("Building Tree... Time: $(Dates.now())"))
        end
        obsvalcols = 1
        if sett.randomw > 0
            @warn("Running single tree with positive random weight:$(sett.randomw)")
        end

    # build a simple tree
        fitted_values_tree = zero(numerator)    
        T_Uint8_or_UInt16 = find_max_type(features)
        if prnt
            @time tree = build_tree!(trnidx, validx, dtmtable.candMatWOMaxValues, dtmtable.mappings, sett, numerator, denominator, weight, features, fitted_values_tree, T_Uint8_or_UInt16)
        else
            tree = build_tree!(trnidx, validx, dtmtable.candMatWOMaxValues, dtmtable.mappings, sett, numerator, denominator, weight, features, fitted_values_tree, T_Uint8_or_UInt16)
        end
        
        resulting_model = tree
        prnt && println("Deriving Trn Estimates... Time: $(Dates.now())")
        leaves_of_tree = create_leaves_array(tree)
        leaf_number_vector = zeros(Int, length(fitted_values_tree))
        # todo/tbd possibly improve this by only applying the function once!
        # notably fitted_values_tree is overwritten here, however the values do NOT change.
        apply_tree_by_leaf!(fitted_values_tree, leaf_number_vector, trnidx, tree.rootnode, features)
        if length(validx) > 0
            apply_tree_by_leaf!(fitted_values_tree, leaf_number_vector, validx, tree.rootnode, features)
        end        
        prnt && println("Preparing output... Time: $(Dates.now())")
         # 1-dim Variable Importance
        var_imp1d_str_arr, var_imp2d_str_arr, onedimintvec_unused, twodimintvec_unused = variable_importance(leaves_of_tree, sett.df_name_vector, sett.number_of_num_features)
        nleaves = Int(size(leaves_of_tree, 1))
        res = hcat(key, trnidx_one_zero_full_length, fitted_values_tree)
        # Create output statistics        
            # this is not done in an efficient manner yet..... we are copying data to get est[trnidx] etc...
        @time xlData, overallstats = createSingleTreeExcel(trnidx, validx, sett, tree, leaves_of_tree, fitted_values_tree, leaf_number_vector, numerator, denominator, weight)
        resulting_model.exceldata = xlData
        if sett.writeStatistics
            writeStatisticsFile!(statsfileExcel, xlData, filelistWithFilesToBeZipped)
        end
        z = permutedims(overallstats, (2, 1))        
        dfStats = DataFrame(convert(Array{Float64,2}, z[2:size(z, 1),:]),Symbol[Symbol(x) for x in view(z, 1, :)])
        DataFrames.rename!(dfStats, Symbol[Symbol(x) for x in view(z, 1, :)])
        resulting_model.modelstats = dfStats
        if sett.saveResultAsJLDFile
            println("Saving results to jld2 file: \n $(jldresultsfile)")
            isfile(jldresultsfile) && rm(jldresultsfile)
            @time JLD2.@save jldresultsfile xlData res leaves_of_tree trnidx validx fitted_values_tree tree var_imp1d_str_arr var_imp2d_str_arr onedimintvec_unused twodimintvec_unused
            push!(filelistWithFilesToBeZipped, jldresultsfile)
        end

         # print_tree(tree)
        dot_graph = graph(resulting_model)
        model_setting_string = convert(String, string(model_setting_string, "\n", "Julia end time: \n$(Dates.now()) \n"))
        if sett.write_dot_graph
            println("Writing DOT tree structure to file: \n $(dot_graph_TXT_file)")
            isfile(dot_graph_TXT_file) && rm(dot_graph_TXT_file)
            f = open(dot_graph_TXT_file, "w");
            write(f, "\n\n", dot_graph, "\n\n");
            close(f);
            draw_dot_graph(sett.graphvizexecutable, dot_graph_TXT_file, dot_graph_visualization_file)
            push!(filelistWithFilesToBeZipped, dot_graph_TXT_file)
            isfile(dot_graph_visualization_file) && push!(filelistWithFilesToBeZipped, dot_graph_visualization_file)
        end
        if sett.writeTree
            # write settings
            println("Writing tree structure to file: \n $(tree_file)")                
            isfile(tree_file) && rm(tree_file)
            f = open(tree_file, "w")
            write(f, model_setting_string)
            write(f, "\n\n", dot_graph, "\n\n")
            close(f);
            write_tree(dtmtable.candMatWOMaxValues, tree, sett.number_of_num_features, var_imp1d_str_arr, var_imp2d_str_arr, 0, tree_file, sett.df_name_vector, dtmtable.mappings)
            push!(filelistWithFilesToBeZipped, tree_file)            
        end
        if sett.writeSasCode
                # write SAS Code
            if isa(tree.rootnode, Leaf)
                @warn("DTM: Tree is a single leaf. No SAS Code was produced.")
            else 
                println("Writing SAS Code: \n $(sas_code_file)")
                isfile(sas_code_file) && rm(sas_code_file)
                writeSasCode(dtmtable.candMatWOMaxValues, tree, sett.number_of_num_features, sas_code_file, sett.df_name_vector, model_setting_string, dtmtable.mappings, 1.0)
                push!(filelistWithFilesToBeZipped, sas_code_file)
            end
        end        
        if sett.writeResult
            println("Exporting Data: \n $(this_outfile)")
            isfile(this_outfile) && rm(this_outfile)
            DelimitedFiles.writedlm(this_outfile, res, ',')
            push!(filelistWithFilesToBeZipped, this_outfile)
        end
        # return filelistWithFilesToBeZipped,resulting_model
    end  # end build tree
    
    # output for boosting & bagging
    if in(sett.model_type, ["boosted_tree","bagged_tree"])
    
        if "boosted_tree" == sett.model_type
            @timeConditional(sett.print_details,begin
                xlData, estimatedRatio, vectorOfLeafNumbers, vectorOfLeafArrays, rawObservedRatioPerScore, est_matrixFromScores, stats, estimateUnsmoothed, estimateSmoothed, estimateFromRelativities, resulting_model = boosted_tree(dtmtable, sett)
            end)
            model_setting_string = convert(String, string(model_setting_string, "\n", "Julia end time: \n$(Dates.now()) \n "))
            # this line is for  the sas,vba and csharp code
            estimatesPerScoreForCode = sett.smoothEstimates == "0" ? rawObservedRatioPerScore : resulting_model.ScoreToSmoothedEstimate
        end
      
        if sett.writeSasCode
                # write SAS Code
            println("Writing SAS Code: \n $(sas_code_file)")
            isfile(sas_code_file) && rm(sas_code_file)
            if sett.roptForcedPremIncr
                warn("This does currently not work")
            else
                writeSasCode(estimatesPerScoreForCode, dtmtable.candMatWOMaxValues, resulting_model, sett.number_of_num_features, sas_code_file, sett.df_name_vector, model_setting_string, dtmtable.mappings)
            end
            push!(filelistWithFilesToBeZipped, sas_code_file)
        end
        if sett.writeCsharpCode                
            println("Writing C# Code: \n $(csharp_code_file)")
            isfile(csharp_code_file) && rm(csharp_code_file)                
            writeCsharpCode(vectorOfLeafArrays, estimatesPerScoreForCode, dtmtable.candMatWOMaxValues, resulting_model, csharp_code_file, model_setting_string, dtmtable.mappings, 0, sett)
            push!(filelistWithFilesToBeZipped, csharp_code_file)
        end 
        if sett.writeVbaCode
            println("Writing VBA Code: \n $(vba_code_file)")
            isfile(vba_code_file) && rm(vba_code_file)
            writeVbaCode(vectorOfLeafArrays, estimatesPerScoreForCode, dtmtable.candMatWOMaxValues, resulting_model, vba_code_file, model_setting_string, dtmtable.mappings, 0, sett)
            push!(filelistWithFilesToBeZipped, vba_code_file)
        end            
       # end #boosting
        #= 
            elseif "bagged_tree"==sett.model_type        
                #bagging
                @timeConditional(sett.print_details,begin
                @warn("bagging is not yet updated for dtm2! This will most likely fail!")
                    xlData,estimatedRatio,vectorOfLeafNumbers,vectorOfLeafArrays,rawObservedRatioPerScore,est_matrixFromScores,stats,estimateUnsmoothed,estimateSmoothed,estimateFromRelativities,resulting_model=bagged_tree(dtmtable,sett)
                end)            
                model_setting_string=convert(String,string(model_setting_string, "\n", "Julia end time: \n$(Dates.now()) \n "))
            #else
                #error("DTM: Unknown Model Type: $(sett.model_type) (Within Bagging/Boosting loop)")
            end #end bagging =#
        
    # after boosting bagging  export etc
    # here we are in the "section":     #if in(sett.model_type,["bagged_tree","boosted_tree"])
    # output for boosting & bagging#write tree

    # resulting_model=resultEnsemble

        if sett.writeResult
            res = hcat(key, trnidx_one_zero_full_length, resulting_model.scores, estimateUnsmoothed, estimateSmoothed, estimatedRatio)
        end    
    # save results as *.JLD2 file
        if sett.saveResultAsJLDFile
            println("Saving results to jld2 file: \n $(jldresultsfile)")
            isfile(jldresultsfile) && rm(jldresultsfile)
            @time JLD2.@save jldresultsfile xlData dtmtable trnidx validx vectorOfLeafNumbers vectorOfLeafArrays rawObservedRatioPerScore est_matrixFromScores stats estimateUnsmoothed estimateSmoothed estimateFromRelativities resulting_model
            push!(filelistWithFilesToBeZipped, jldresultsfile)
        end

        if sett.writeTree
                # write settings
            println("Writing tree structure to file: \n $(tree_file)")
            isfile(tree_file) && rm(tree_file)
            f = open(tree_file, "w");write(f, model_setting_string);close(f);
            write_tree(dtmtable.candMatWOMaxValues, resulting_model, sett.number_of_num_features, 0, tree_file, sett.df_name_vector, dtmtable.mappings)
            push!(filelistWithFilesToBeZipped, tree_file)
        end
    # DOT Graphs for boosted trees
        if sett.write_dot_graph
            println("Writing DOT tree structure to file: \n $(dot_graph_TXT_file)")
            fnbase,ext = splitext(dot_graph_TXT_file)
            for i=1:sett.iterations
                #txt filename
                dot_graph_TXT_file = string(fnbase,".tree.",i,ext)
                #dot graph pdf filename 
                dot_graph_visualization_file = string(fnbase,".tree.",i,".pdf")

                #write each tree
                tree = gettree(resulting_model,i)
                dot_graph = graph(tree)
                
                isfile(dot_graph_TXT_file) && rm(dot_graph_TXT_file)                
                f = open(dot_graph_TXT_file, "w");
                    write(f, "\n\n", dot_graph, "\n\n");
                close(f);
                #dot_graph_visualization_file = "C:\\Users\\BERNHA~1.KON\\AppData\\Local\\Temp\\jl_ta4Y2r\\dtmresult.dot.pdf"
                draw_dot_graph(sett.graphvizexecutable, dot_graph_TXT_file, dot_graph_visualization_file)
                push!(filelistWithFilesToBeZipped, dot_graph_TXT_file)
                isfile(dot_graph_visualization_file) && push!(filelistWithFilesToBeZipped, dot_graph_visualization_file)
            end
        end
    # Write statistics
        if sett.writeStatistics
            writeStatisticsFile!(statsfileExcel, xlData, filelistWithFilesToBeZipped)
        end
    # Write result to CSV
        if sett.writeResult
            println("Exporting Data: \n $(this_outfile)")
            isfile(this_outfile) && rm(this_outfile)
            DelimitedFiles.writedlm(this_outfile, res, ',')
            push!(filelistWithFilesToBeZipped, this_outfile)
        end
    # write estimates matrix
        if sett.writeIterationMatrix    && (length(vectorOfLeafNumbers) < length(key))
            @warn "DTM: The iteration matrices wil not be written to file. length(key)=$(length(key)), but length(vectorOfLeafNumbers)=$(length(vectorOfLeafNumbers)). \r\n"
            @warn "DTM: Please check the values of sett.writeIterationMatrix and sett.prroduceEstAndLeafMatrices"
            sett.prroduceEstAndLeafMatrices || @warn("DTM: You may want to set sett.prroduceEstAndLeafMatrices to true")
        else
            if sett.writeIterationMatrix        
                # Matrix with Leaf Numbers
                leafnrfile = string(path_and_fn_wo_extension, ".leafnumbers.csv")
                vectorOfLeafNumbersMOD = hcat(key, vectorOfLeafNumbers[:,2:end])
                println("Exporting LeafNumber Matrix: \n $(leafnrfile)")
                isfile(leafnrfile) && rm(leafnrfile)
                DelimitedFiles.writedlm(leafnrfile, vectorOfLeafNumbersMOD, ',')
                push!(filelistWithFilesToBeZipped, leafnrfile)
                # MATRIX 1
                estmat_outile2 = string(path_and_fn_wo_extension, "_iteration_matrixFromScores.csv")
                res_estmatrixFromScores = hcat(key, est_matrixFromScores)
                println("Exporting Estimates Matrix (based on Scores): \n $(estmat_outile2)")
                isfile(estmat_outile2) && rm(estmat_outile2)
                DelimitedFiles.writedlm(estmat_outile2, res_estmatrixFromScores, ',')
                push!(filelistWithFilesToBeZipped, estmat_outile2)
                # MATRIX 2
                estmat_outile = string(path_and_fn_wo_extension, "_iteration_matrix.csv")
                res_estmatrix = hcat(key, resulting_model.iterationmatrix)
                println("Exporting Estimates Matrix: \n $(estmat_outile)")                        
                isfile(estmat_outile) && rm(estmat_outile)
                DelimitedFiles.writedlm(estmat_outile, res_estmatrix, ',')
                push!(filelistWithFilesToBeZipped, estmat_outile)
            end
        end
      
    end # end of:  elseif in(sett.model_type,["boosted_tree","bagged_tree"])
# else
#    error("Unknown Model Type: $(sett.model_type)")
# end #end distinction between three model types

# this code applies to all model types
   # Create ZIP file
     # Push Log file to file list
    elapsed_until_this_point = (-t0RunModel0 + time_ns()) / 1e9
    prnt && println("Modelling finished. Time: $(Dates.now()) - Total time was $(round(elapsed_until_this_point, digits=1))s = $(round(elapsed_until_this_point / 60, digits=1))m")
    if sett.boolCreateZipFile
        logfile = string(path_and_fn_wo_extension, ".log")
        tmplogfile = string(path_and_fn_wo_extension, ".zip.log");isfile(tmplogfile) && rm(tmplogfile)
        try
            isfile(logfile) && isreadable(logfile) && cp(logfile, tmplogfile)
        catch error_while_copying_log
            @show error_while_copying_log
        end
        isfile(tmplogfile) && push!(filelistWithFilesToBeZipped, tmplogfile)
        zipFilename = string(path_and_fn_wo_extension, ".7z")
        println("Creating zip file with results \n $(zipFilename)")
        @time createZipFile(zipFilename, filelistWithFilesToBeZipped)
    end
  # result needs to have the same type, irrespective if we ran a tree or a boosted tree
  # return true

    return filelistWithFilesToBeZipped, resulting_model
end


function run_model(ARGS;nrows::Int=-1)
   # the first call of this function takes about 15 seconds (probably the first call of the module takes that much time....)
    failed_return_value = String["false"] # should be of the same type as the result of run_model
    try
        # first function which is called by run.jl:
        settingsFilename, dataFilename, datafolder, outfilename, outfileStringOnly = return_file_and_folders(ARGS)
        result_of_runmodel, resulting_model = run_model_main(settingsFilename, dataFilename, convert(String, datafolder), outfilename, outfileStringOnly, nrows=nrows)
        if typeof(result_of_runmodel) == Array{String,1}
            return result_of_runmodel
        else
            @warn("Unexpected Model result:")
            @show typeof(result_of_runmodel)
            @show result_of_runmodel
            dump(result_of_runmodel)
            return failed_return_value
        end
    catch model_error
        @show model_error
        @warn("An error occurred during modelling phase. See above.")
        return failed_return_value
    end
end

function dtm(dtmtable::DTMTable, sett::ModelSettings;file::String=joinpath(mktempdir(), defaultModelNameWtihCSVext))
    return run_model_actual(dtmtable::DTMTable, sett::ModelSettings, file::String)
end