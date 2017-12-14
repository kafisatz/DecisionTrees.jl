#consider the functions
#prepare_dataframe_for_dtm -> prepares an arbitrary df from the dtm models
#prepare_non_jld_data

function run_model(ARGS;nrows::Int=-1)
   #the first call of this function takes about 15 seconds (probably the first call of the module takes that much time....)
   failed_return_value=String["false"] #should be of the same type as the result of run_model
   try
		tic() #this tic is for the time measurment until right before the zip file is created.
		tic()		
		#first function which is called by run.jl:
		settingsFilename,dataFilename,datafolder,outfilename,outfileStringOnly=return_file_and_folders(ARGS)
		result_of_runmodel,resulting_model=run_model_main(settingsFilename,dataFilename,convert(String,datafolder),outfilename,outfileStringOnly,nrows=nrows)
		if typeof(result_of_runmodel)==Array{String,1}
			return result_of_runmodel
		else
			warn("Unexpected Model result:")
			@show typeof(result_of_runmodel)
			@show result_of_runmodel
			dump(result_of_runmodel)
			return failed_return_value
		end
   catch model_error
		@show model_error
		warn("An error occurred during modelling phase. See above.")
		return failed_return_value
   end
end

function prepare_dataframe_for_dtm!(dfin::DataFrame;treat_as_categorical_variable::Vector{String}=Vector{String}(),numcol::String="",denomcol::String="",weightcol::String="",trnvalcol::String="",valpct::Float64=0.3,keycol::String="",independent_vars::Vector{String}=Vector{String}())
    denomcol=uppercase(denomcol)
    keycol=uppercase(keycol)
    trnvalcol=uppercase(trnvalcol)
    weightcol=uppercase(weightcol)
    independent_vars=uppercase.(independent_vars)

    sz=size(dfin,1)
    dfnames=names(dfin)

    @assert ((valpct<1)&&(valpct>0))
    @assert length(numcol)>0 "Error: No argument numcol provided."
    @assert in(Symbol(numcol),dfnames) "Numerator column not found. We searched for the column $(numcol) in the vector $(string.(dfnames))"
    numDA=dfin[Symbol(numcol)]

    if denomcol==""
        info("Using constant denominator for each row (.==1)") #of the form ones(size(yourIntputDataFrame,1))")
        denomDA=ones(sz)
    else
        @assert in(Symbol(denomcol),dfnames) "Denominator column not found. We searched for the column $(denomcol) in the vector $(string.(dfnames))"
        denomDA=dfin[Symbol(denomcol)]
    end

    if weightcol==""
        info("Using constant weight for each row (.==1)") #of the form ones(size(yourIntputDataFrame,1))")
        weightDA=ones(sz)
    else
        @assert in(Symbol(weightcol),dfnames) "Weight column not found. We searched for the column $(weightcol) in the vector $(string.(dfnames))"
        weightDA=dfin[Symbol(weightcol)]
    end

    if trnvalcol==""
        info("Using random choice of Training-Validation column with porportion $(valpct) for validation")
        trnvalDA=map(x->x<valpct ? 1.0 : 0.0,rand(sz))
    else
        @assert in(Symbol(trnvalcol),dfnames) "Training-Validation column not found. We searched for the column $(trnvalcol) in the vector $(string.(dfnames))"
        trnvalDA=dfin[Symbol(trnvalcol)]
        #tbd ensure trnvalDA has only 0 and 1 entries
        tmp_found_vals=float.(sort(unique(trnvalDA))[:])
        tmp_expected_vals=Float64[0 1][:]
        if !isequal(tmp_found_vals,tmp_expected_vals)
            @show tmp_found_vals
            @assert false  "Error: Training-Validation column contained entries which are not equal to 1 or 0. Abort."
        end
        @assert isequal(tmp_found_vals,tmp_expected_vals)
    end

    if keycol==""
        info("Using canonical choice of indentifier column 1:size(yourIntputDataFrame,1)")
        irkeyDA=collect(1:sz)
    else
        @assert in(Symbol(keycol),dfnames) "Identifier column not found. We searched for the column $(keycol) in the vector $(string.(dfnames))"
        irkeyDA=dfin[Symbol(keycol)]
    end

#construct resulting DF
#requirements on the first 5 columns: irkey (=string identifier), numerator, denominator, weight, training_validation_bool (0 or 1)
    #'main' columns
    dfres=DataFrame()
    dfres[:irkey]=irkeyDA
    dfres[:numerator]=numDA
    dfres[:denominator]=denomDA
    dfres[:weight]=weightDA
    dfres[:training_validation_bool]=trnvalDA

#explanatory columns
    if length(independent_vars)==0
        independent_vars=string.(dfnames)
    end
    #delete main columns if they appear in the indep columns
    for remove_this_col in [numcol denomcol weightcol trnvalcol keycol][:]
        deleteat!(independent_vars, findin(independent_vars, [remove_this_col]))
    end
    n_indep=length(independent_vars)

    @assert issubset(Symbol.(treat_as_categorical_variable),dfnames)
    #n_char_features=0
    char_features=Vector{Symbol}()
    num_features=Vector{Symbol}()
    for nn = 1:n_indep
        this_is_a_string=false
        this_name=Symbol(independent_vars[nn])
        @assert in(Symbol(this_name),dfnames) "Error: explanatory column $(this_name) not found in data."
        if in(this_name,Symbol.(treat_as_categorical_variable))
            #change variable type in original dataframe
            if !is_categorical_column(dfin[this_name])
                info("DTM: Converting the column $(this_name) from $(eltype(dfin[this_name])) to string.")
                dfin[this_name]=string.(dfin[this_name])
                this_is_a_string=true
            else
                info("DTM: Variable $(this_name) is already categorical in the input data. No type conversion performed.")
            end
        end

        #if this_is_a_string||(eltype(dfin[this_name])<:AbstractString)
        #if (eltype(dfin[this_name])<:AbstractString)
        if is_categorical_column(dfin[this_name])
            push!(char_features,this_name)
        else
            push!(num_features,this_name)
        end

    end

 @assert (length(char_features)+length(num_features)==n_indep) #this should always be the case!
 n_char=length(char_features)
 n_num=length(num_features)

 #@show num_features
 #@show char_features
#warn("this is written in a terrible manner, we should improve this!!. It takes ages for larger data sets")
    if n_num>0
        #dfres=hcat(dfres,deepcopy(dfin[num_features]))
        dfres[num_features]=dfin[num_features]
    end
    if n_char>0
        #dfres=hcat(dfres,deepcopy(dfin[char_features]))
        dfres[char_features]=dfin[char_features]
    end
    #dftmp=DataFrame()

    #sort such that training observations are at the beginning of the file
    #this is not necessary!
    #sort!(dfres,cols = [:training_validation_bool, :irkey],rev=true)

    this_sett=ModelSettings()

    this_sett.df_name_vector=names(dfres)[1+global_const_shift_cols:end]
    updateSettingsMod!(this_sett,number_of_num_features=n_num,number_of_char_features=n_char) #model_type="build_tree",minw=.91,randomw=0,mf=.05,subsampling_features_prop=1.0,niter=10,subsampling_prop=1.00)
    info("Data initialization finished:")
    @show size(dfres)
    @show n_num,n_char
    @show char_features
    @show num_features
    return dfres,this_sett
end

function is_categorical_column(x)
    return eltype(x) <:AbstractString
end

function run_model_main(settingsFilename::String,dataFilename::String,datafolder::String,outfilename::String,outfileStringOnly::String;nrows::Int=-1)
	#this function distinguishes two cases: either a *.JLD file is loaded or a *.CSV file is converted to a *.jld file
	@assert length(settingsFilename)>3
	@assert length(dataFilename)>3
	const_shift_cols=global_const_shift_cols #this was 6 in earlier versions where the settings column was column 6, now the settings column is separate
	@assert isfile(settingsFilename) "Settings file $(settingsFilename) not found. Abort."
	if !(lowercase(dataFilename)[end-2:end]=="sql")&&!isfile(dataFilename);error("Data file $(dataFilename) not found.");end; # Did you provide a *.csv instead of a *.jld2 file (or vice-versa)?");end;
	#the next function is mainly here to avoid that there is a result from an old run, whis is interpreted as the result of the current run (which may have failed for instance)
		#we should check how much time it needs and maybe disable it eventually
		deleteAllOutputFiles(datafolder,outfileStringOnly)
	#Read Settings
		settingsArray,number_of_num_features=readSettings(settingsFilename)
	if size(settingsArray,1)>2
		#warn("Experimental: Settings File has more than two rows. Running multiple models")
		if lowercase(reverse(dataFilename)[1:4])!="dlj."
			#prep data and create *.jld file #NOTE: the fact that we save the *.jld file to disk is obsolete here (could be improved in the future) todo/tbd
			info("Prepping data for multirow run")
			settingsArrayTWOROWS=settingsArray[1:2,:]
			key,trn,numeratortrn,denominatortrn,weighttrn,trn_numfeatures,trn_charfeatures_PDA,keyval,val,val_numfeatures,val_charfeatures_PDA,numeratorval,denominatorval,weightval,mappings,sett,num_levels,char_levels,all_levels,all_levels_as_string_vector,names_and_levels,candMatWOMaxValues=prepare_non_jld_data(dataFilename,settingsArrayTWOROWS,datafolder,outfilename,outfileStringOnly,const_shift_cols,nrows,number_of_num_features)
			dataFilename,extp=splitext(dataFilename)
			dataFilename=string(dataFilename,".jld2")
			info("this will probably fail with a *.db or mysql  input file/reference")
			@assert isfile(dataFilename) "Something went wrong: dataFilename should be a prepped *.jld2 file here"
		end
		println("Loading Julia save file: \n $(dataFilename)")
		dictOfVariables=load(dataFilename);
		telapsed=toq()
		info("Time to read data: $(telapsed)")
		return run_model_multirow_settings(dataFilename,settingsArray,dictOfVariables,datafolder,outfilename,outfileStringOnly,const_shift_cols)
	end

	#check if JLD or CSV was provided
	if lowercase(reverse(dataFilename)[1:4])=="dlj." #dataFilename[end-4:end]
		println("Loading Julia save file: \n $(dataFilename)")
		dictOfVariables=load(dataFilename);
		telapsed=toq()
		info("Time to read data: $(telapsed)")
		return run_model_jld(dataFilename,settingsArray,dictOfVariables,datafolder,outfilename,outfileStringOnly,const_shift_cols)
	#non jld file was provided
	else
		key,trn,numeratortrn,denominatortrn,weighttrn,trn_numfeatures,trn_charfeatures_PDA,keyval,val,val_numfeatures,val_charfeatures_PDA,numeratorval,denominatorval,weightval,mappings,sett,num_levels,char_levels,all_levels,all_levels_as_string_vector,names_and_levels,candMatWOMaxValues=prepare_non_jld_data(dataFilename,settingsArray,datafolder,outfilename,outfileStringOnly,const_shift_cols,nrows,number_of_num_features)
		#key,trn,numeratortrn,denominatortrn,weighttrn,trn_numfeatures,trn_charfeatures_PDA,keyval,val,val_numfeatures,val_charfeatures_PDA,numeratorval,denominatorval,weightval,mappings,sett,num_levels,char_levels,all_levels,all_levels_as_string_vector,names_and_levels,candMatWOMaxValues=runmodel_non_jld(dataFilename,settingsArray,datafolder,outfilename,outfileStringOnly,const_shift_cols,nrows)
        telapsed=toq()
        sett.print_details&&(!sett.preppedJLDFileExists)&&info("Time to prepare data: $(telapsed)\n")
        return run_model_actual(key,trn,numeratortrn,denominatortrn,weighttrn,trn_numfeatures,trn_charfeatures_PDA,keyval,val,val_numfeatures,val_charfeatures_PDA,numeratorval,denominatorval,weightval,mappings,datafolder,outfilename,outfileStringOnly,const_shift_cols,sett,num_levels,char_levels,all_levels,all_levels_as_string_vector,names_and_levels,candMatWOMaxValues,dataFilename)
		#return runmodel_non_jld(dataFilename,settingsArray,datafolder,outfilename,outfileStringOnly,const_shift_cols,nrows)
	end
end

#function if dictionary from *.jld2 was loaded
function run_model_jld(dataFilename::String,settingsArray::Array{String,2},di::Dict,datafolder::String,outfilename::String,outfileStringOnly::String,const_shift_cols::Int)
	tic()
	oldsettings=di["oldsettings"]
	sett=ModelSettings() #default model settings
	updateSettings!(dataFilename,sett,settingsArray,copy(oldsettings.ncolsdfIndata),const_shift_cols,copy(oldsettings.df_name_vector))
    telapsed=toq()
    sett.print_details&&(!sett.preppedJLDFileExists)&&info("Time to prepare data: $(telapsed)\n")
    return run_model_actual(di["key"],di["trn"],di["numeratortrn"],di["denominatortrn"],di["weighttrn"],di["trn_numfeatures"],di["trn_charfeatures_PDA"],di["keyval"],di["val"],di["val_numfeatures"],di["val_charfeatures_PDA"],di["numeratorval"],di["denominatorval"],di["weightval"],	di["mappings"],datafolder,outfilename,outfileStringOnly,const_shift_cols,sett,di["num_levels"],di["char_levels"],di["all_levels"],di["all_levels_as_string_vector"],di["names_and_levels"],di["candMatWOMaxValues"],dataFilename)
end

function prep_data_from_df(df_userinput::DataFrame,sett::ModelSettings,fn_with_ext::String)
#warn: this function has a lot of redundancies with another function in this file ->prep_data_from_df
tic()

	datafilename,ext=splitext(fn_with_ext) #actually one can provide fn_with_ext="c:\\temp\\my folder\\out" (so no extension is necessary)
	datafolder,outfileStringOnly=splitdir(datafilename)
	if length(outfileStringOnly)==0
		@show fn_with_ext
		@show datafilename
		@show ext
		@show datafolder
		@show outfileStringOnly
		error("fn_with_ext should be of the form C:\\some\\folder\\my_output.txt (the extension is irrelevant) \n You provided: $(fn_with_ext) \nAbort.")
	end
	dfIndata=df_userinput
	println("Preparing Data for analysis...")

	@assert check_for_missing_data(dfIndata,global_const_shift_cols) "Error dataset contains missing (NA) values! ABORT."
   
	key=convert(Vector{String},string.(dfIndata[:irkey]))
	trn_val_idx=convert(Vector{UInt8},UInt8.(dfIndata[:training_validation_bool]))
	@assert sort(unique(trn_val_idx))==[0x00,0x01] "DTM: training validation column contains values not equal to 0 or 1." 
	#construct data
		numerator=convert(Array{Float64,1},dfIndata[:numerator])
		denominator=convert(Array{Float64,1},dfIndata[:denominator])
		weight=convert(Array{Float64,1},dfIndata[:weight])
		
	if any(denominator.<=0) 
	  warn("denominator[trn] contains observations with value 0.0")
	  foundobs=findfirst(denominator,0.0)
	  foundk=key[foundobs]
	  println("First observation: key=$(foundk), rownumber_trn=$(foundobs), value=$(denominator[foundobs])")
	  println("Denominator needs to be nonzero (negative values are allowed)")
	  println("Consider discarding or changing these values/rows.")
	  error("Abort.")
	end
	const strNegativeRatioWarning="Negative Ratios may lead to errors during the modelling process for models which rely on multiplicative residuals! \nIt is suggested that you adjust the data such that no negative values occur!"
	if any(x->x<0.0,denominator)
	  warn("There are training obsevations with negative values in the denominator column!")
	  warn(strNegativeRatioWarning)
	end
	#obsratiotrn=numeratortrn./denominatortrn
	if any(numerator.<0)
		warn("There are negative observed ratios!")
		warn(strNegativeRatioWarning)
    if typeof(sett.crit) == PoissonDevianceSplit
      error("DTM: Negative numerator values are not allowed for Poisson Deviance!. Abort. \n Crit = $(sett.critt)")
    end
	end

	df_names=sett.df_name_vector #names(dfIndata)
	#this DataFrame will contain both categorical and numerical features!
	features=DataFrame() 
	
	print("Preparing numeric variables...")
	#this function modifies features -> the num features are added as pdas!
	candMatWOMaxValues=add_coded_numdata!(dfIndata,sett,trn_val_idx,sett.max_splitting_points_num,features)
	println(" done.")
	#IMPORTANT convention  numerical data 'comes first' (i.e. first x columns are numerical), then categorical
	print("Preparing character variables...")
	for i=1:sett.number_of_char_features
		if !(eltype(dfIndata[global_const_shift_cols+i+sett.ishift+sett.number_of_num_features])<:AbstractString) 
			@show eltype(dfIndata[global_const_shift_cols+i+sett.ishift+sett.number_of_num_features])
			@show dfIndata[1:min(10,size(dfIndata,1)),global_const_shift_cols+i+sett.ishift+sett.number_of_num_features]
			error("Character variable $(i) = $(sett.df_name_vector[i+sett.number_of_num_features]) is not <:AbstractString in the dataframe.")
			#warn("It may be that the variable was exported as character (by SAS) even though it is infact an integer. Please check!")			
		end		
	  try
		colnumber=global_const_shift_cols+sett.ishift+sett.number_of_num_features+i
		thisname=Symbol(df_names[i+sett.number_of_num_features]) #df_names[colnumber]
		thiscol_as_utf8=dfIndata[thisname] #view(dfIndata,:,colnumber)
		if eltype(thiscol_as_utf8)!=String
			thiscol_as_utf8=convert(Array{String,1},thiscol_as_utf8)
		end
		features[thisname]=PooledArray(thiscol_as_utf8)
		if length(features[thisname].pool)>255
			@show nlevels=length(features[thisname].pool)
			warn("DTM: Variable $(string(thisname)) has more than 255 levels, namely $(nlevels)\nThe tree may take very long to be constructed.")			
		end
	  catch this_error
		@show this_error
		error("DTM: Unable to initialize character variable $(i) = $(sett.df_name_vector[i+sett.number_of_num_features]) (see error message above).")
	  end	  
	end
	mappings=deepcopy(map(i->features[i].pool,sett.number_of_num_features+1:size(features,2)))
	
	println(" done.")
	wlo,whi=extrema(weight)
	if wlo!=whi
		print("Max weight is $(whi), min weight is $(wlo)")
		warn("Currently the automatic choice of splitting points does not consider the weight distribution!")
	end
	
	pools=map(i->features[i].pool,1:size(features,2)) 
	pool_lengths=length.(pools)
	for i=1:length(pool_lengths)
		if pool_lengths[i]==1
			info("DTM: Column $(i) = $(sett.df_name_vector[i]) has zero splitting points.")  
		end
	end
		
	trnidx,validx=createIntegerIndices(trn_val_idx)
	dtmtable=DTMTable(key,trnidx,validx,numerator,denominator,weight,features,candMatWOMaxValues,mappings)
	
	if sett.boolSaveJLDFile
        jldfile=string(fileroot(datafilename),".jld2")
        println("Saving prepped data to file:\n $(jldfile)")
		isfile(jldfile)&&rm(jldfile)
		#@time save(jldfile,"dtmtable",dtmtable,"candMatWOMaxValues",candMatWOMaxValues,"features",features,"mappings",mappings,"key",key,"trn_val_idx",trn_val_idx,"numerator",numerator,"denominator",denominator,"weight",weight,"oldsettings",sett)
		@time save(jldfile,"dtmtable",dtmtable,"candMatWOMaxValues",candMatWOMaxValues,"mappings",mappings,"sett",sett)
		println("JLD2 file saved. Starting modelling... \n")
	end

	return dtmtable
end


function prepare_non_jld_data(dataFilename::String,settingsArray::Array{String,2},datafolder::String,outfilename::String,outfileStringOnly::String,const_shift_cols::Int,nrows::Int,number_of_num_features::Int)
local eltypev,dfIndata
	if lowercase(dataFilename)[end-2:end]==".db"
		info("todo: Test *.DB Sqlite functionality.")
		info("the current version may work, but it will likely fail if there are any non standard characters in the data (i.e. non \'1-9,a-Z\')")
		warn("This will not work properly if there are any non-standard characters (such as umlauts)!") #the sqlite driver has issues with ? and so on
		println("Importing data from SQLiteDB file: \n $(dataFilename)")
		try
			info("check if disabling gc speeds this up")
			dfIndata=readDFfromSQLiteDB(dataFilename,const_shift_cols)
		catch sqllitereaderr
			@show sqllitereaderr
			error("Unable to read from SQLiteDB. See above.")
		end
	elseif lowercase(dataFilename)[end-2:end]=="sql"
		print("Importing data from SQL db.")
		warn("If this uses to  much memory implement limit 100,300 !!")
		warn("like this we can read stepwise and append to the df in julia!")
		try
			dfIndata=readDFfromSQLDB()
		catch eri
			@show eri
			error("Unable to read from SQL DB. See above.")
		end
		println(" done.")
	else
		#CSV was provided, Reading Data
			println("Importing CSV data: \n $(dataFilename)")
        #Read only a few rows to determine AbstractString/Numeric columns
		  try
			eltypev=define_eltypevector(DataFrames.readtable(dataFilename,nrows=100),const_shift_cols,number_of_num_features)
		  catch eri
			warn("Readtable failed! Check the size of the *.CSV: $(dataFilename)") #sometimes the file only consists of the header (e.g. if the SAS export failed), then readtable will fail.
			fz=0
			fz=filesize(dataFilename)/1024
			warn("Filesize is $(fz) KB")
			@show eri
			@assert false
		  end
		#Read whole dataset
		try
			dfIndata=readtable(dataFilename,eltypes=eltypev,nrows=nrows) #this currently uses by far to much memory on a 14m row file. The limitation is about a 7m row file with 77 columns
		catch readerror
			@show readerror
			error("Unable to read file $(dataFilename)")
		end
	end

	#at this point the data exists as a dataframe
		telapsed=toq()
		println("Time to read data: $(telapsed)")

	key,trn,numeratortrn,denominatortrn,weighttrn,trn_numfeatures,trn_charfeatures_PDA,keyval,val,val_numfeatures,val_charfeatures_PDA,numeratorval,denominatorval,weightval,mappings,sett,num_levels,char_levels,all_levels,all_levels_as_string_vector,names_and_levels,candMatWOMaxValues=prep_data_from_df(dataFilename,settingsArray,dfIndata,datafolder,outfilename,outfileStringOnly,const_shift_cols)
return 	key,trn,numeratortrn,denominatortrn,weighttrn,trn_numfeatures,trn_charfeatures_PDA,keyval,val,val_numfeatures,val_charfeatures_PDA,numeratorval,denominatorval,weightval,mappings,sett,num_levels,char_levels,all_levels,all_levels_as_string_vector,names_and_levels,candMatWOMaxValues
end

function run_model_actual(dtmtable::DTMTable,sett::ModelSettings,fn::String)
#this function is called with prepped julia data, this is the core modelling function (all previous ones are for preparational tasks only)
tic() #total modelling time , prnt&&println("Modelling finished. Time: $(now()) - Total time was $(round(elapsed_until_this_point,1))s = $(round(elapsed_until_this_point/60,1))m")

path_and_fn_wo_extension,ext=splitext(fn)
if in(lowercase(ext),[".txt",".sas",".jl",".pdf"])
else
	@show fn
	@show path_and_fn_wo_extension,ext
	error("DTM: Unexpected file extension $(ext)")	
end
datafolder,outfileStringOnly=splitdir(path_and_fn_wo_extension)

filelistWithFilesToBeZipped=String[]
csharp_code_file=string(path_and_fn_wo_extension,".cs")
vba_code_file=string(datafolder,"\\",outfileStringOnly,".vba.txt")
sas_code_file=string(path_and_fn_wo_extension,".sas")
tree_file=string(path_and_fn_wo_extension,".txt")
dot_graph_TXT_file=string(path_and_fn_wo_extension,".dot.txt")
dot_graph_visualization_file=string(path_and_fn_wo_extension,".dot.pdf")
this_outfile=convert(String,string(path_and_fn_wo_extension,".csv"))
statsfile=string(path_and_fn_wo_extension,".stats.csv")
jldresultsfile=string(path_and_fn_wo_extension,".jld2")
statsfileExcel=string(path_and_fn_wo_extension,".xlsx")

#for now trnindx and validx should always be sorted!
if !issorted(dtmtable.trnidx) 
	info("DTM: trnidx was not sorted. Sorting trnindex")
	sort!(dtmtable.trnidx)
end
if !issorted(dtmtable.validx) 
	info("DTM: validx was not sorted. Sorting validx")
	sort!(dtmtable.validx)
end

key=dtmtable.key
trnidx=dtmtable.trnidx
validx=dtmtable.validx
features=dtmtable.features
weight=dtmtable.weight
numerator=dtmtable.numerator
denominator=dtmtable.denominator

#define 'version', we store the sha1 hash of the current commit in sett.version
sett.version=get_sha1()

#initialize random number generator
	#we may want to make this independent from the data (trn), then again different data should/will lead to a different result.
#hash of DF is not working properly, thus we need to iteration over the columns of the DF (this can be simpplified later on.)
s00=hash(trnidx,hash(validx,hash(sett.seed,hash(numerator,hash(denominator,hash(weight))))))
for x=1:size(features,2)
	s00=hash(features[x],s00)
end
srandInt=floor(Int,1/3*hash(931,s00))
srand(srandInt)
	#srand(floor(Int,1/3*hash(931,hash(features,hash(trnidx,hash(validx,hash(sett.seed,hash(numerator,hash(denominator,hash(weight))))))))))
	#@show rand();@show rand();@show rand(); #check reproducability
if sett.subsampling_prop <1.0 
	warn("BK: subsampling_prop <1. This can be terribly slow in the current implementation (need to improve this)!")
end

#there are at least two ways to derive estimates for hold out data: either loop through the leaves of the tree or loop over the rows of the validation data setting
#also there are two types of hold out estimates: either the ones derived from training data (more meaningful case) or the ones derived from validation labels
sumwtrn=sum(view(weight,trnidx))
minweighttmp = sett.minw < 0.0 ? sumwtrn*(-max(-1.0,sett.minw)) : sett.minw
@assert minweighttmp>=0
sumwtrn>=2*minweighttmp||warn("DTM: Specified minweight is larger than half of the total training weight!! No split is possible!")

prnt=sett.print_details
general_settings=convert(String,string("Write tree to txt file: $(sett.bool_write_tree), statsByVariables=$(join(sett.statsByVariables,','))"))
#general_settings=convert(String,string(general_settings,))
	char_levels=map(x->length(features[x].pool),sett.number_of_num_features+1:sett.number_of_num_features+sett.number_of_char_features)
	num_levels=map(x->length(features[x].pool),1:sett.number_of_num_features)
		if size(sett.moderationvector,1)>1;if (sett.adaptiveLearningRate!=1.0 && sett.adaptiveLearningRate!=0.0);println("######################################## \nWarning: (bk) Adaptive Learning has been disabled as a moderation vector has been provided. \n######################################## ");end;sett.adaptiveLearningRate=1.0;end;
      if prnt
		  println(stars)
		  println("---General Settings------------------------------------------------------------")
		  println(general_settings)
		  print_some_settings(sett,["write_sas_code","write_iteration_matrix","write_result","write_statistics","boolCreateZipFile","write_csharp_code","statsRandomByVariable","boolSaveJLDFile"])
		  println("---Data Summary----------------------------------------------------------------")
		  println("Number of records - Trn: $(length(trnidx)) Val: $(length(validx))")
		  println("Number of numeric/character variables: $(sett.number_of_num_features)/$(sett.number_of_char_features)")
      end
	  if size(num_levels,1) >0
        prnt&&println(string("Number of levels for numeric variables: (Maximal # of linear steps: $(sett.max_splitting_points_num))"))
        num_level_string=join(sort(num_levels, rev=true),",")
        prnt&&println(string(" ",num_level_string))
      else
        prnt&&println("The data does not contain any numeric attributes.")
        num_level_string=""
	  end
	  if size(char_levels,1) >0
        prnt&&println(string("Number of levels for character variables: "))
        #println(arr_to_str(sort(char_levels, rev=true)))
        char_level_string=join(sort(char_levels, rev=true),",")
        prnt&&println(string(" ",char_level_string))
      else
        char_level_string=""
        prnt&&println("The data does not contain any character attributes.")
      end
	  prnt&&println("catSortByThreshold: $(sett.catSortByThreshold), catSortBy: $(sett.catSortBy), smoothEstimates: $(sett.smoothEstimates)")

  #create a string which contains all model settings (this will be written to output file(s))
      model_setting_string=convert(String,string("Number of records - Trn: $(length(trnidx)) Val: $(sum(validx)) \nNumber of levels for numeric variables: (Maximal # of linear steps: $(sett.max_splitting_points_num)) "))
      model_setting_string=convert(String,string(model_setting_string, "\n",num_level_string, "\nNumber of levels for character variables: \n", char_level_string ))
      model_setting_string=convert(String,string(model_setting_string, "\n", general_settings))
      model_setting_string=convert(String,string(model_setting_string, "\n", "Type: $(sett.model_type) \nMinweight: $(sett.minw) \nRnd: $(sett.randomw) \nNumMaxSplitPoints: $(sett.max_splitting_points_num) \nCriterion: $(sett.crit) \nIterations: $(sett.niter) \nModeration Vector: $(join(sett.moderationvector,",")) \nNScores: $(sett.nscores) \nBoolStartAtMean: $(sett.BoolStartAtMean) "))
      model_setting_string=convert(String,string(model_setting_string, "\n", "Datafolder= \n",datafolder,"\nOutfilename= \n",outfileStringOnly," "))
      model_setting_string=convert(String,string(model_setting_string, "\n\n", "CSV column names (in the order in which Julia processes them):\n",join(sett.df_name_vector,",")))
      #permutedims is "transpose" in julia 0.5
	  #model_setting_string=convert(String,string(model_setting_string, "\n\n", "Names:Levels: \n",join(permutedims(names_and_levels, (2, 1)),":")))
      model_setting_string=convert(String,string(model_setting_string, "\n", "Folder of Julia code (decisiontree.jl): \n", dirname(@__FILE__)))
	  model_setting_string=convert(String,string(model_setting_string, "\n", "Julia started at: \n$(juliaStartedat) "))
      model_setting_string=convert(String,string(model_setting_string, "\n", "Julia (Modelling) start time: \n$(now()) "))

prnt&&println("---Model Settings--------------------------------------------------------------")
      if sett.boolRankOptimization && "build_tree"==sett.model_type
		  error("Error: (bk) I do not think that _RankOpt is currently working with a single tree....")
	  end

	 local resulting_model,resultEnsemble

	  if "build_tree"==sett.model_type
		  if prnt
			  println("Type: $(sett.model_type), Minweight: $(sett.minw), Rnd: $(sett.randomw)")
			  println("NumMaxSplitPoints: $(sett.max_splitting_points_num), Criterion: $(sett.crit)")
			  println(stars)
			  println(string(""))
			  println(string("Building Tree... Time: $(now())"))
		  end
		  obsvalcols=1
		  if sett.randomw>0
			warn("Running single tree with positive random weight:$(sett.randomw)")
		  end

	#build a simple tree
		fitted_values_tree=zeros(numerator)		
		@time tree=build_tree!(trnidx,validx,dtmtable.candMatWOMaxValues,dtmtable.mappings,sett,numerator,denominator,weight,features,fitted_values_tree)
		resulting_model=tree
		println("Deriving Trn Estimates... Time: $(now())")
		leaves_of_tree=create_leaves_array(tree)
        leaf_number_vector=zeros(Int,length(fitted_values_tree))
        #todo/tbd possibly improve this by only applying the function once!
        #notably fitted_values_tree is overwritten here, however the values do NOT change.
        apply_tree_by_leaf!(fitted_values_tree,leaf_number_vector,trnidx,tree.rootnode,features)
        if length(validx)>0
            apply_tree_by_leaf!(fitted_values_tree,leaf_number_vector,validx,tree.rootnode,features)
		end		
		 prnt&&println("Preparing output... Time: $(now())")
         #1-dim Variable Importance
          var_imp1d_str_arr,var_imp2d_str_arr,onedimintvec_unused,twodimintvec_unused=variable_importance(leaves_of_tree,sett.df_name_vector,sett.number_of_num_features)
          nleaves=Int64(size(leaves_of_tree,1))
		 res=hcat(key,fitted_values_tree)
		#Create output statistics		
			#this is not done in an efficient manner yet..... we are copying data to get est[trnidx] etc...
			@time xlData,overallstats=createSingleTreeExcel(trnidx,validx,sett,tree,leaves_of_tree,fitted_values_tree,leaf_number_vector,numerator,denominator,weight)
			resulting_model.exceldata=xlData
			if sett.write_statistics
				writeStatisticsFile!(statsfileExcel,xlData,filelistWithFilesToBeZipped)
			end
			#z=transpose(overallstats) #this should work, but it somehow doesnt.... https://discourse.julialang.org/t/transpose-not-working-on-array-any-2/6510
			z=permutedims(overallstats,(2,1))
			dfStats=DataFrame(convert(Array{Float64,2},z[2:size(z,1),:]))
			names!(dfStats,Symbol[Symbol(x) for x in view(z,1,:)])
			resulting_model.modelstats=dfStats
			if sett.boolSaveResultAsJLDFile
				println("Saving results to jld2 file: \n $(jldresultsfile)")
				isfile(jldresultsfile)&&rm(jldresultsfile)
				@time @save jldresultsfile xlData res leaves_of_tree trnidx validx fitted_values_tree tree var_imp1d_str_arr var_imp2d_str_arr onedimintvec_unused twodimintvec_unused
				push!(filelistWithFilesToBeZipped,jldresultsfile)
			end

		 #print_tree(tree)
		 dot_graph=graph(resulting_model)
         model_setting_string=convert(String,string(model_setting_string, "\n", "Julia end time: \n$(now()) \n"))
		 if sett.write_dot_graph
			println("Writing DOT tree structure to file: \n $(dot_graph_TXT_file)")
			isfile(dot_graph_TXT_file)&&rm(dot_graph_TXT_file)
			f=open(dot_graph_TXT_file,"w");
				write(f,"\n\n",dot_graph,"\n\n");
			close(f);
			draw_dot_graph(sett.graphvizexecutable,dot_graph_TXT_file,dot_graph_visualization_file)
			push!(filelistWithFilesToBeZipped,dot_graph_TXT_file)
			isfile(dot_graph_visualization_file)&&push!(filelistWithFilesToBeZipped,dot_graph_visualization_file)
		end
		if sett.bool_write_tree			
			#write settings
			println("Writing tree structure to file: \n $(tree_file)")				
				isfile(tree_file)&&rm(tree_file)
				f=open(tree_file,"w")
					write(f,model_setting_string)
					write(f,"\n\n",dot_graph,"\n\n")
				close(f);
				@time write_tree(dtmtable.candMatWOMaxValues,tree,sett.number_of_num_features,var_imp1d_str_arr,var_imp2d_str_arr,0,tree_file,sett.df_name_vector,dtmtable.mappings)
				push!(filelistWithFilesToBeZipped,tree_file)			
            end
			if sett.write_sas_code
				#write SAS Code
				println("Writing SAS Code: \n $(sas_code_file)")
				isfile(sas_code_file)&&rm(sas_code_file)
				@time write_sas_code(dtmtable.candMatWOMaxValues,tree,sett.number_of_num_features,sas_code_file,sett.df_name_vector,model_setting_string,dtmtable.mappings,1.0)
				push!(filelistWithFilesToBeZipped,sas_code_file)
			end
		#Write result to CSV
		#writecsv(this_outfile,res)
		#prrr=sortperm(float.(res)[:,1])
		#res[:,2]=round.(res[:,2],13)
		#writecsv(string(this_outfile,"sorted"),res[prrr,:])
			if sett.write_result
				println("Exporting Data: \n $(this_outfile)")
				isfile(this_outfile)&&rm(this_outfile)
				@time writecsv(this_outfile,res)
				push!(filelistWithFilesToBeZipped,this_outfile)
            end
	#end build tree
	elseif in(sett.model_type,["bagged_tree","boosted_tree"])
		if prnt
			println("Type: $(sett.model_type),Iterations: $(sett.niter),Minweight: $(sett.minw)")
			if size(sett.moderationvector,1)>1;println("Moderationvector = $(join(sett.moderationvector,",")) ");else println("Moderation Factor: $(sett.mf)"); end;
			println("Rnd: $(sett.randomw), NumMaxSplitPoints: $(sett.max_splitting_points_num), NScores: $(sett.nscores), Criterion: $(sett.crit)")
			println("boolRandomizeOnlySplitAtTopNode: $(sett.boolRandomizeOnlySplitAtTopNode), sett.boolRankOptimization: $(sett.boolRankOptimization)")
			println("subsampling_features_prop: $(sett.subsampling_features_prop), subsampling_prop: $(sett.subsampling_prop)")
			println(stars)
		end

		if "boosted_tree"==sett.model_type
			@timeConditional(sett.print_details,begin
				xlData,estimatedRatio,vectorOfLeafNumbers,vectorOfLeafArrays,rawObservedRatioPerScore,est_matrixFromScores,stats,estimateUnsmoothed,estimateSmoothed,estimateFromRelativities,resultEnsemble=boosted_tree(dtmtable,sett)
			end)
			
			model_setting_string=convert(String,string(model_setting_string, "\n", "Julia end time: \n$(now()) \n "))
			#this line is for  the sas,vba and csharp code
			estimatesPerScoreForCode = sett.smoothEstimates=="0" ? rawObservedRatioPerScore : resultEnsemble.ScoreToSmoothedEstimate
			
			if sett.write_sas_code
				#write SAS Code
				println("Writing SAS Code: \n $(sas_code_file)")
				isfile(sas_code_file)&&rm(sas_code_file)
				if sett.roptForcedPremIncr
					error("this does currently not work")
				else
					@time write_sas_code(estimatesPerScoreForCode,dtmtable.candMatWOMaxValues,resultEnsemble,sett.number_of_num_features,sas_code_file,sett.df_name_vector,model_setting_string,dtmtable.mappings)
				end
				push!(filelistWithFilesToBeZipped,sas_code_file)
			end
			if sett.write_csharp_code				
				println("Writing C# Code: \n $(csharp_code_file)")
				isfile(csharp_code_file)&&rm(csharp_code_file)				
				@time write_csharp_code(vectorOfLeafArrays,estimatesPerScoreForCode,dtmtable.candMatWOMaxValues,resultEnsemble,csharp_code_file,model_setting_string,dtmtable.mappings,0,sett)
				push!(filelistWithFilesToBeZipped,csharp_code_file)
			end 
			if sett.write_vba_code
				println("Writing VBA Code: \n $(vba_code_file)")
				isfile(vba_code_file)&&rm(vba_code_file)
				@time write_vba_code(vectorOfLeafArrays,estimatesPerScoreForCode,dtmtable.candMatWOMaxValues,resultEnsemble,vba_code_file,model_setting_string,dtmtable.mappings,0,sett)
				push!(filelistWithFilesToBeZipped,vba_code_file)
			end
		#end #boosting
		elseif "bagged_tree"==sett.model_type
			#bagging
			@timeConditional(sett.print_details,begin
			warn("bagging is not yet updated for dtm2! This will most likely fail!")
				xlData,estimatedRatio,vectorOfLeafNumbers,vectorOfLeafArrays,rawObservedRatioPerScore,est_matrixFromScores,stats,estimateUnsmoothed,estimateSmoothed,estimateFromRelativities,resultEnsemble=bagged_tree(dtmtable,sett)
				#xlData,vectorOfLeafNumbersTrn,vectorOfLeafArrays,rawObservedRatioPerScore,est_matrixFromScores,est_matrixFromScoresVAL,stats,scoresval,est_matrix_val,estimateUnsmoothedVal,estimateSmoothedVal,estimateFromRelativitiesVal,estimateUnsmoothedTrn,estimateSmoothedTrn,estimateFromRelativitiesTrn,resultEnsemble=bagged_tree(mappings,candMatWOMaxValues,sett,numeratortrn,denominatortrn,weighttrn,trn_numfeatures,trn_charfeatures_PDA,numeratorval,denominatorval,weightval,val_numfeatures,val_charfeatures_PDA)
			end)
			#
			model_setting_string=convert(String,string(model_setting_string, "\n", "Julia end time: \n$(now()) \n "))
		else
			error("Unknown Model Type: $(sett.model_type) (Within Bagging/Boosting loop)")
		end #end bagging or boosting

resulting_model=resultEnsemble

if sett.write_result
	res=hcat(key,resultEnsemble.scores,estimateUnsmoothed,estimateSmoothed,estimatedRatio)
end	
#save results as *.JLD2 file
if sett.boolSaveResultAsJLDFile
	println("Saving results to jld2 file: \n $(jldresultsfile)")
	isfile(jldresultsfile)&&rm(jldresultsfile)
	@time @save jldresultsfile xlData dtmtable trnidx validx vectorOfLeafNumbers vectorOfLeafArrays rawObservedRatioPerScore est_matrixFromScores stats estimateUnsmoothed estimateSmoothed estimateFromRelativities resultEnsemble
	push!(filelistWithFilesToBeZipped,jldresultsfile)
end

    #after boosting bagging  export etc
#here we are in the "section": 	#if in(sett.model_type,["bagged_tree","boosted_tree"])
#output for boosting & bagging#write tree
	if sett.bool_write_tree
				#write settings
				println("Writing tree structure to file: \n $(tree_file)")
				isfile(tree_file)&&rm(tree_file)
				f=open(tree_file,"w");write(f,model_setting_string);close(f);
				@time write_tree(dtmtable.candMatWOMaxValues,resultEnsemble,sett.number_of_num_features,0,tree_file,sett.df_name_vector,dtmtable.mappings)
				push!(filelistWithFilesToBeZipped,tree_file)
	end
	#Write statistics
	if sett.write_statistics
			writeStatisticsFile!(statsfileExcel,xlData,filelistWithFilesToBeZipped)
		 end
	#Write result to CSV
		if sett.write_result
			println("Exporting Data: \n $(this_outfile)")
			isfile(this_outfile)&&rm(this_outfile)
			@time writecsv(this_outfile,res)
			push!(filelistWithFilesToBeZipped,this_outfile)
		end
	#write estimates matrix
		if sett.write_iteration_matrix
			#Matrix with Leaf Numbers
			#currently we only have the training mat
			leafnrfile=string(path_and_fn_wo_extension,".leafnumbers.csv")
			#the first column is 0 until here
			#we now fill it with the irkeys
			vectorOfLeafNumbersTrnMOD=hcat(key,vectorOfLeafNumbersTrn[:,2:end])
			println("Exporting LeafNumber Matrix: \n $(leafnrfile)")
				isfile(leafnrfile)&&rm(leafnrfile)
				@time writecsv(leafnrfile,vectorOfLeafNumbersTrnMOD)
				push!(filelistWithFilesToBeZipped,leafnrfile)
			#MATRIX 1
				estmat_outile2=string(path_and_fn_wo_extension,"_iteration_matrixFromScores.csv")
				res_estmatrixFromScores=hcat(vcat(key,keyval),vcat(est_matrixFromScores,est_matrixFromScoresVAL))
				println("Exporting Estimates Matrix (based on Scores): \n $(estmat_outile2)")
				isfile(estmat_outile2)&&rm(estmat_outile2)
				@time writecsv(estmat_outile2,res_estmatrixFromScores)
				push!(filelistWithFilesToBeZipped,estmat_outile2)
			#MATRIX 2
				estmat_outile=string(path_and_fn_wo_extension,"_iteration_matrix.csv")
				res_estmatrix=hcat(vcat(key,keyval),vcat(resultEnsemble.iterationmatrix,est_matrix_val))
				println("Exporting Estimates Matrix: \n $(estmat_outile)")
					#println("Saving iteration matrix to jld file; this is only a backup in case the CSV export does not work (which sometimes happens with very large files)")
					#@time save(string(estmat_outile,".jld"),"res_estmatrix",res_estmatrix)
				isfile(estmat_outile)&&rm(estmat_outile)
				@time writecsv(estmat_outile,res_estmatrix)
				push!(filelistWithFilesToBeZipped,estmat_outile)
		end
else
	error("Unknown Model Type: $(sett.model_type)")
end #end distinction between three model types

#NOTE: the option boolCreateZipFile should be disabled when code is run through grails
#as the *.7z file will be created by the listener (because the *.log file will only exist when run_model terminated

#this code applis to all model types
   #Create ZIP file
     #Push Log file to file list
	 #println("\n")
	 elapsed_until_this_point=toq()
	 prnt&&println("Modelling finished. Time: $(now()) - Total time was $(round(elapsed_until_this_point,1))s = $(round(elapsed_until_this_point/60,1))m")
	 if sett.boolCreateZipFile
		 #strBasename,ext=splitext(basename(dataFilename))
		 #logfile=string(datafolder,"\\",strBasename,".log")
		 logfile=string(path_and_fn_wo_extension,".log")
		 tmplogfile="c:\\temp\\julia_logfile.log";isfile(tmplogfile)&&rm(tmplogfile)
		 try
			isfile(logfile)&&isreadable(logfile)&&cp(logfile,tmplogfile)
		 catch error_while_copying_log
			@show error_while_copying_log
		 end
		 isfile(tmplogfile)&&push!(filelistWithFilesToBeZipped,tmplogfile)
		 zipFilename=string(path_and_fn_wo_extension,".7z")
		 println("Creating zip file with results \n $(zipFilename)")
		 @time createZipFile(zipFilename,filelistWithFilesToBeZipped)
	end
  #result needs to have the same type, irrespective if we ran a tree or a boosted tree
  #return true

 return filelistWithFilesToBeZipped,resulting_model
end

"""
If run_model_actual is called with an array of settings, the model is run on all settings.
The model type (e.g. boosting/tree/bagging) needs to be the same for each entry of settings
"""
function run_model_actual(dtmtable::DTMTable,setts::Vector{ModelSettings},fn::String)
	if length(unique(map(x->x.model_type,setts)))!=1 
		uq_model_types = unique(map(x->x.model_type,setts))
		error("DTM: model_type must be the same for all elements of the settings vector. You provided these model_types $(uq_model_types)")
	end

	path_and_fn_wo_extension="some"
	header=Vector{String}(1)
	header_settings=Vector{String}(1)
	allstats=Array{Float64,2}
	allsettings=Array{Any,2}
	i=1
	nsetts=length(setts)
	for sett in setts 
		path_and_fn_wo_extension,ext=splitext(fn)
		path_and_fn_wo_extension_mod=string(path_and_fn_wo_extension,"_",i)
		fnmod=string(path_and_fn_wo_extension_mod,ext)
		somestrings,model=run_model_actual(dtmtable,sett,fnmod)
		desc,numbrs,desc_settingsvec,settingsvec=get_stats(model)
		if i==1
			header=deepcopy(desc)
			header_settings=deepcopy(desc_settingsvec)
			allstats=Array{Float64,2}(length(numbrs),nsetts)			
			allsettings=Array{Any,2}(length(settingsvec),nsetts)			
		else			
			if !all(header.==desc)
				warn("Model run $(i) returned an unexpected results vector:")
				@show desc
				@show header
			end
			if !all(header_settings.==desc_settingsvec)
				warn("Model run $(i) returned an unexpected results vector:")
				@show desc_settingsvec
				@show header_settings
			end
		end
		allstats[:,i].=deepcopy(numbrs)
		allsettings[:,i].=deepcopy(settingsvec)
		i+=1
	end	

	statsdf=DataFrame(transpose(allstats))
	names!(statsdf,Symbol.(header))
	settsdf=DataFrame(permutedims(allsettings,[2,1]))
	names!(settsdf,Symbol.(header_settings))

	@show fld,namestr=splitdir(path_and_fn_wo_extension)
	filen=string(fld,"multistats.xlsx")
	if isfile(filen)
		info("Deleting $(filen).")
		rm(filen)
	end

	return statsdf,settsdf
end

function run_legacy(mld::String,fld::String)
    @assert isdir(fld)
    sett_file=string(fld,"\\",mld,".settings.csv")
    @assert isfile(sett_file)
    csv_file=string(fld,"\\",mld,".csv")
    @assert isfile(csv_file)
    args=[sett_file csv_file string(mld,"_out_")]
    rs=run_model(args)
    return rs
end
