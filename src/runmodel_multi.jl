
function run_model_multirow_settings(dataFilename::String,multirowSettingsArray::Array{String,2},di::Dict,datafolder::String,outfilename::String,outfileStringOnlyOrig::String,const_shift_cols::Int)
	warn("Experimental: Multirow settings file. Currently this is only supported for Boosting and normal Trees!")
	oldsettings=di["oldsettings"]
	sett=ModelSettings() #default model settings	
	nmodels=size(multirowSettingsArray,1)-1
	@assert nmodels>1
	#this is a hack to suppress a warning in updateSettings! which depends on the ending of dataFilename (here we won't read any data anyway)	
		dataFilenameMOD=string(dataFilename[1:end-2],"jld") 	
	#Parallelized Run
	#NOTE: it may make sense to submit a parameter to calcIterationsPerCore2 which is larger than "nrpocs()" 
	#the pmap function will distribute the work evenly among all processes anyway
	#this approach (with pmap) should be much more efficient if the work is not even for all iterations!
				
		#send data to all workers: without this command we will have a lot of overhead to send the data to the processes; but the data is constant and only needs to be sent ONCE!)
			sendto_module(DTM,workers(),local_data_dict=deepcopy(di))
		
		warn("defaulted_modelstats_df! this may not work!")
		
		df_modelstats_loc="R:\\TariffWatch New Repo\\algorithms\\Julia\\Code\\defaulted_modelstats_df.jld"
		loaded_dict = JLD.load(df_modelstats_loc);
		global const defaulted_modelstats_df = loaded_dict["defaulted_modelstats_df"]
		@show defaulted_modelstats_df
		
		#fist iteration (with header)
		v0=@spawn run_model_multirow_perCore(defaulted_modelstats_df,[1 1],deepcopy(sett),deepcopy(oldsettings),dataFilenameMOD,multirowSettingsArray,local_data_dict,datafolder,outfilename,outfileStringOnlyOrig,const_shift_cols)		
		#check: we need to avoid that the data is transferred for every call here! (not sure how to check this) Initial timings suggest that we have no issue
		pmapresult=pmap(ii->
			run_model_multirow_perCore(defaulted_modelstats_df,[1 ii],deepcopy(sett),deepcopy(oldsettings),dataFilenameMOD,multirowSettingsArray,local_data_dict,datafolder,outfilename,outfileStringOnlyOrig,const_shift_cols)
			,2:nmodels)
		v0fetched=fetch(v0)		
		info("Aggregating results...")
		allmodelStats=reduce(vcat, v0fetched, pmapresult)
		info("Writing output...")
	#old approach:
		#models_per_core=calcIterationsPerCore2(size(multirowSettingsArray,1)-1,max(1,nprocs()-1)) #I think we need to submit nprocs()-1 here, or does the main task (=worker1?) also do any work? todo/tbd need to test this.			
		#allmodelStats = @parallel (vcat) for ii=1:sum(models_per_core[:,1].>0) #only run on the cores where there is actual work to do				
		#	run_model_multirow_perCore([1 ii],deepcopy(sett),deepcopy(oldsettings),
		#			dataFilenameMOD,multirowSettingsArray,di,datafolder,outfilename,outfileStringOnlyOrig,const_shift_cols)
		#end
	#note, we may want to reorder allmodelStats.... but maybe it is not necessary
	#write to Excel file 				
		outf=string(datafolder,"\\",outfileStringOnlyOrig,".allstats.xlsx")
		multirowSettingsArray2=hcat(["Model Number";collect(1:size(multirowSettingsArray,1)-1)],multirowSettingsArray)
		info("Writing overall model stats file \n $(outf)") 
			if size(allmodelStats,1)<1048576 #maxrows of Excel 1<<20
				write_any_array_to_excel_file(allmodelStats,multirowSettingsArray2,outf)
			else
				writecsv(string(datafolder,"\\",outfileStringOnlyOrig,".allstats.csv"),allmodelStats)
				writecsv(string(datafolder,"\\",outfileStringOnlyOrig,".allmodelsettings.csv"),allmodelStats)			
			end
		info("Done. Time - $(now())")
	return String["true"],"No Model was returned as this was a multirow run" #need to think about what we return here....
end
	
function run_model_multirow_perCore(defaulted_modelstats_df::DataFrame,iteridx::Array{Int,2},sett::ModelSettings,oldsettings::ModelSettings,dataFilename::String,multirowSettingsArray::Array{String,2},di::Dict,datafolder::String,outfilename::String,outfileStringOnlyOrig::String,const_shift_cols::Int)	
	@assert size(iteridx)==(1,2)
	local allmodelStats
	local thesem
	allmodelStats=deepcopy(defaulted_modelstats_df)
	try
		for i=iteridx[2]:iteridx[2]+iteridx[1]-1
			try
				tic()
				tic()
				info("Process: $(myid()) - Running setting entry $(i)")
				outfileStringOnly=string(outfileStringOnlyOrig,i)
				outfilename=convert(String,string(datafolder,"\\",outfileStringOnly))
				settingsArray=deepcopy(multirowSettingsArray[[1,i+1],:])
				updateSettings!(dataFilename,sett,settingsArray,copy(oldsettings.ncolsdfIndata),const_shift_cols,copy(oldsettings.df_name_vector))		
				@assert in(sett.model_type,["build_tree","boosted_tree"]) "Multirow settings are only working for boosting and simple trees (for now)"
				result_of_runmodel,resulting_model=run_model_actual(di["key"],di["trn"],di["numeratortrn"],di["denominatortrn"],di["weighttrn"],di["trn_numfeatures"],di["trn_charfeatures_PDA"],di["keyval"],di["val"],di["val_numfeatures"],di["val_charfeatures_PDA"],di["numeratorval"],di["denominatorval"],di["weightval"],	di["mappings"],datafolder,outfilename,outfileStringOnly,const_shift_cols,sett,di["num_levels"],di["char_levels"],di["all_levels"],di["all_levels_as_string_vector"],di["names_and_levels"],di["candMatWOMaxValues"],dataFilename)
				#get model statistics
				thesem=get_model_stats(resulting_model,defaulted_modelstats_df,i)
			catch
				thesem=deepcopy(defaulted_modelstats_df)		
			end
		append!(allmodelStats,thesem)
		end		
	return allmodelStats
	catch eri
		println(string(eri))
		warn("Model(s) $(iteridx) failed!")
		return defaulted_modelstats_df
	end				
end	

#grid search functions

function grid_search(main_ARGS,minw_list,randomw_list,subsampling_features_prop_list,smoothEstimates_list,niter_list,mf_list,nscores_list,boolRandomizeOnlySplitAtTopNode_list,subsampling_prop_list;metric="Sum MRAE VAL (Raw Estimates)")
	srand(1234)
	settingsFilename,dataFilename,datafolder,outfilename,outfileStringOnly=return_file_and_folders(main_ARGS)
	@assert lowercase(reverse(dataFilename)[1:4])=="dlj." "Grid search requires a *.jld file as input. Abort!"
	
	#settingsFilename="A:\\JuliaTestData\\boostmed.settings.csv"
	settingsArray,number_of_num_features=readSettings(settingsFilename)	
	info("Loading data....")
	@time di=load(dataFilename);
	oldsettings=di["oldsettings"]
	sett=ModelSettings() #default model settings
	const_shift_cols=global_const_shift_cols

	updateSettings!(dataFilename,sett,settingsArray,copy(oldsettings.ncolsdfIndata),const_shift_cols,copy(oldsettings.df_name_vector))
	settlist=Array{ModelSettings}(0) 
	#push!(settlist,deepcopy(sett)) #we do not want/need the initial setting (do we?)	
	#NOTE: it is best to have the minw loop at the lowest level because models with a small minw value will require most memory and these model should thus not be run in parallel but rather sequentially
	for this_randomw in randomw_list,this_subsampling_features_prop in subsampling_features_prop_list,this_smoothEstimates in smoothEstimates_list,this_niter in niter_list,this_mf in mf_list,this_nscores in nscores_list,this_boolRandomizeOnlySplitAtTopNode in boolRandomizeOnlySplitAtTopNode_list,this_subsampling_prop in subsampling_prop_list, this_minw in minw_list
		snew=deepcopy(sett)
		updateSettingsMod!(snew,minw=this_minw,randomw=this_randomw,mf=this_mf,subsampling_features_prop=this_subsampling_features_prop,smoothEstimates=this_smoothEstimates,niter=this_niter,mf=this_mf,nscores=this_nscores,boolRandomizeOnlySplitAtTopNode=this_boolRandomizeOnlySplitAtTopNode,subsampling_prop=this_subsampling_prop)
		push!(settlist,deepcopy(snew))
	end
	info("Number of models to be run: $(length(settlist)). Time $(now())")	
	
	#send data to all workers: without this command we will have a lot of overhead to send the data to the processes; but the data is constant and only needs to be sent ONCE!)
		sendto_module(DTM,workers(),local_data_dict=di)
	
	t0=now()
	#fist iteration 
		i=1		
		v0=@spawn run_mode_grids_single(DataFrame(),t0,i,length(settlist),deepcopy(settlist[i]),dataFilename,local_data_dict,datafolder,outfilename,outfileStringOnly,const_shift_cols)
		v0fetched=fetch(v0) #we want to ensure that this iteration runs properly, before starting work on all workers()
	ncols_stats=size(v0fetched,2)
	@assert ncols_stats>4 "this should not have happened"
	defaulted_modelstats_df=deepcopy(v0fetched)
	deleterows!(defaulted_modelstats_df,1:size(defaulted_modelstats_df,1))
	
	pmapresult=pmap(i->
		      run_mode_grids_single(defaulted_modelstats_df,t0,i,length(settlist),deepcopy(settlist[i]),dataFilename,local_data_dict,datafolder,outfilename,outfileStringOnly,const_shift_cols)
				,2:length(settlist))				
	
	info("Aggregating results...")
		allmodelStats=deepcopy(v0fetched)
		for x in pmapresult
			append!(allmodelStats,x)
		end
		#allmodelStats=reduce(vcat, v0fetched, pmapresult)
		local bestmodels,best10models_idx,best_model_per_iteration
		try
			#find best model according to metrics
			#info("remove this");@time save("c:\\temp\\allmodelStats.jld","allmodelStats",allmodelStats);
			bestmodels=best_models(allmodelStats,metric)
			best10models_idx=bestmodels[1:min(size(bestmodels,1),10),1][:]			
			best_model_per_iteration=best_model_per_iter(allmodelStats,metric)			
		catch
			warn("bestmodels could not be determined!")
			bestmodels=repmat(Any["-1" "-1"],3)
			best_model_per_iteration=repmat(Any["-1" "-1"],3)			
			best10models_idx=[1]
		end
	info("Writing output...")
	#write to Excel file 				
	outf=string(datafolder,"\\",outfileStringOnly,".allstats.xlsx")
	multirowSettingsArray=mapfoldl(x->convert(Array,x)[2,:], vcat, settlist)
	#add header
		multirowSettingsArray=vcat(convert(Array,settlist[1])[1,:],multirowSettingsArray)
	#add index column
		multirowSettingsArray2=hcat(["Model Number";collect(1:size(multirowSettingsArray,1)-1)],multirowSettingsArray)	 
	if size(allmodelStats,1)<1048576 #maxrows of Excel 1<<20
		info("Writing overall model stats file \n $(outf)")
		write_any_array_to_excel_file(vcat(transpose([string(x) for x in names(allmodelStats)]),convert(Array,allmodelStats)),multirowSettingsArray2,bestmodels,best_model_per_iteration,outf)
	else
		outfc=string(datafolder,"\\",outfileStringOnly,".allstats.csv")
		rows_to_keep=findin(allmodelStats[1],best10models_idx) #keep only the best performing models
		unshift!(rows_to_keep,1) #header row must be kept
		allmodelStats_reduced=allmodelStats[rows_to_keep,:]
		info("Writing overall model stats file \n $(outf)")
			write_any_array_to_excel_file(vcat(transpose([string(x) for x in names(allmodelStats_reduced)]),convert(Array,allmodelStats_reduced)),multirowSettingsArray2,bestmodels,best_model_per_iteration,outf)
		info("Writing overall model stats file \n $(outfc)")
			writecsv(outfc,vcat(transpose([string(x) for x in names(allmodelStats)]),convert(Array,allmodelStats)))
		writecsv(string(datafolder,"\\",outfileStringOnly,".allmodelsettings.csv"),multirowSettingsArray2)
		writecsv(string(datafolder,"\\",outfileStringOnly,".bestmodels.csv"),bestmodels)
	end
	info("Done. Time - $(now())")			
	return String["true"],"No Model was returned as this was a multirow run" #need to think about what we return here....
end
	
function run_mode_grids_single(defaulted_modelstats_df::DataFrame,t0::DateTime,i::Int,nmodels::Int,sett::ModelSettings,dataFilename::String,di::Dict,datafolder::String,outfilename::String,outfileStringOnly::String,const_shift_cols::Int)	
	tic();tic();tic();
	try
		outfileStringOnly=string(outfileStringOnly,i)
		outfilename=string(outfilename,i)
		result_of_runmodel,resulting_model=run_model_actual(di["key"],di["trn"],di["numeratortrn"],di["denominatortrn"],di["weighttrn"],di["trn_numfeatures"],di["trn_charfeatures_PDA"],di["keyval"],di["val"],di["val_numfeatures"],di["val_charfeatures_PDA"],di["numeratorval"],di["denominatorval"],di["weightval"],	di["mappings"],datafolder,outfilename,outfileStringOnly,const_shift_cols,sett,di["num_levels"],di["char_levels"],di["all_levels"],di["all_levels_as_string_vector"],di["names_and_levels"],di["candMatWOMaxValues"],dataFilename)
	#get model statistics
		dfModelstats=get_model_stats(resulting_model,defaulted_modelstats_df,i)
		eta=eta_after(t0,i,nmodels);eta_m=round(eta/60,2)	
		telapsed=toq();telapsed=round(telapsed/60,1);
		println("$(now()+Dates.Second(round(eta))) <- ETA -  Id=$(myid()) - Finished #$(i) - ETA: $(round(eta,1))s = $(eta_m)m - Time: $(telapsed)m");
		return dfModelstats
	catch eri
		println(string(eri))
		warn("Model(s) $(i) failed!")
		return defaulted_modelstats_df
	end
end

function get_model_stats(resulting_model::DTModel,defaulted_modelstats_df::DataFrame,i::Int)	
	thesemodelstats=deepcopy(resulting_model.modelstats)	
	nm=DataFrame(ones(Int,size(thesemodelstats,1),1).*i)	
	names!(nm,Symbol[symbol("ModelNumber")])	
	return  thesemodelstats=hcat(nm,thesemodelstats)
end