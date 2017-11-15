include("analyze_pw_only_data_init_only.jl")

niters=[200]
mfs=[0.05 0.03] #test 0.05 and 0.03 test mf 0.1 0.05 and 0.15  #0.13 is no good
randomws=[0.0]
minws=[-0.08 -0.07]
subsampling_features_props=[0.6666]
subsampling_props=[1.0]
smoothEstimatess=[1] #smooth =1 seem SLIGHTLY BETTER
my_count=0

for niter in niters
for mf in mfs #[.03 .05] #unclear what about .08?
for randomw in randomws #[0.0 .02]  #preference 0
for minw in minws #[-.075] #preference -.75
for subsampling_features_prop in subsampling_features_props #[1.0] # .7 .5]
for subsampling_prop in subsampling_props #[.6666 1.0] # yet unclear which is the best paramter here.... #preference .6666
for smoothEstimates in smoothEstimatess #[1 0] #we have not checke this yet!
		try
			outfileStringOnly=string("boost_it",niter,"_mf",mf,"_rand",randomw,"_minw",minw,"_subsp",subsampling_features_prop,"_subsampl",subsampling_prop,"_smooth",smoothEstimates);
			updateSettingsMod!(this_sett,minw=minw,randomw=randomw,mf=mf,niter=niter,subsampling_features_prop=subsampling_features_prop,subsampling_prop=subsampling_prop,smoothEstimates=smoothEstimates,boolRandomizeOnlySplitAtTopNode="true",model_type="boosted_tree",scorebandsstartingpoints="1,50,100,150,200,300,400,500,600,700,800,850,900,950",showTimeUsedByEachIteration=1,moderationvector="$(mf)");
			@time result_of_runmodel,resulting_model=run_model_actual(key,trn,numeratortrn,denominatortrn,weighttrn,trn_numfeatures,trn_charfeatures_PDA,keyval,val,val_numfeatures,val_charfeatures_PDA,numeratorval,denominatorval,weightval,mappings,datafolder,outfilename,outfileStringOnly,DTM.global_const_shift_cols,sett,num_levels,char_levels,all_levels,all_levels_as_string_vector,names_and_levels,candMatWOMaxValues,dataFilename);
			nothing
		catch
			warn("BK: this model did not work: $(my_count)")
		end
	my_count+=1
end
end
end
end
end
end

finished=now()
totaltime=finished-time_x_started
totaltime_m=totaltime/convert(typeof(totaltime), Base.Dates.Minute(1))
@show totaltime_m

@show my_count;
