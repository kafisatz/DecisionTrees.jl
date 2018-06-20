
function dtm(dtmtable::DTMTable,settingsVector::Vector{ModelSettings};file::String=joinpath(mktempdir(),defaultModelNameWtihCSVext))    
    #set rnd state to make this function reporducible (irrespective of trn/val idx)
    #srand should not depend on sett.seed as we do not 'store' the original seed in the resulting Excel file.
    nSettings=length(settingsVector)

    s=hash(dtmtable.numerator,hash(dtmtable.denominator,hash(dtmtable.weight)))
    for x=1:size(dtmtable.features,2)
        s=hash(dtmtable.features[x],s)    
    end
    #intDatahash = Int(.25*hash(2231,s)) 
    intDatahash = Int(.25*hash(2231,hash(dtmtable.features,hash(dtmtable.numerator,hash(dtmtable.denominator,hash(dtmtable.weight))))))
    srand(intDatahash)
    
   #initialize variables
        path_and_fn_wo_extension="some"
        path_and_fn_wo_extension,ext=splitext(file)
        header=Vector{String}(1)
        header_settings=Vector{String}(1)
        i=1        

        #run an initial dummy model to get the 'structure' (ie header and size) of the output statistics
            dummySett=deepcopy(settingsVector[1])
            dummySett.print_details=false;dummySett.bool_write_tree=false;dummySett.write_sas_code=false;dummySett.write_iteration_matrix=false;dummySett.write_result=false;dummySett.write_statistics=false;dummySett.boolCreateZipFile=false;dummySett.write_csharp_code=false;dummySett.write_vba_code=false;dummySett.boolSaveJLDFile=false;dummySett.boolSaveResultAsJLDFile=false;dummySett.showProgressBar_time=false;dummySett.boolProduceEstAndLeafMatrices=false;dummySett.write_dot_graph=false;
            path_and_fn_wo_extension_mod=string(path_and_fn_wo_extension,"_0")
            fnmod=string(path_and_fn_wo_extension_mod,ext)
            dummySett.niter=min(10,dummySett.niter) #is this a good idea? I think we should be fine.
            somestrings,model=run_model_actual(dtmtable,dummySett,fnmod)
            desc,numbrs,desc_settingsvec,settingsvec,selectedPerformanceMeasure=get_stats(model,perfMeasure=settingsVector[1].performanceMeasure)
            #add first column
                pushfirst!(numbrs,i)
                pushfirst!(settingsvec,string(i))
                pushfirst!(desc,"Number")
                pushfirst!(desc_settingsvec,"Number")
            header=deepcopy(desc)
            header_settings=deepcopy(desc_settingsvec)
            allstats=Array{Float64,2}(length(numbrs),nSettings)			
            allsettings=Array{Any,2}(length(settingsvec),nSettings)
            allmodels=Vector{Any}(nSettings)
        #define default stats if the model fails
            defaulted_stats=deepcopy(numbrs)
            defaulted_settings=deepcopy(settingsvec)            
            emptyModel=deepcopy(EmtpyDTModel())
            fill!(defaulted_stats,0.0)

        #create a dictionary which contains all data
            di=Dict("ext"=>ext,"emptyModel"=>emptyModel,"defaulted_settings"=>defaulted_settings,"header_settings"=>header_settings,"header"=>header,"intDatahash"=>intDatahash,"settingsVector"=>deepcopy(settingsVector),"dtmtable"=>dtmtable,"path_and_fn_wo_extension"=>path_and_fn_wo_extension,"defaulted_stats"=>defaulted_stats)
        #send data to all workers: without this command we will have a lot of overhead to send the data to the processes; but the data is constant and only needs to be sent ONCE!)
        #sendto(workers(),local_data_dict=deepcopy(di))
        sendto_module(DecisionTrees,workers(),local_data_dict=deepcopy(di))
        
        #run all models in parallel
        #warn("this is currently rather slow as the local_data_dict might be transferred to each worker for each iteration -> improve this!.... ? global const variable....?")
            pmapresult=pmap(iLoop -> singleRunDtm(iLoop,local_data_dict),1:length(settingsVector))
                       
        #2. run models   
        for pmap_entry in pmapresult# this_sample in cvsampler
            numbrs=pmap_entry[1]
            settingsvec=pmap_entry[2]
            i=pmap_entry[3]
            mdl=pmap_entry[4]
            #append stats 
                allstats[:,i].=deepcopy(numbrs)
                allsettings[:,i].=deepcopy(settingsvec)
                allmodels[i]=mdl
        end
    
        #3. aggregate some statistics
        statsdf=DataFrame(transpose(allstats))    
        settsdf=DataFrame(permutedims(allsettings,[2,1]))
        names!(settsdf,Symbol.(header_settings))
    
        fld,namestr=splitdir(path_and_fn_wo_extension)
        filen=string(path_and_fn_wo_extension,"_multistats.xlsx")
        if isfile(filen)
            @info "Deleting $(filen)."
            rm(filen)
        end
        
        #calc some 'stats of stats'
        col_rng_stats=2:size(statsdf,2)    
        means=map(x->mean(statsdf[x]),col_rng_stats)
        medians=map(x->median(statsdf[x]),col_rng_stats)
        mins=map(x->minimum(statsdf[x]),col_rng_stats)
        maxs=map(x->maximum(statsdf[x]),col_rng_stats)
        stds=map(x->std(statsdf[x]),col_rng_stats)
        sums=map(x->sum(statsdf[x]),col_rng_stats)
    
        stats_of_stats=hcat(means,medians,mins,maxs,stds,sums)
        titles=["Mean","Median","Min","Max","Std Dev","Sum"]
        stats_of_stats=hcat(titles,transpose(stats_of_stats))
        stats_of_stats=vcat(repmat([""],1,size(stats_of_stats,2)),stats_of_stats)
        
        allstats_with_stats=vcat(transpose(allstats),stats_of_stats)
        #NOTE! statsdf is RE defined here!
        statsdf=DataFrame(allstats_with_stats)
        names!(statsdf,Symbol.(header))
        
        #define Exceldata
        sh1=ExcelSheet("settings",settsdf)
        sh2=ExcelSheet("stats",statsdf)
        xld=ExcelData(ExcelSheet[sh1,sh2],Array{Chart,1}(0))
        #write data
        try 		
            println("Writing summary statistics to file:\r\n$(filen)")
            @time write_statistics(xld,filen,true,false)		
        catch e
            @show e       
            warn("DTM: Failed to create Excel Statistics file. \r\n $(filen)")
        end                
        
        return statsdf,settsdf,allmodels
    end
    
        
function singleRunDtm(i::Int,local_data_dict::Dict)	
    di=local_data_dict
    defaulted_settings=di["defaulted_settings"]
    defaulted_stats=di["defaulted_stats"]      
    emptyModel=di["emptyModel"]  
    
    try         
        ext=di["ext"]
        header_settings=di["header_settings"]       
        header=di["header"]       
        intDatahash=di["intDatahash"]
        settingsVector=di["settingsVector"]                
        dtmtable=di["dtmtable"]
        path_and_fn_wo_extension=di["path_and_fn_wo_extension"]
        
        sett=settingsVector[i]
        
        #set seed
            sett.seed = Int(hash(intDatahash,hash(99812,hash(i)))/3)
        #define path
            path_and_fn_wo_extension_mod=string(path_and_fn_wo_extension,"_",i)
            fnmod=string(path_and_fn_wo_extension_mod,ext)
        #run model
            somestrings,model=run_model_actual(dtmtable,sett,fnmod)            
            desc,numbrs,desc_settingsvec,settingsvec,selectedPerformanceMeasure=get_stats(model,perfMeasure=sett.performanceMeasure)
        #add index as first column        
            pushfirst!(numbrs,i)
            pushfirst!(settingsvec,string(i))
            pushfirst!(desc,"Number")
            pushfirst!(desc_settingsvec,"Number")
        #check if header and settings header are as expected
            if !all(header[2:end].==desc[2:end])
                warn("Model run $(i) returned an unexpected results vector:")
                @show desc
                @show header
            end            
            if !all(header_settings[2:end].==desc_settingsvec[2:end])
                warn("Model run $(i) returned an unexpected results vector:")
                @show desc_settingsvec
                @show header_settings
            end            
        return numbrs,settingsvec,i,model
    catch eri
        warn("DTM: Model $(i) failed.")
        println(eri)
        return defaulted_stats,defaulted_settings,i,emptyModel
    end

    return nothing
end	 #singleRunDtm


function get_model_stats(resulting_model::DTModel,defaulted_modelstats_df::DataFrame,i::Int)	
	thesemodelstats=deepcopy(resulting_model.modelstats)	
	nm=DataFrame(ones(Int,size(thesemodelstats,1),1).*i)	
	names!(nm,Symbol[symbol("ModelNumber")])	
	return  thesemodelstats=hcat(nm,thesemodelstats)
end


function createGridSearchSettings(sett::ModelSettings;args...)
    settingsVector=[sett]
	seen=[]
	for (smember,valuelist) in args
		if in(smember,seen)
			error("You provided multiple values for the field $(smember)")
		end
        push!(seen,smember)
        try 
            sz=size(valuelist)
            @assert length(sz)==1 "DTM: expected a vector for $(smember). You provided $(typeof(valuelist)). \r\n$(valuelist)"
            n=sz[1]
            @assert n>0
            
            settingsVectorTemplate=deepcopy(settingsVector)
            
            for i=1:n
                #set value for this 'batch'  
                v=valuelist[i]        
                thisSettingsBatch=deepcopy(settingsVectorTemplate)      
                for kk=1:length(thisSettingsBatch)
                    s=thisSettingsBatch[kk]
                    #updateSettingsMod!(s,smember=)
                    oldValue=getfield(s,smember)
                    if typeof(v)!=typeof(oldValue)
                        valConverted=convertFromString(oldValue,string(v)) #tbd maybe string(v) could/should be replaced with v here
                        setfield!(s,smember,valConverted)
                    else
                        setfield!(s,smember,v)
                    end		
                    res=checkIfSettingsAreValid(s)
                    @assert res==nothing "DTM: Some settings are invalid! Abort."
                end
                if i==1
                    settingsVector=deepcopy(thisSettingsBatch)
                else 
                    append!(settingsVector,thisSettingsBatch)
                end    
            end
        catch e            
            @show smember,size(valuelist)
            @show valuelist
            @show e
            throw("DTM: unable to create settings vector.")
        end
	end
	@assert length(unique(map(x->x.performanceMeasure,settingsVector)))==1 "DTM: Perf measure must be the same for the entire settings Vector"
	println("DTM: Size of settings vector is $(length(settingsVector))")
    return settingsVector
end