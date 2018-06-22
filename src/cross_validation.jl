export dtm,dtm_debug,dtm_single_threaded

#note: these functions (in cross_validation.jl have major redundancies with runmodel_multi.jl
#we may want to improve this

function dtm_single_threaded(dtmtable::DTMTable,sett::ModelSettings,cvo::CVOptions;file::String=joinpath(mktempdir(),defaultModelNameWtihCSVext))
    @info "DTM: (bk) You may want to consider the multi threaded verison of this function: simply call dtm(...)"
#if folds<0 then we consider n disjoint training sets

#currently we have trhee possible CVs
#k disjoint samples 
#    KfoldDisjoint(length(weight),k) 
#random sampling, k samples of size samplesize
#    RandomSub(length(weight),samplesize,k)
#k samples where the k validation data sets are a partition
#    Kfold(length(weight),k)
    
    #during CV  we will overwrite trnidx and validx!
    #store orig indices
    trnidx_orig=deepcopy(dtmtable.trnidx)
    trnsize_orig=length(trnidx_orig)
    fullsize=length(dtmtable.weight)
    validx_orig=deepcopy(dtmtable.validx)
    minw_orig=deepcopy(sett.minw)
    seed_orig=deepcopy(sett.seed)
    minw_prop=minw_orig/trnsize_orig
    
#set rnd state to make this function reporducible (irrespective of trn/val idx)
#srand should not depend on sett.seed as we do not 'store' the original seed in the resulting Excel file.
	intDatahash = floor(Int,.25*hash(2231,hash(dtmtable.features,hash(dtmtable.numerator,hash(dtmtable.denominator,hash(dtmtable.weight))))))
    srand(intDatahash)
    
#1. sample Data
    size_which_is_sampled=0
    #local cvsampler
    if cvo.use_all_data
        #this is the default case
        #trn and val datat is considered for the "sampled trn data"
        size_which_is_sampled=fullsize
    else 
        size_which_is_sampled=trnsize_orig
    end
    if cvo.training_proportion>0
        @info "DTM: training_proportion > 0 provided. Performing random sampling without replacement ($(cvo.folds) samples)"
        sample_size=max(1,round(cvo.training_proportion*size_which_is_sampled))
        cvsampler=RandomSub(size_which_is_sampled,sample_size,cvo.folds)
    else
        if cvo.folds<0
            @info "DTM: Performing 'DISJOINT training sets' k-fold cross validation (k=$(cvo.folds))"
            cvsampler=KfoldDisjoint(size_which_is_sampled,abs(cvo.folds))
        else
            #we have folds>=0
            @info "DTM: Performing k-fold cross validation (disjoint validation sets, k=$(cvo.folds))"
            cvsampler=Kfold(size_which_is_sampled,cvo.folds)
        end
    end

    #initialize variables
            path_and_fn_wo_extension="some"
            path_and_fn_wo_extension,ext=splitext(file)
            header=Vector{String}(undef,1)
            header_settings=Vector{String}(undef,1)
            allstats=Array{Float64,2}
            allsettings=Array{Any,2}
            allmodels=Vector{Any}
            i=1
            n_folds=abs(cvo.folds)

    #2. run models   
    for this_sample in cvsampler
        sett.seed = floor(Int,0.31*hash(intDatahash,hash(99812,hash(i))))
        if cvo.use_all_data
            #size_which_is_sampled=fullsize            
            this_trnidx=this_sample
            dtmtable.validx=setdiff(1:size_which_is_sampled,this_sample)
        else 
            #size_which_is_sampled=trnsize_orig
            #here the training data is a subset of the trn orig data
            this_trnidx=trnidx_orig[this_sample]
            #here we leave the val data as it is
            #notably trn and val "do not fill out" all data, there is third part of the data which remains unused here
        end
        dtmtable.trnidx=this_trnidx
        if sett.minw>0 #if sett.minw is specified as a negative number then it will be interpreted as a percentage anyway (no need to update it)
            sett.minw=minw_prop*length(dtmtable.trnidx) #note minw is adjusted for each run!
        end
        
        path_and_fn_wo_extension_mod=string(path_and_fn_wo_extension,"_",i)
        fnmod=string(path_and_fn_wo_extension_mod,ext)
        somestrings,model=run_model_actual(dtmtable,sett,fnmod)        
        desc,numbrs,desc_settingsvec,settingsvec,selectedPerformanceMeasure=get_stats(model,perfMeasure=sett.performanceMeasure)
        
        #add index as first column        
            pushfirst!(numbrs,i)
            pushfirst!(settingsvec,string(i))
            pushfirst!(desc,"Number")
            pushfirst!(desc_settingsvec,"Number")        
        if i==1
            header=deepcopy(desc)
            header_settings=deepcopy(desc_settingsvec)
            allstats=Array{Float64,2}(length(numbrs),n_folds)			
            allsettings=Array{Any,2}(length(settingsvec),n_folds)
            allmodels=Vector{Any}(n_folds)
        else			
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
        end
        allstats[:,i].=deepcopy(numbrs)
        allsettings[:,i].=deepcopy(settingsvec)
        allmodels[i]=model      
        i+=1
    end # this_sample in cvsampler

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
    stds=map(x->StatsBase.std(statsdf[x]),col_rng_stats)
    sums=map(x->sum(statsdf[x]),col_rng_stats)

    stats_of_stats=hcat(means,medians,mins,maxs,stds,sums)
    titles=["Mean","Median","Min","Max","Std Dev","Sum"]
    stats_of_stats=hcat(titles,transpose(stats_of_stats))
    stats_of_stats=vcat(repeat([""],1,size(stats_of_stats,2)),stats_of_stats)
    
    allstats_with_stats=vcat(transpose(allstats),stats_of_stats)
    #NOTE! statsdf is RE defined here!
    statsdf=DataFrame(allstats_with_stats)
    names!(statsdf,Symbol.(header))
    
    #define Exceldata
    sh1=ExcelSheet("settings",settsdf)
    sh2=ExcelSheet("stats",statsdf)
    xld=ExcelData(ExcelSheet[sh1,sh2],Array{Chart,1}(undef,0))
    #write data
    @info "DTM: Writing multistats to file:\r\n$(filen)"
    try 		
		@time write_statistics(xld,filen,true,false)		
	catch e
        @show e       
		warn("DTM: Failed to create Excel Statistics file. \r\n $(filen)")
	end        
    #restore indices
    dtmtable.trnidx=deepcopy(trnidx_orig)
    dtmtable.validx=deepcopy(validx_orig)  
    sett.minw=deepcopy(minw_orig)
    sett.seed = seed_orig 
    
    return statsdf,settsdf,allmodels 
end

function dtm(dtmtable::DTMTable,sett::ModelSettings,cvo::CVOptions;file::String=joinpath(mktempdir(),defaultModelNameWtihCSVext))
    #warn("This requires testing especially when Distributed.nprocs()>2")
    #if folds<0 then we consider n disjoint training sets
    
    #currently we have three possible CVs
    #k disjoint samples 
    #    KfoldDisjoint(length(weight),k) 
    #random sampling, k samples of size samplesize
    #    RandomSub(length(weight),samplesize,k)
    #k samples where the k validation data sets are a partition
    #    Kfold(length(weight),k)
        
    #during CV  we will overwrite trnidx and validx!
    #store orig indices
    trnidx_orig=deepcopy(dtmtable.trnidx)
    trnsize_orig=length(trnidx_orig)
    fullsize=length(dtmtable.weight)
    validx_orig=deepcopy(dtmtable.validx)
    minw_orig=deepcopy(sett.minw)
    seed_orig=deepcopy(sett.seed)
    minw_prop=minw_orig/trnsize_orig

    #set rnd state to make this function reporducible (irrespective of trn/val idx)
    #srand should not depend on sett.seed as we do not 'store' the original seed in the resulting Excel file.
    s=hash(dtmtable.numerator,hash(dtmtable.denominator,hash(dtmtable.weight)))
    for x=1:size(dtmtable.features,2)
        s=hash(dtmtable.features[x],s)    
    end
    intDatahash = Int(.25*hash(2231,s)) 
     # the next line does not work because (at the time of writing) the DataFrames Pkg has not implemented hash(x::DataFrame,u::UInt) correctly.
     # intDatahash = Int(.25*hash(2231,hash(dtmtable.features,hash(dtmtable.numerator,hash(dtmtable.denominator,hash(dtmtable.weight))))))
    srand(intDatahash)
    
    #1. sample Data
        size_which_is_sampled=0
        #local cvsampler
        if cvo.use_all_data
            #this is the default case
            #trn and val datat is considered for the "sampled trn data"
            size_which_is_sampled=fullsize
        else 
            size_which_is_sampled=trnsize_orig
        end
        if cvo.training_proportion>0
            @info "DTM: training_proportion > 0 provided. Performing random sampling without replacement ($(cvo.folds) samples)"
            sample_size=max(1,round(cvo.training_proportion*size_which_is_sampled))
            cvsampler=RandomSub(size_which_is_sampled,sample_size,cvo.folds)
        else
            if cvo.folds<0
                @info "DTM: Performing 'DISJOINT training sets' k-fold cross validation (k=$(cvo.folds))"
                cvsampler=KfoldDisjoint(size_which_is_sampled,abs(cvo.folds))
            else
                #we have folds>=0
                @info "DTM: Performing k-fold cross validation (disjoint validation sets, k=$(cvo.folds))"
                cvsampler=Kfold(size_which_is_sampled,cvo.folds)
            end
        end
    
        #initialize variables
        path_and_fn_wo_extension="some"
        path_and_fn_wo_extension,ext=splitext(file)
        header=Vector{String}(undef,1)
        header_settings=Vector{String}(undef,1)
        i=1
        n_folds=abs(cvo.folds)

        #run an initial dummy model to get the 'structure' (ie header and size) of the output statistics
            dummySett=deepcopy(sett)
            dummySett.print_details=false;dummySett.bool_write_tree=false;dummySett.write_sas_code=false;dummySett.write_iteration_matrix=false;dummySett.write_result=false;dummySett.write_statistics=false;dummySett.boolCreateZipFile=false;dummySett.write_csharp_code=false;dummySett.write_vba_code=false;dummySett.boolSaveJLDFile=false;dummySett.boolSaveResultAsJLDFile=false;dummySett.showProgressBar_time=false;dummySett.boolProduceEstAndLeafMatrices=false;dummySett.write_dot_graph=false;                
            path_and_fn_wo_extension_mod=string(path_and_fn_wo_extension,"_0")
            fnmod=string(path_and_fn_wo_extension_mod,ext)
            somestrings,model=run_model_actual(dtmtable,dummySett,fnmod)
            desc,numbrs,desc_settingsvec,settingsvec,selectedPerformanceMeasure=get_stats(model,perfMeasure=sett.performanceMeasure)
            #add first column
                pushfirst!(numbrs,i)
                pushfirst!(settingsvec,string(i))
                pushfirst!(desc,"Number")
                pushfirst!(desc_settingsvec,"Number")
            header=deepcopy(desc)
            header_settings=deepcopy(desc_settingsvec)
            allstats=Array{Float64,2}(undef,length(numbrs),n_folds)			
            allsettings=Array{Any,2}(undef,length(settingsvec),n_folds)
            allmodels=Vector{Any}(undef,n_folds)
        #define default stats if the model fails
            defaulted_stats=deepcopy(numbrs)
            defaulted_settings=deepcopy(settingsvec)
            fill!(defaulted_stats,0.0)

        #create a dictionary which contains all data
            di=Dict("ext"=>ext,"minw_prop"=>minw_prop,"size_which_is_sampled"=>size_which_is_sampled,"defaulted_settings"=>defaulted_settings,"header_settings"=>header_settings,"header"=>header,"cvsampler"=>cvsampler,"intDatahash"=>intDatahash,"sett"=>deepcopy(sett),"cvo"=>cvo,"dtmtable"=>dtmtable,"path_and_fn_wo_extension"=>path_and_fn_wo_extension,"defaulted_stats"=>defaulted_stats)
        #send data to all workers: without this command we will have a lot of overhead to send the data to the processes; but the data is constant and only needs to be sent ONCE!)        
        sendto_module(DecisionTrees,Distributed.workers(),local_data_dict=deepcopy(di))
        
        #run all models in parallel
        warn("this is currently terribly slow as the local_data_dict might be transferred to each worker for each iteration -> improve this!.... ? global const variable....?")
            pmapresult=Distributed.pmap(iLoop -> run_cvsample_on_a_process(iLoop,local_data_dict),1:length(cvsampler))
                       
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
        stds=map(x->StatsBase.std(statsdf[x]),col_rng_stats)
        sums=map(x->sum(statsdf[x]),col_rng_stats)
    
        stats_of_stats=hcat(means,medians,mins,maxs,stds,sums)
        titles=["Mean","Median","Min","Max","Std Dev","Sum"]
        stats_of_stats=hcat(titles,transpose(stats_of_stats))
        stats_of_stats=vcat(repeat([""],1,size(stats_of_stats,2)),stats_of_stats)
        
        allstats_with_stats=vcat(transpose(allstats),stats_of_stats)
        #NOTE! statsdf is RE defined here!
        statsdf=DataFrame(allstats_with_stats)
        names!(statsdf,Symbol.(header))
        
        #define Exceldata
        sh1=ExcelSheet("settings",settsdf)
        sh2=ExcelSheet("stats",statsdf)
        xld=ExcelData(ExcelSheet[sh1,sh2],Array{Chart,1}(undef,0))
        #write data
        try 		
            @time write_statistics(xld,filen,true,false)		
        catch e
            @show e       
            warn("DTM: Failed to create Excel Statistics file. \r\n $(filen)")
        end        
        #restore indices
        dtmtable.trnidx=deepcopy(trnidx_orig)
        dtmtable.validx=deepcopy(validx_orig)  
        sett.minw=deepcopy(minw_orig)
        sett.seed = seed_orig 
        
        return statsdf,settsdf,allmodels 
    end
        
function run_cvsample_on_a_process(i::Int,local_data_dict::Dict)	
    di=local_data_dict
    defaulted_settings=di["defaulted_settings"]
    defaulted_stats=di["defaulted_stats"]        
    
    try         
        ext=di["ext"]
        minw_prop=di["minw_prop"]
        size_which_is_sampled=di["size_which_is_sampled"]
        header_settings=di["header_settings"]       
        header=di["header"]       
        cvsampler=di["cvsampler"]       
        intDatahash=di["intDatahash"]
        sett=di["sett"]
        cvo=di["cvo"]
        dtmtable=di["dtmtable"]
        path_and_fn_wo_extension=di["path_and_fn_wo_extension"]
                
        #find ith sample
            jj=1
            local this_sample
            for a_sample in cvsampler            
                jj+=1
                if jj>=i
					this_sample=a_sample # this line is new; on 0.7alpha Julia was complainig that I use the value of the iterator variable which may not exist outside the loop in later verions (although I clearly defined this_sample outside the loop....)
                    break
                end
            end
        #define trnidx
            if cvo.use_all_data
                #size_which_is_sampled=fullsize            
                this_trnidx=this_sample
                dtmtable.validx=setdiff(1:size_which_is_sampled,this_sample)
            else                 
                #here the training data is a subset of the trn orig data
                this_trnidx=trnidx_orig[this_sample]
                #here we leave the val data as it is;notably trn and val "do not fill out" all data, there is third part of the data which remains unused here
            end
            dtmtable.trnidx=this_trnidx
        #update minw
            if sett.minw>0 #if sett.minw is specified as a negative number then it will be interpreted as a percentage anyway (no need to update it)
                sett.minw=minw_prop*length(dtmtable.trnidx) #note minw is adjusted for each run!
            end
        #set seed
            sett.seed =  floor(Int,0.31*hash(intDatahash,hash(99812,hash(i))))
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
        warn("DTM: CV model $(i) failed.")
        println(eri)
        return defaulted_stats,defaulted_settings,i,deepcopy(EmtpyDTModel())
    end

    return nothing
end	 #run_cvsample_on_a_process