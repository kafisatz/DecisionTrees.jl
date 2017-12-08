export dtm 

function dtm(dtmtable::DTMTable,sett::ModelSettings,fn::String)
    return run_model_actual(dtmtable::DTMTable,sett::ModelSettings,fn::String)
end

function dtm(dtmtable::DTMTable,sett::ModelSettings,fn::String,cvo::CVOptions)
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
    minw_prop=minw_orig/trnsize_orig
    
#1. sample Data
    size_which_is_sampled=0
    local cvsampler
    if cvo.use_all_data
        #this is the default case
        #trn and val datat is considered for the "sampled trn data"
        size_which_is_sampled=fullsize
    else 
        size_which_is_sampled=trnsize_orig
    end
    if cvo.training_proportion>0
        info("DTM: training_proportion > 0 provided. Performing random sampling without replacment ($(cvo.folds) samples)")
        sample_size=max(1,round(cvo.training_proportion*size_which_is_sampled))
        cvsampler=RandomSub(size_which_is_sampled,sample_size,cvo.folds)
    else
        if cvo.folds<0
            info("DTM: Performing 'DISJOINT' kfold cross validation")
            cvsampler=KfoldDisjoint(size_which_is_sampled,abs(cvo.folds))
        else
            #we have folds>=0
            info("DTM: Performing kfold cross validation")
            cvsampler=Kfold(size_which_is_sampled,cvo.folds)
        end
    end

    #initialize variables
            path_and_fn_wo_extension="some"
            path_and_fn_wo_extension,ext=splitext(fn)
            header=Vector{String}(1)
            header_settings=Vector{String}(1)
            allstats=Array{Float64,2}
            allsettings=Array{Any,2}
            i=1
            n_folds=abs(cvo.folds)
    
    #2. run models   
    for this_sample in cvsampler
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
        dtmtable.trnidx=deepcopy(this_trnidx)
        sett.minw=minw_prop*length(dtmtable.trnidx) #note minw is adjusted for each run!
        
        path_and_fn_wo_extension_mod=string(path_and_fn_wo_extension,"_",i)
        fnmod=string(path_and_fn_wo_extension_mod,ext)
        somestrings,model=run_model_actual(dtmtable,sett,fnmod)
        desc,numbrs,desc_settingsvec,settingsvec=get_stats(model)
        if i==1
            header=deepcopy(desc)
            header_settings=deepcopy(desc_settingsvec)
            allstats=Array{Float64,2}(length(numbrs),n_folds)			
            allsettings=Array{Any,2}(length(settingsvec),n_folds)			
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
    end # this_sample in cvsampler

    #3. aggregate some statistics

    statsdf=DataFrame(transpose(allstats))
    names!(statsdf,Symbol.(header))
    settsdf=DataFrame(permutedims(allsettings,[2,1]))
    names!(settsdf,Symbol.(header_settings))

    fld,namestr=splitdir(path_and_fn_wo_extension)
    filen=string(fld,"multistats.xlsx")
    if isfile(filen)
        info("Deleting $(filen).")
        rm(filen)
    end
    
    #restore indices
    dtmtable.trnidx=deepcopy(trnidx_orig)
    dtmtable.validx=deepcopy(validx_orig)  
    sett.minw=deepcopy(minw_orig)
    
    return statsdf,settsdf
    #return nothing
end