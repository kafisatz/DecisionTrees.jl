#this_dir="C:\\Users\\bernhard.konig\\Documents\\ASync\\irobot\\TariffWatch New Repo\\algorithms\\Julia\\Code\\";
#this_dir="R:\\TariffWatch New Repo\\algorithms\\Julia\\Code\\dev\\src"
#unshift!(LOAD_PATH,string(this_dir,"\\DTM2_dev\\src"))
#unshift!(LOAD_PATH,string(this_dir,"\\dev\\src"))
using Compat
using BenchmarkTools
using Revise    
using DecisionTrees
using DataFrames,JLD2
datafolder="C:\\temp\\";outfilename="C:\\temp\\DTM2Model";dataFilename="c:\\temp\\DTM2Model.jld";outfileStringOnly="DTM2Model";fn=string(datafolder,outfileStringOnly)
tmpfilei="C:\\temp\\i.csv";isfile(tmpfilei)&&rm(tmpfilei)
#=
    include("C:\\Users\\bernhard.konig\\Documents\\ASync\\home\\Code\\Julia\\DecisionTrees.jl\\testing\\define_df_prepped.jl");
    df_prepped_orig=deepcopy(df_prepped);
    #df_prepped=df_prepped_orig[1:2000,:];
    df_prepped=deepcopy(df_prepped_orig)
    @time dtmtable=DTM2.prep_data_from_df(df_prepped,this_dt2msett,fn);
=#

di9=load("C:\\Users\\bernhard.konig\\Documents\\ASync\\irobot\\data_tmp\\prepped_dtm2table.jld");dtmtable=di9["dtmtable"];this_dt2msett=di9["this_dt2msett"];

settstr="""model_type="boosted_tree",write_iteration_matrix="false",bool_write_tree="false",write_sas_code="false",write_statistics="false",write_iteration_matrix="false",minw=-.1,randomw=0.0,mf=.1,subsampling_features_prop=1,niter=2,subsampling_prop=1.00,graphvizexecutable="C:\\\\Program Files (x86)\\\\Graphviz2.38\\\\bin\\\\dot.exe",write_result=false"""
eval(Meta.parse(string("DTM2.updateSettingsMod!(this_dt2msett,",settstr,")")))
@time st2,dtm2model=DTM2.run_model_actual(dtmtable,this_dt2msett,fn);

settstr="""model_type="boosted_tree",boolProduceEstAndLeafMatrices=false,niter=50,bool_write_tree="false",write_sas_code="false",write_statistics="false",write_iteration_matrix="false",minw=-.08,randomw=0.0,mf=.1,subsampling_features_prop=1,subsampling_prop=1.00,graphvizexecutable="C:\\\\Program Files (x86)\\\\Graphviz2.38\\\\bin\\\\dot.exe",write_result=false"""
eval(Meta.parse(string("DTM2.updateSettingsMod!(this_dt2msett,",settstr,")")))

Profile.clear_malloc_data()
 st2,dtm2model=DTM2.run_model_actual(dtmtable,this_dt2msett,fn); #2.1s roughly
info("allocation run done.")
quit();



#=
run(`c:\\users\\bernhard.konig\\alloc.bat`)

tmpfilei="C:\\temp\\i.csv";isfile(tmpfilei)&&rm(tmpfilei)
mmfld="C:\\Users\\bernhard.konig\\Documents\\ASync\\irobot\\TariffWatch New Repo\\algorithms\\Julia\\Code\\DTM2_dev\\src\\"
using Coverage
mmm=analyze_malloc(mmfld);idx=map(x->x.bytes>0,mmm);mm2=mmm[idx];
sumKB=sum(map(x->x.bytes/1024, mm2));pct1=map(x->x.bytes/1024/sumKB,mm2);pct2=deepcopy(reverse(pct1));DTM2.incrementalToCumulative!(pct2);reverse!(pct2);
mfiles=map(x->x.filename,mm2);
mfiles_dironly=map(x->splitdir.(x)[1], mfiles);
mfiles_fileonly=map(x->splitdir.(x)[2], mfiles);
kbs=map(x->string(mm2[x].bytes/1024),eachindex(mm2));

#read source code
uq_files=unique(mfiles)
linesdict=Dict{String,Vector{String}}()
ff=uq_files[2]
for i=1:length(uq_files)
    ff=uq_files[i]
    lns=readlines(ff[1:end-4])
    linesdict[splitdir(ff)[2][1:end-4]]=deepcopy(lns)
end

ab=map(x->hcat(string(mm2[x].linenumber),linesdict[mfiles_fileonly[x][1:end-4]][mm2[x].linenumber],mfiles_fileonly[x],kbs[x],string(pct1[x]),string(pct2[x])), eachindex(mm2));writecsv("C:\\temp\\i.csv",ab);
@show "Done. Consider i.csv"
=#