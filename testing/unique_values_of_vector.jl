using CSV
using Test 
using DecisionTrees 

elt=[Int64,	Float64,	Float64,	Float64,	Float64,	Float64,	Int64,	String,	String,	String,	Int64,	String,	String,	String,	String,	String,	String,	String,	String,	String,	String,	Int64,	String,	Int64,	String,	String,	Int64,	String,	String,	Int64,	String,	String,	Int64,	String,	String,	Int64,	String,	Int64,	String,	String,	Int64,	Int64,	String,	Int64,	Int64,	String,	Int64,	Int64,	String,	Int64,	String,	String,	String,	Int64,	String,	String,	String,	String,	Int64,	Int64,	Int64,	String,	String,	String,	String,	String,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Float64,	Int64,	Float64,	Float64,	Float64,	Float64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64]

fi="data1medium.csv"
fi="data1small.csv"
@time df_tmp=CSV.read(joinpath(datadir,fi),allowmissing=:none,types=elt,categorical=false,rows_for_type_detect=10000);

selected_explanatory_vars=["PLZ_WOHNORT","ART_DES_WOHNEIGENTUM","GEBURTSDATUM","FAMILIENSTAND","NATIONALITAET","GESCHLECHT","FINANZIERUNGSART","STADT","KENNZEICHEN"]

#for x in selected_explanatory_vars
#    @show x,eltype(df_tmp[Symbol(x)])
#end

dtmtable,sett,df_prepped=prepare_dataframe_for_dtm!(df_tmp,treat_as_categorical_variable=["PLZ_WOHNORT"],weightcol="EXPOSURE",numcol="LOSS20HALF",denomcol="PREMIUM66",independent_vars=selected_explanatory_vars);
sett.minw=-.2


using StatsBase
using BenchmarkTools

v=df_prepped[:PLZ_WOHNORT]
#version 1
@benchmark uq1 = countmap(v)


#version 2
function me(y)
    u=unique(y)
    d=Dict([(i,count(x->x==i,y)) for i in u])
    return d
end
@benchmark uq2=me(v) #this is about 40% faster than countmap. Then again


#sorting / cutoff
function frequencyRel(x)
    d=proportionmap(x)
    srt=sortperm(collect(values(d)))
    return d,srt
end


function getCounts(x;infreq::Float64=-0.01)
    d=countmap(x);    
    vals=collect(values(d));
    sumk=sum(vals);
    if sumk==length(x)
        #ok
    else
       @warn("DTM: Totals do not match. sumv=$(sumk) but length(x)=$(length(x))") 
    end    
    srt=sortperm(vals,rev=true);
    vals=vals[srt];
    valsp=vals./sumk;
    ks=collect(keys(d))[srt];
    return vals,valsp,ks
end
a,b,c=getCounts(v)
0