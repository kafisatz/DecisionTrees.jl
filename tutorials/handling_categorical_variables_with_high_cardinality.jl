using Revise
using CSV 
using DecisionTrees  

elt=[Int64,	Float64,	Float64,	Float64,	Float64,	Float64,	Int64,	String,	String,	String,	Int64,	String,	String,	String,	String,	String,	String,	String,	String,	String,	String,	Int64,	String,	Int64,	String,	String,	Int64,	String,	String,	String,	String,	String,	Int64,	String,	String,	String,	String,	Int64,	String,	String,	Int64,	Int64,	String,	Int64,	Int64,	String,	Int64,	Int64,	String,	Int64,	String,	String,	String,	Int64,	String,	String,	String,	String,	Int64,	Int64,	Int64,	String,	String,	String,	String,	String,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64,	Float64,	Int64,	Float64,	Float64,	Float64,	Float64,	Int64,	Int64,	Int64,	Int64,	Int64,	Int64]
datafile="data\\data1small.csv"
@assert isfile(datafile);
@time df_tmp=CSV.read(datafile,types=elt,rows_for_type_detect=20000,allowmissing=:none,categorical=false); 

selected_explanatory_vars=["ART_DES_WOHNEIGENTUM","GEBURTSDATUM","FAMILIENSTAND","NATIONALITAET","PLZ_WOHNORT","GESCHLECHT","FINANZIERUNGSART","STADT","KENNZEICHEN"]

for x in selected_explanatory_vars
    @show x,eltype(df_tmp[Symbol(x)])
end

#consider the unique values of PLZ_WOHNORT
uq=unique(df_tmp[:PLZ_WOHNORT])
length(uq)

#we can analze the frequency with this helper function
counts,freqs,vals,keep30=getCounts(df_tmp[:PLZ_WOHNORT],threshold=30)
#keep is the 30 largest values

#alternatively we can keep all values with are more frequent, than x%
#Let us keep all values more frequent than 1/1000
counts,freqs,vals,keep=getCounts(df_tmp[:PLZ_WOHNORT],threshold=-0.002)

#the we can easily modify the data and overwrite the infrequent values with a custom value

thiscol=df_tmp[:PLZ_WOHNORT]
DecisionTrees.mapToOther!(thiscol,keep30,"some rare PLZ")
unique(df_tmp[:PLZ_WOHNORT])
