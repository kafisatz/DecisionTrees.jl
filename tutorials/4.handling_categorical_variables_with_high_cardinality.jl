#Date: July 1, 2018
#Author: Bernhard König
#Title: Short example on how categorical variables with high cardinalilty can be treated and preprocessed for these algorithms
#This snippet shall merely show the functionality of the two very simple Julia funcitons mapToOther! and getCounts 
#We refer to the general literature on data preprocessing (specifically for high cardinality categorical data) and feature generation

using Distributed
#Distributed.@everywhere using Revise
Distributed.@everywhere import CSV
Distributed.@everywhere import DataFrames
Distributed.@everywhere import DataFrames: DataFrame

Distributed.@everywhere using DecisionTrees    

elt=[Int,	Float64,	Float64,	Float64,	Float64,	Float64,	Int,	String,	String,	String,	Int,	String,	String,	String,	String,	String,	String,	String,	String,	String,	String,	Int,	String,	Int,	String,	String,	Int,	String,	String,	String,	String,	String,	Int,	String,	String,	String,	String,	Int,	String,	String,	Int,	Int,	String,	Int,	Int,	String,	Int,	Int,	String,	Int,	String,	String,	String,	Int,	String,	String,	String,	String,	Int,	Int,	Int,	String,	String,	String,	String,	String,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Int,	Float64,	Int,	Float64,	Float64,	Float64,	Float64,	Int,	Int,	Int,	Int,	Int,	Int]
datafile=joinpath("data","GermanMotorPremiums","data1small.csv")
@assert isfile(datafile);
@time df_tmp=CSV.read(datafile,types=elt,rows_for_type_detect=20000,allowmissing=:none,categorical=false); 

selected_explanatory_vars=["ART_DES_WOHNEIGENTUM","GEBURTSDATUM","FAMILIENSTAND","NATIONALITAET","PLZ_WOHNORT","GESCHLECHT","FINANZIERUNGSART","STADT","KENNZEICHEN"]

for x in selected_explanatory_vars
    @show x,eltype(df_tmp[Symbol(x)])
end

#consider the unique values of PLZ_WOHNORT
uq=unique(df_tmp[:PLZ_WOHNORT])
length(uq)

#we can analyse the frequency with this helper function
counts,freqs,vals,keep30=getCounts(df_tmp[:PLZ_WOHNORT],threshold=30)
#keep is the 30 largest values

#alternatively we can keep all values with are more frequent, than x%
#Let us keep all values more frequent than 1/1000
counts,freqs,vals,keep=getCounts(df_tmp[:PLZ_WOHNORT],threshold=-0.002)

#the we can easily modify the data and overwrite the infrequent values with a custom value

thiscol=df_tmp[:PLZ_WOHNORT]
DecisionTrees.mapToOther!(thiscol,keep30,"some rare PLZ")
unique(df_tmp[:PLZ_WOHNORT])
