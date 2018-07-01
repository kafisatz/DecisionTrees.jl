#Date: July 1, 2018
#Author: Bernhard König
#Title: An application of Regression Trees for Reserving 

#We have used the code in RCode\GenerateIndividualClaimsData.R 
#to generate the data Simulated.Cashflows.csv
#We refer to the work of Wuethrich and Gabrielli for details on the
#simulated claims data, see the following references:

### Source: https://people.math.ethz.ch/~wmario/simulation.html
### Reference: Individual Claims History Simulation Machine. ###
###            SSRN Manuscript ID 3130560.                   ###  
###  Authors: Andrea Gabrielli, Mario V. Wuthrich    ###

import Distributed
import Random  
import DelimitedFiles
Distributed.@everywhere import CSV
Distributed.@everywhere import DataFrames
Distributed.@everywhere import DataFrames: DataFrame,groupby,combine,names!

Distributed.@everywhere using DecisionTrees

##############################
#Read the data
##############################
datafile=joinpath("data","reservingAccident","Simulated.Cashflows.csv")
@assert isfile(datafile) "File not found: $(datafile). Consider the readme.txt in the folder data\\reservingAccident\\"
@time fullData=CSV.read(datafile,rows_for_type_detect=100000,allowmissing=:none,categorical=false);

#The data should match this total for the column Paycum11
@assert sum(fullData[:PayCum11])==1464014219 "Expected 1464014219, have: $(sum(fullData[:PayCum11]))"
#The data is described in the paper mentioned above

#We provide a short summary from the Readme of Wuethrich and Gabrielli
#The variable ClNr is the claim number, this is a unique numerical identifier for each individual claim. 
#LoB is the line of business which should be of factor type (categorical) and may take the labels 1; : : : ; 4. 
#cc is the claims code that should be of factor type (categorical) and may take the labels 1; : : : ; 53 (with possible gaps).
#AY is the year of claims occurrence (accident year) that should be of integer type and may take the values 1994; : : : ; 2005. 
#AQ is the quarter of claims occurrence (accident quarter) that should be of integer type and may take the values 1; : : : ; 4.
#age is the age of the injured that should be of integer type and may take the values 15; : : : ; 70.
#inj_part is the body part injured that should be of factor type (categorical) and may take the labels 1; : : : ; 99 (with gaps).

#set independent variables
selected_explanatory_vars=["LoB","age","cc","inj_part"]
#note: we refrain from using the exposure as explanatory variable, as this is not meaningful in practical applications
#where one aims to predict the claim frequency of a given policy (without knowing how long that policy will remain within a certain portfolio)

#check type of each column
for x in selected_explanatory_vars
    println(x)
    println(eltype(fullData[Symbol(x)]))
    println("first ten values")
    print(fullData[1:10,Symbol(x)])
    println("")
end

#consider unique elements for each variable
for x in 1:size(fullData,2)
    println(x)
    @show size(unique(fullData[:,x]))
end

##############################
#Prepare the data
##############################
#consider the Accient Years
ayears=unique(fullData[:AY])
#For the first model we perform we consider the Cumulative Payment after 1 year and compare it to the Cumulative Payment after 2 years
dataWOLastYear=copy(fullData)
#Discard most recent AY
dataWOLastYear=dataWOLastYear[dataWOLastYear[:AY].<2005,:]
#Discard claims which had 0 cumulative payment by year 1
dataWOLastYear=dataWOLastYear[dataWOLastYear[:PayCum01].>0,:]
dtmtable,sett,dfprepped=prepare_dataframe_for_dtm!(dataWOLastYear,keycol="ClNr",numcol="PayCum00",denomcol="PayCum01",independent_vars=selected_explanatory_vars,treat_as_categorical_variable=["inj_part","cc"]);

#Define a random split into training and validaiton data 
#training shall be 70% of the data
srand(1239)
resample_trnvalidx!(dtmtable,0.7);
#dtmtable.trnidx[1:10]

#Consider the explanatory variables
dtmtable.features
#There are four variables
#One could possibly add the accident year and quarter as addtitional explanatory variables. However, depending on the problem setting this may not be meaningful.

updateSettingsMod!(sett,minw=-0.15,model_type="build_tree")
#run a single tree 
resultingFiles,resM=dtm(dtmtable,sett)

#the chain ladder volume weighted LDF for Y1-Y2 is 1.62458536172099 
#its inverse is 0.615541678241308 

@warn("show scores y1/y2 ldf versus y2/y3 ldf in an xy scatter plot")

function buildTriangle(data;rowIdentifier=:AY,columns=[Symbol(string("PayCum",lpad(i,2,0))) for i=0:11])
    minRow,MaxRow=extrema(fullData[rowIdentifier])
    nRows=MaxRow-minRow+1
    nColumns=length(columns)
    res=zeros(Int,nRows,nColumns)
    have=names(fullData)
    @assert issubset(columns,have)
    for i=1:size(fullData,1)
        for j=1:length(columns)     
            @inbounds thisAY = fullData[rowIdentifier][i]
            @inbounds res[thisAY-minRow+1,j] += fullData[columns[j]][i]
        end  
    end
    return res
end

@warn("add mini CL code here")

@time triangle=buildTriangle(fullData)
DelimitedFiles.writedlm("C:\\temp\\1.csv",triangle,',')
