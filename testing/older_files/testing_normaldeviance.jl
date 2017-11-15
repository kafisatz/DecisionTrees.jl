#using Gallium
using DTM

#Gallium.breakpoint(calculateSplitValue)

loc="C:\\jdata\\datafolder\\mortgage\\"
rg=[string(loc,"jmortgage.settings_normd.csv") string(loc,"jmortgage.CSV") string("out_jmortgage")];
resbool=DTM.run_model(rg) #small version to compile functions

using JLD
dix=JLD.load(string(loc,"c_temp_myf.jld"))


a=dix[",a"];
labellist=dix["labellist"];
sumnumerator=dix["sumnumerator"];
sumdenominator=dix["sumdenominator"];
sumweight=dix["sumweight"];
countlistfloat=dix["countlistfloat"];
minweight=dix["minweight"];
subs=dix["subs"];
moments_per_pdaclass=dix["moments_per_pdaclass"];


"sumweight",sumweight,"countlistfloat",countlistfloat,"minweight",minweight,"subs",subs,"moments_per_pdaclass",moments_per_pdaclass)
