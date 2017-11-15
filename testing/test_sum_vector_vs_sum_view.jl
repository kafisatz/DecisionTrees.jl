#this_dir="C:\\Users\\bernhard.konig\\Documents\\ASync\\irobot\\TariffWatch New Repo\\algorithms\\Julia\\Code\\dev\\src";cd(this_dir);push!(LOAD_PATH,this_dir)
using Compat
using BenchmarkTools
using DataFrames, StatsBase
n=20_000_000;
r=rand(n);

idx=sample(1:n,Int(ceil(n*0.7,1)),replace=false,ordered=true);
#idx=[200,2000,5000005];

v=view(r,idx); #
varr=r[idx]; 
@btime view(r,idx); #time to define view
@btime r[idx]; #time to subset vector (copy! / allocation)
@btime sum(v) #time to sum the view
@btime sum(varr) #time to sum the vector

function mys(a)
    z=a[1]
    for j in a
        z+=j
    end
    z
end

@btime mys(v)

pda_sub=pda[idx]
pdaM=pdaMod(pda_sub);
szz=size(v,1)
numeratortot=210.*rand(n);
denominatortot=map(x->x<1.2 ? 1.1 : x,((5.*rand(n)).^3));
weighttot=rand(n);
numerator=numeratortot[idx] #210.*rand(szz);
denominator=denominatortot[idx] #210.*rand(szz);
weight=weighttot[idx] #210.*rand(szz);