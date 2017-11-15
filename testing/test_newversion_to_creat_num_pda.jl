unshift!(LOAD_PATH,"C:\\Users\\bernhard.konig\\Documents\\ASync\\irobot\\TariffWatch New Repo\\algorithms\\Julia\\Code\\DTM2_dev\\src")
using DTM2

n=2000;
x=rand(n,20);

this_column=view(x,:,3)
max_splitting_points_num=6

candlist=DTM2.define_candidates(this_column,max_splitting_points_num)
if max_splitting_points_num>500
    warn("max_splitting_points_num=$(max_splitting_points_num) is larger than 500. That may not be a good idea")
    if max_splitting_points_num>2000
        error("DTM: max_splitting_points_num too large max_splitting_points_num=$(max_splitting_points_num)")
    end
end
(length(candlist)==1)&&(info("Numeric column $(i):$(header[i]) (numeric) has zero splitting points.")) #throw an error when not splitting point exists, was this intended? does this ever happen?
sort!(candlist,alg=QuickSort) #although it should already be sorted
mini,maxi=extrema(this_column) #it IS CRUCIAL that the maximum is taken from the whole vector here! (and not only the training part)
push!(candMatWOMaxValues,deepcopy(candlist))
candlist[end]=max(candlist[end],maxi) 

this_column
candlist
mapped_vec=zeros(this_column)
for ii=1:length(this_column)
    #x=this_column[1]
    mapped_idx=searchsortedfirst(candlist,this_column[ii])
    @show xnew=candlist[mapped_idx]
    mapped_vec[ii]=xnew
end

mapped_vec
sort(unique(mapped_vec))==candlist
