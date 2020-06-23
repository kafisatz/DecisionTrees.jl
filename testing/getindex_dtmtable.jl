idx = [4,5,8,10]
one_to_n = collect(1:10)
trnidx = [2,4,8,10]
validx = [1,3,5,6,7,9]

# found=findall((in)(idx),trnidx)
# integersNew=trnidx[found]
# trnidxNEW=findall((in)(integersNew),idx)
found = findall((in)(idx), trnidx)
found = trnidx[found]
trnidxNEW = findall((in)(found), idx)

found = findall((in)(idx), validx)
integersNew = validx[found]
validxNEW = findall((in)(integersNew), idx)


# function Base.getindex(r::DTMTable,idx::Vector{T};compact_features=false) where T <: Number
	@assert !compact_features 
	key = r.key[idx]
numerator = r.numerator[idx]    
denominator = r.denominator[idx]    
weight = r.weight[idx]    
	features = r.features[idx,:]

	# trn and val idx:
	one_to_n = collect(1:length(r.weight))
	rows_to_keep = idx # one_to_n[idx]        
	foundtrn = findall((in)(rows_to_keep), r.trnidx)
	foundval = findall((in)(rows_to_keep), r.validx)	
	new_trnidx = r.trnidx[foundtrn]
new_validx = r.validx[foundval]
    
	(length(new_validx) != 0) || @warn("DTM: Validation index of the DTMTable has length zero! You may need to redefine validx. The function resample_trnvalidx! might be helpful.")
	(length(new_trnidx) != 0) || @warn("DTM: Training index of the DTMTable has length zero! You may need to redefine trnidx. The function resample_trnvalidx! might be helpful.")
	
#	
#	dt=DTMTable(key,new_trnidx,new_validx,numerator,denominator,weight,features,deepcopy(r.candMatWOMaxValues),deepcopy(r.mappings))
#   return dt
# end