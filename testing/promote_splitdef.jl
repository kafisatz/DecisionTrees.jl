

x=Splitdef(1,2,:mi,UInt8[1 2 33][:],2.3,3.1,99.2)

z=Splitdef{UInt16}(x)

#todo check the difference between immutable and type!! are we using immutables correctly here?
struct Splitdef{T<:Unsigned} #Splitdef(feature_column_id,feature_column_id2,fname,tmp_result[2],tmp_result[1],tmp_result[3],tmp_result[4])]
	featid::Int64
	featid_new_positive::Int64
    featurename::Symbol
    subset::Array{T,1}
    splitvalue::Float64 #Depends on the criterion function used (e.g. _difference)
    weightl::Float64
    weightr::Float64
end