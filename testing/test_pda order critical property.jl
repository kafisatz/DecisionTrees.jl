

#test if the PooldDataArray is constructed the way I think it is
# f.pool
ft=dtmtable.features[2]
ft.pool
unique(ft)

for j=1:size(dtmtable.features,2)
    ft=dtmtable.features[j]
    for i=1:length(ft.pool)
        @assert all(ft[ft.refs.==UInt8(i)].==ft.pool[i])
    end
end