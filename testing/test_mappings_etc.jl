
ff=dtmtable.features;
mappings=dtmtable.mappings;
candMatWOMaxValues=dtmtable.candMatWOMaxValues;

for i=1:size(ff,2)
if eltype(ff[i])<:String
    @show all(ff[i].pool.==mappings[i])
else
    @show al=all(ff[i].pool.==candMatWOMaxValues[i-this_sett.number_of_char_features])
    if !al
        @show names(ff)[i]
        @show ff[i].pool
        @show candMatWOMaxValues[i-this_sett.number_of_char_features]
    end
end
end

this_sett.number_of_char_features