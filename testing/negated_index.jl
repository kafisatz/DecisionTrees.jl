# allows equivalent of A[-x] in R, assumes x is vector,
# assumes, without checking, that to_remove is sorted!!!!!
# assumes, without checking, that to_keep is precisely least n-length(to_remove) long!
@inbounds function inv_index!(to_keep, to_remove, n::Integer)
    if isempty(to_remove)
        to_keep[1:n] = 1:n
        return to_keep
    elseif length(to_remove) == n
        return to_keep # zero something out?
    end
    prev_r = 0
    to_keepidx = 0
    for ri in 1:length(to_remove)
        r = to_remove[ri]
        if prev_r + 1 != r
            next_to_keep_idx = to_keepidx + (r - prev_r) - 1
            to_keep[(to_keepidx + 1) : next_to_keep_idx] = (prev_r + 1) : (r - 1)
            to_keepidx = next_to_keep_idx
        end
        prev_r = r
    end
    if last(to_remove) != n
        to_keep[(to_keepidx + 1) : (n - length(to_remove))] = (last(to_remove) + 1) : n
    end
    return to_keep
end

using StatsBase
idxpool = zeros(Int64, n);
n=20_000_000;
x=rand(n);

nidx=sample(1:n,1)[1]
to_invert = sort(sample(1:n, nidx, replace=false))
@time idxl=falses(n)
@time for jj in to_invert
    idxl[jj]=true
end
r = n - length(to_invert)
@time ra=inv_index!(view(idxpool, 1:r), to_invert, n); #.28s 232mb
le=x[to_invert]


ff=dtmtable.features[:BRAND_NAME_JY];
po=ff.pool
subset=collect(UInt8(1):UInt8(length(po)))
for i=length(po):-2:1
    deleteat!(subset,i)
end

@btime DTM2.lrIndices(dtmtable.trnidx,ff,subset)
l,r=DTM2.lrIndices(dtmtable.trnidx,ff,subset)


@time ra=inv_index!(view(idxpool, 1:r), to_invert, n); #.28s 232mb
le=x[to_invert]