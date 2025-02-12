# sorting.jl

using BenchmarkTools

nn = 600000

raw_rel = rand(nn);
raw_rel1 = deepcopy(raw_rel);
raw_rel2 = deepcopy(raw_rel);
srt = zeros(Int, nn);

srt = zeros(Int, nn);
@benchmark sortperm!(srt, raw_rel1, alg=QuickSort)
srt2 = zeros(Int, nn);
@benchmark DecisionTrees.my_sortperm_you_should_regularily_check_if_sortperm_in_base_has_become_more_efficient!(srt2, raw_rel) # we could use RadixSort here, it uses about 40m of memory for 5m obs, AND I am unsure how easy it is to get the permutation (compared to sorting in place). But the algorithm is about 25% faster for random floats than QuickSort.    
nothing


    
function tt(n)
    t1 = 0.0
    t2 = 0.0
    for i = 1:n
        srt = zeros(Int, nn);
        a = @timed sortperm!(srt, raw_rel1, alg=QuickSort)
        srt2 = zeros(Int, nn);
        b = @timed my_sortperm_you_should_regularily_check_if_sortperm_in_base_has_become_more_efficient!(srt2, raw_rel) # we could use RadixSort here, it uses about 40m of memory for 5m obs, AND I am unsure how easy it is to get the permutation (compared to sorting in place). But the algorithm is about 50% faster for random floats than QuickSort.    
        t1 += a[2]
        t2 += b[2]
    end
    return t1, t2
end