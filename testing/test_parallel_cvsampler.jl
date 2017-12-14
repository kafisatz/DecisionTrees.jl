import DecisionTrees.Kfold


cvs=Kfold(10,4)


function mura(c)   
    return c,3,99.1
end

pmapres=pmap(x->mura(x),cvs)


mura(c)
collect(c)