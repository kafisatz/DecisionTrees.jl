n = [3.36, 3.692, 2.05, 22.851, 0.389, 0.351, 0.919, 0.813, 0.055, 6.6];
d =  [1.136, 13.6192, 12.05, 2.851, 9, 10.351, 0.111919, 10.813, 1.055, 1.6];
leftIdx=1:4
rightIdx=5:length(n)
nl=n[leftIdx]
nr=n[rightIdx]
dl=d[leftIdx]
dr=d[rightIdx]

w=copy(d)
wl=w[leftIdx]
wr=w[rightIdx]

include(joinpath("testing","customVariance_fns.jl"))

cv=CustomVariance(n,d)
cv.meanratio
cv.m2
cv.var
cvl=CustomVariance(nl,dl)
cvr=CustomVariance(nr,dr)

merge!(cvl,cvr)

unmerge!(cvl,cvr)
CustomVariance(nl,dl)
#var want 53.3375356025
#m2 want 160.0126068077387

sum(n)/sum(d)
sum(d.*(n./d-sum(n)/sum(d)).^2)
sum(nl)/sum(dl)
sum(nr)/sum(dr)
0