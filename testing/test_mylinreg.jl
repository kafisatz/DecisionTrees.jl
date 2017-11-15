(x, w, thisRange) = (
    x=[1.0, 2.0, 3.0]
    w= [0.75, 1.0, 0.75]
    thisRange= [0.0961381, 0.156914, 0.224675]

    DTM2.mylinreg(x,thisRange,w)
    @code_warntype DTM2.mylinreg(x,thisRange,w)
    @btime DTM2.mylinreg(x,thisRange,w)
    @btime DTM2.mylinreg(x,thisRange,w)
    @btime DTM2.mylinreg_alt(x,thisRange,w)

(typeof(x), typeof(w), typeof(thisRange)) = (Array{Float64,1}, Array{Float64,1}, Array{Float64,1})
