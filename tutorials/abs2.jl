y=rand();
x=rand();

using BenchmarkTools

f(x,y)=(x-y)^2
f2(x,y)=(x-y)*(x-y)

@btime abs2(x-y)
@btime f(x,y)
@btime f2(x,y)
