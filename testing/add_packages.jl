for x in ["OnlineStats",
"StatsBase",
"Distributions",
"DataFrames",
"ProgressMeter",
"Iterators",
"ArrayViews",
"HDF5",
"JLD",
"SQLite",
"ProfileView",
"Coverage",
"PyCall"]
    @show x
    Pkg.add(x)
end
Pkg.update();
using WinRPM
WinRPM.update();
