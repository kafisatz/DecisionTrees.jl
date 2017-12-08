#First, install the most recent Julia version

Pkg.update() #Upadate all packages
#Pkg.add("MySQL") #this may not work
#warn("You need to install mysql-connector-c-6.1.6-winx64 !!") #-> resources folder
#warn("ensure that C:\Program Files\MySQL\MySQL Connector C 6.1\lib is in your path such that libmysql.dll is found")
#Pkg.clone("https://github.com/JuliaDB/MySQL.jl")
Pkg.add("RDatasets") #I do not think we really need this package
Pkg.add("CSV")
Pkg.add("ProfileView")
Pkg.add("Coverage")
Pkg.add("PooledArrays")

Pkg.add("Revise")
Pkg.add("BenchmarkTools")
Pkg.add("Compat")
Pkg.add("StatsBase")
Pkg.add("DataFrames")
Pkg.add("ProgressMeter")
Pkg.add("Iterators")
Pkg.add("OnlineStats")
Pkg.add("JLD2")
Pkg.add("SQLite")
Pkg.add("Loess")
Pkg.add("PyCall") #Will be used to export to *.xlsx and create graphs

#= 
#if pycall is not working, you may want to try this...
warn("Follow this step by step:")
warn("this may be outdated, newer versions of python will probably work too")
info("FIRST Install Anaconda with Python 2.7")
info("Then start an anaconda command prompt and type this")
# "conda install mingw libpython" #Comment: mingw seems to be a requirement (it will be installed automatically with libpython)
info("then add an environment Variable PYTHON which points to the python.exe")
info("then perform a reboot of the computer and add the PyCall package")
#After the reboot in Julia ENV["PYTHON"] should point to the executable

#ENV["PYTHON"]="C:\\JuliaPro-0.5.1.1\\Python\\python.exe"
Pkg.build("PyCall")
warn("install package xlsxwriter for python")
warn("install package pandas for python")

info("ensure that python, the correct verison, is in the path")
warn("run this command in windows to install pandas "python -m pip install pandas")

#For debugging purposes only

Pkg.add("SortingAlgorithms") #I do not think we really need this package
Pkg.add("ZipFile") #I do not think we really need this package
#Pkg.add("Winston")

#Pkg.add("Dates") #this will be included in base in Julia v"0.4"

using WinRPM
WinRPM.update()
Pkg.update() #Upadate all packages

#Also I would suggest to open *.jl files with Notepad++
#and install the Autohotkey plugin
=#

#=
	#you should test this
	Pkg.build("PyCall")
	using PyCall
=#
