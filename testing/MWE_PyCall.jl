#using DecisionTrees 
using PyCall
using DataFrames

statsfile = joinpath(mktempdir(),"abc.xlsx")
#statsfile = raw"C:\temp\8.xlsx"

function me(statsfile) 
data = DataFrame(x=collect(1:13), y=collect(1:13).*2)

@assert isdir(splitdir(statsfile)[1])
isfile(statsfile) && rm(statsfile)

pyModPandas = PyCall.pyimport_conda("pandas", "pandas")

writer = pyModPandas.ExcelWriter(statsfile, engine="xlsxwriter")
sheet = "sheet1"
# create python dataframe	
dataDict = Dict("x" => data.x,"y"=>data.y)
pyDF = PyCall.pycall(pyModPandas.DataFrame, PyCall.PyObject, dataDict, columns=propertynames(data))
PyCall.pycall(pyDF."to_excel", PyCall.PyAny, writer, header=false, index=false, sheet_name=sheet, startrow=1, startcol=1)
#writer.close()

chartDict = Dict(
"set_y_axis" => Dict("name"=>"Observed Ratio", "major_gridlines"=>Dict{Any, Any}("visible"=>true)),
"set_legend" => Dict("position"=>"right"),
"set_size"   => Dict("y_scale"=>1.5, "x_scale"=>2.4),
"series1"    => Dict("name"=>"=sheet1!\$B\$2", "values"=>"=sheet1!\$B\$3:\$B\$10", "categories"=>"=sheet1!\$C\$3:\$C\$10"),
"add_chart"  => Dict("type"=>"column"),
"set_title"  => Dict("name"=>"Observed Ratio per Segment")
)

workbook = writer.book
	
#crt = Chart("sheet1",chartDict,"D1") 

#writer.close()

chart = PyCall.pycall(workbook."add_chart", PyCall.PyAny, chartDict["add_chart"])

PyCall.pycall(chart."add_series", PyCall.PyAny, chartDict["series1"])

# set other properties
for x in ["set_y_axis","set_legend","set_size","set_title"]
    PyCall.pycall(chart[x], PyCall.PyAny, chartDict[x])	# TBD: unclear how to fix this line... (deprecated syntax...)
end
worksheet = writer.sheets["sheet1"]
PyCall.pycall(worksheet."insert_chart", PyCall.PyAny, "D1", chart)

writer.close()
@show statsfile
return nothing 
end

me(statsfile) 