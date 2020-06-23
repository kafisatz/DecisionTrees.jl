using DataFrames 
using PyCall

mutable struct Chart	
	sheet::String
	chartDict::Dict{AbstractString,Dict{AbstractString,Any}}
	location::String
end

mutable struct ExcelSheet
	name::String #I think this name is currently not used anywhere
	data::DataFrame
end

mutable struct ExcelData
	sheets::Array{ExcelSheet,1}
    charts::Array{Chart,1}
    function ExcelData()
        return new(Array{ExcelSheet}(undef,0),Array{Chart}(undef,0))
    end
    function ExcelData(a,b)
        return  new(a,b)
    end
end

cdict=Dict{AbstractString,Dict{AbstractString,Any}}(
  "set_y_axis" => Dict{AbstractString,Any}("name"=>"Observed Ratio","major_gridlines"=>Dict{Any,Any}("visible"=>true)),
  "set_legend" => Dict{AbstractString,Any}("position"=>"right"),
  "set_size"   => Dict{AbstractString,Any}("y_scale"=>1.5,"x_scale"=>2.4),
  "series1"    => Dict{AbstractString,Any}("name"=>"=ModelStatistics!\$A\$1","values"=>"=ModelStatistics!\$F\$3:\$F\$6","categories"=>"=ModelStatistics!\$A\$3:\$A\$6"),
  "add_chart"  => Dict{AbstractString,Any}("type"=>"column")
)

aChart=Chart("sheet1",cdict,"K2")

xlData=ExcelData()
sheet1df=DataFrame(rand(10,10)))
xlSheet=ExcelSheet("sheet1",sheet1df)

push!(xlData.sheets,sheet1df)
push!(xlData.charts,aChart)

ct.sheet
julia> resm.exceldata.charts[1].sheet


function writeStatistics(excelData::ExcelData, statsfile::T, write_header::Bool, write_index::Bool) where {T <: AbstractString} # ,charts::Array{Chart,1})		
	# writing an Excel file seems very slow if the file already exists!
	isfile(statsfile) && rm(statsfile)
	
	writer = writeDFtoExcel(excelData, statsfile, 0, 0, write_header, write_index)
	workbook = writer.book
	# Plot charts	
	for c in excelData.charts
		sheetWhereChartIsLocated = c.sheet		
		worksheet = writer.sheets[sheetWhereChartIsLocated]
		addChartToWorkbook!(workbook, worksheet, c.chartDict, c.location);
	end
	# save (=write) Excel file and close it	
	writer.save()
	return nothing
end



function addChartToWorkbook!(workbook::PyCall.PyObject, worksheet::PyCall.PyObject, chartDict::Dict{AbstractString,Dict{AbstractString,Any}}, location::AbstractString) # ;properties=["set_x_axis", "set_y_axis","set_legend"])
	chart = PyCall.pycall(workbook."add_chart", PyCall.PyAny, chartDict["add_chart"])
	stopboolean = true
	i = 1
	local thiskey
	while stopboolean
		thiskey = string("series", i)
		if haskey(chartDict, thiskey)
			PyCall.pycall(chart."add_series", PyCall.PyAny, chartDict[thiskey])
		else
			stopboolean = false
			break # if there is no series2 then we assume there is no series 3 to n either
		end
		i += 1
	end
	# check if this is a combined chart
	thiskey = "combChart"
	i2 = 1
	stopboolean = true
	if haskey(chartDict, thiskey)
		second_chart = PyCall.pycall(workbook."add_chart", PyCall.PyAny, chartDict[thiskey])
		while stopboolean
			thiskey = string("series_comb", i2)
			if haskey(chartDict, thiskey)
				PyCall.pycall(second_chart."add_series", PyCall.PyAny, chartDict[thiskey])
			else
				stopboolean = false
				break # if there is no series2 then we assume there is no series 3 to n either
			end
			i2 += 1
		end
		# combine the charts
		PyCall.pycall(chart."combine", PyCall.PyAny, second_chart)		
	end
    # set other properties
    
    xlData = ExcelData(Array{ExcelSheet}(undef, 0), Array{Chart}(undef, 0))
    
    #chartDict::Dict{AbstractString,Dict{AbstractString,Any}}
    chart = PyCall.pycall(workbook."add_chart", PyCall.PyAny, chartDict["add_chart"])

	for x in keys(chartDict)
		fieldsWhichAreAlreadySet = [convert(String, string("series", zz)) for zz in 1:i]
		resevedKeywords = ["combChart","series_comb1","series_comb2","series_comb3", "series_comb4"] # currently limited to 1+4 series for combined charts
		somefields=["add_series","add_chart"]		
		if !in(x, somefields)
            # TBD: unclear how to fix this line...
            # TBD: unclear how to fix this line... (deprecated syntax...)
            PyCall.pycall(chart[x], PyCall.PyAny, chartDict[x])	 # TBD: unclear how to fix this line...
		end			
	end
	PyCall.pycall(worksheet."insert_chart", PyCall.PyAny, location, chart)
	# writer[:save]()
end