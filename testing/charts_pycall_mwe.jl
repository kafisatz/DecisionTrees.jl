
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
	for x in keys(chartDict)
		fieldsWhichAreAlreadySet = [convert(String, string("series", zz)) for zz in 1:i]
		resevedKeywords = ["combChart","series_comb1","series_comb2","series_comb3", "series_comb4"] # currently limited to 1+4 series for combined charts
		append!(fieldsWhichAreAlreadySet, ["add_series","add_chart"])
		append!(fieldsWhichAreAlreadySet, resevedKeywords)		
		if !in(x, fieldsWhichAreAlreadySet)
            # TBD: unclear how to fix this line...
            # TBD: unclear how to fix this line... (deprecated syntax...)
            PyCall.pycall(chart[x], PyCall.PyAny, chartDict[x])	 # TBD: unclear how to fix this line...
		end			
	end
	PyCall.pycall(worksheet."insert_chart", PyCall.PyAny, location, chart)
	# writer[:save]()
end