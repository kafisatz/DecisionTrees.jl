#tbd todo add single quotes around sheetnames for all chart formulas
export write_statistics

function create_dataframe(arr::Array{U,2},header::Array{T,1}) where {T <: AbstractString,U <: Any}
	res=convert(DataFrame,arr)
	hsym=Symbol[x for x in header]
	names!(res,hsym)
	return res
end

function create_custom_dict(df::DataFrame)	
	header=names(df)
	d=Dict{AbstractString,Array{Any,1}}()
	for i=1:length(header)			
		d[string(header[i])]=df[i]
	end
	return d
end

function addChartToWorkbook!(workbook::PyObject,worksheet::PyObject,chartDict::Dict{AbstractString,Dict{AbstractString,Any}},location::AbstractString) #;properties=["set_x_axis", "set_y_axis","set_legend"])
	chart = pycall(workbook["add_chart"],PyAny,chartDict["add_chart"])
	stopboolean=true
	i=1
	local thiskey
	while stopboolean
		thiskey=string("series",i)
		if haskey(chartDict,thiskey)
			pycall(chart["add_series"],PyAny,chartDict[thiskey])
		else
			stopboolean=false
			break #if there is no series2 then we assume there is no series 3 to n either
		end
		i+=1
	end
	#check if this is a combined chart
	thiskey="combChart"
	i2=1
	stopboolean=true
	if haskey(chartDict,thiskey)
		second_chart = pycall(workbook["add_chart"],PyAny,chartDict[thiskey])
		while stopboolean
			thiskey=string("series_comb",i2)
			if haskey(chartDict,thiskey)
				pycall(second_chart["add_series"],PyAny,chartDict[thiskey])
			else
				stopboolean=false
				break #if there is no series2 then we assume there is no series 3 to n either
			end
			i2+=1
		end
		#combine the charts
		pycall(chart["combine"],PyAny,second_chart)		
	end
	#set other properties
	for x in keys(chartDict)
		fieldsWhichAreAlreadySet=[convert(String,string("series",zz)) for zz in 1:i]
		resevedKeywords=["combChart","series_comb1","series_comb2","series_comb3", "series_comb4"] #currently limited to 1+4 series for combined charts
		append!(fieldsWhichAreAlreadySet,["add_series","add_chart"])
		append!(fieldsWhichAreAlreadySet,resevedKeywords)		
		if !in(x,fieldsWhichAreAlreadySet)
			pycall(chart[x],PyAny,chartDict[x])	
		end			
	end
	pycall(worksheet["insert_chart"],PyAny,location, chart)
	#writer[:save]()
end

function write_statistics(excelData::ExcelData,statsfile::T,write_header::Bool,write_index::Bool) where {T <: AbstractString} #,charts::Array{Chart,1})		
	#writing an Excel file seems very slow if the file already exists!
	isfile(statsfile)&&rm(statsfile)
	
	writer=writeDFtoExcel(excelData,statsfile,0,0,write_header,write_index)
	workbook = writer[:book]
	#Plot charts	
	for c in excelData.charts
		sheetWhereChartIsLocated=c.sheet		
		worksheet = writer[:sheets][sheetWhereChartIsLocated]
		addChartToWorkbook!(workbook,worksheet,c.chartDict,c.location);
	end
	#save (=write) Excel file and close it	
	writer[:save]()
	return nothing
end

function writeDFtoExcel(excelData::ExcelData,existingFile::T,row::Int,col::Int,write_header::Bool,write_index::Bool) where {T <: AbstractString}
#http://search.cpan.org/~jmcnamara/Excel-Writer-XLSX/lib/Excel/Writer/XLSX.pm
	#@assert isfile(existingFile)
	@assert min(row,col)>=0
	#write all data    
	#writer = pyModPandas.ExcelWriter(existingFile, engine = "xlsxwriter")
    #pyModPandas[:ExcelWriter]
    writer=pyModPandas[:ExcelWriter](existingFile, engine = "xlsxwriter")
    
	for xlSheet in excelData.sheets
		df=xlSheet.data
	    sheet=xlSheet.name
		#create python dataframe	
		dataDict = create_custom_dict(df)
		dataDict = create_custom_dict(df)
		#pyDF = pyModPandas.DataFrame(dataDict,columns=names(df)) #this was working under 0.4 but in 0.5 is converted to a julia dict (instead of being a python DF)
		pyDF=pycall(pyModPandas[:DataFrame], PyObject, dataDict,columns=names(df))		
		pycall(pyDF["to_excel"],PyAny,writer, header=write_header,index=write_index, sheet_name = sheet,startrow=row, startcol=col, encoding="utf-8")  #index=false suppress the rowcount		
	end
	#DataFrame.to_excel(excel_writer, sheet_name='Sheet1', na_rep='', float_format=None, columns=None, header=True, index=True, index_label=None, startrow=0, startcol=0, engine=None, merge_cells=True, encoding=None, inf_rep='inf')
	#writer[:save]()
	return writer
	#book = pyModOpenpyxl.load_workbook(excel_file)
	#snames=pycall(book["get_sheet_names"],PyAny)
	#pyeval("df.to_excel(writer,sheet_name='Sheet1',startrow=5,startcol=3)",df=pyDF,writer=writer)
	#writer[:sheets]=pyeval("dict((ws.title, ws) for ws in book.worksheets)",book=book,writer=writer)
end

function excelLetter(x::Int)
	@assert x>0
	if x>702
		error("Only columns up to 702 are supported")
	else
		d,r=divrem(x,26)
		if r==0
			r+=26
			d-=1
		end		
			if d>0
				letter=string(Char(64+d))
			else
				letter=""			
			end
			letter=string(letter,string(Char(64+r)))
		#end
	end
	return convert(String,letter)
end

function defineRelativityChart(sheetWhereChartIsLocated::T,dataSheet::T,location::T,rows::Int,headercol1::Int,headerrow1::Int;headerrow2::Int=0,headercol2::Int=0,datarow2::Int=0,xtitle="Score Band",ytitle="Relativity",xscale=1.4,title="",datarow::Int=0,valuescol::Int=0,categoriescol::Int=0,valuescol2::Int=0,yscale::Float64=NaN) where {T <: AbstractString}
	chartDict=Dict{AbstractString,Dict{AbstractString,Any}}()
	chartDict["add_chart"]=Dict{Any,Any}("type"=>"column")
	if categoriescol==0
		categoriescol=1
	end
	if datarow==0
		datarow=headerrow1+1
	end
	if valuescol==0
		valuescol=headercol1
	end
	
	#line1
		nameref=string("=",dataSheet,"!\$",excelLetter(headercol1),'$',headerrow1)
		valueref=string("=",dataSheet,"!\$",excelLetter(valuescol),'$',datarow,":\$",excelLetter(valuescol),'\$',datarow+rows-1)
		categoriesref=string("=",dataSheet,"!\$",excelLetter(categoriescol),'$',datarow,":\$",excelLetter(categoriescol),'\$',datarow+rows-1)	
		series=Dict{Any,Any}("name"=>nameref,"categories"=>categoriesref,"values"=>valueref)
		#series1=Dict{Any,Any}("name"=>nameref,"categories"=>categoriesref,"values"=>valueref,"line"=>{"width"=>1.00})
		chartDict["series1"]=deepcopy(series)
	#line2
		if headerrow2==0
			headerrow2=headerrow1
		end		
		if headercol2==0
			headercol2=headercol1+1
		end		
		if datarow2==0
			datarow2=datarow
		end
		if valuescol2==0
			valuescol2=headercol2
		end
		
		nameref=string("=",dataSheet,"!\$",excelLetter(headercol2),'$',headerrow2)
		valueref=string("=",dataSheet,"!\$",excelLetter(valuescol2),'$',datarow2,":\$",excelLetter(valuescol2),'\$',datarow2+rows-1)
		categoriesref=string("=",dataSheet,"!\$",excelLetter(categoriescol),'$',datarow2,":\$",excelLetter(categoriescol),'\$',datarow2+rows-1)	
		series=Dict{Any,Any}("name"=>nameref,"categories"=>categoriesref,"values"=>valueref)
		chartDict["series2"]=deepcopy(series)
	#options
	chartDict["set_x_axis"]=Dict{Any,Any}("name"=>xtitle, "date_axis"=>false)
	chartDict["set_y_axis"]=Dict{Any,Any}("name"=>ytitle, "major_gridlines"=>Dict{Any,Any}("visible"=>true))
	chartDict["set_legend"]=Dict{Any,Any}("position"=>"right")
	if length(title)!=0
		chartDict["set_title"]=Dict{Any,Any}("name"=>title)
	end
	#NOTE: the default chart size is about 14.41 rows in Excel
	if isnan(yscale)
		yscale=1.0/14.41*(min(30,rows)+2)
	end
	chartDict["set_size"]=Dict{Any,Any}("x_scale"=>xscale, "y_scale"=>yscale)
	
	resChart=Chart(convert(String,sheetWhereChartIsLocated),chartDict,convert(String,uppercase(location)))
	return resChart
end

function defineScoreChart(sheetWhereChartIsLocated::T,dataSheet::T,location::T,rows::Int,scoreCol::Int,weightCol::Int,observedCol::Int,SmoothedCol::Int) where {T <: AbstractString}
	#headercol1::Int,headerrow1::Int)
	chartDict=Dict{AbstractString,Dict{AbstractString,Any}}()
	chartDict["add_chart"]=Dict{Any,Any}("type"=>"line")
	categoriescol=scoreCol
	headerrow=1
	datarow=headerrow+1
		
	categoriesref=string("=",dataSheet,"!\$",excelLetter(categoriescol),'$',datarow,":\$",excelLetter(categoriescol),'\$',datarow+rows-1)	
	
	#line1
	valuescol=observedCol
		nameref=string("=",dataSheet,"!\$",excelLetter(valuescol),'$',headerrow)
		valueref=string("=",dataSheet,"!\$",excelLetter(valuescol),'$',datarow,":\$",excelLetter(valuescol),'\$',datarow+rows-1)		
		series=Dict{Any,Any}("name"=>nameref,"categories"=>categoriesref,"values"=>valueref)
		chartDict["series1"]=deepcopy(series)
	#line2
	valuescol=SmoothedCol
		nameref=string("=",dataSheet,"!\$",excelLetter(valuescol),'$',headerrow)
		valueref=string("=",dataSheet,"!\$",excelLetter(valuescol),'$',datarow,":\$",excelLetter(valuescol),'\$',datarow+rows-1)		
		series=Dict{Any,Any}("name"=>nameref,"categories"=>categoriesref,"values"=>valueref)
		chartDict["series2"]=deepcopy(series)
	#line3 on Secondary Axis
	valuescol=weightCol
		nameref=string("=",dataSheet,"!\$",excelLetter(valuescol),'$',headerrow)
		valueref=string("=",dataSheet,"!\$",excelLetter(valuescol),'$',datarow,":\$",excelLetter(valuescol),'\$',datarow+rows-1)				
		series=Dict{Any,Any}("name"=>nameref,"categories"=>categoriesref,"values"=>valueref,"y2_axis"=>1)
		chartDict["series3"]=deepcopy(series)
		
	#options
	chartDict["set_x_axis"]=Dict{Any,Any}("name"=>"Score", "date_axis"=>false)
	chartDict["set_y_axis"]=Dict{Any,Any}("name"=>"", "major_gridlines"=>Dict{Any,Any}("visible"=>true))
	chartDict["set_y2_axis"]=Dict{Any,Any}("name"=>"Weight")
	chartDict["set_legend"]=Dict{Any,Any}("position"=>"bottom")	
	chartDict["set_size"]=Dict{Any,Any}("x_scale"=>3.4, "y_scale"=>2.4)
	chartDict["set_title"]=Dict{Any,Any}("name"=>"Estimates")
	
	resChart=Chart(convert(String,sheetWhereChartIsLocated),chartDict,convert(String,uppercase(location)))
	return resChart
end

function defineUnivariateChart(sheetWhereChartIsLocated::T,dataSheet::T,location::T,charttitle::String,rows::Int,categoriescol::Int,scoreCol::Int,weightCol::Int,headerrow::Int) where {T <: AbstractString}
	#headercol1::Int,headerrow1::Int)
	chartDict=Dict{AbstractString,Dict{AbstractString,Any}}()
	chartDict["add_chart"]=Dict{Any,Any}("type"=>"column")
	#categoriescol=scoreCol
	#headerrow=1
	datarow=headerrow+1		
	categoriesref=string("=",dataSheet,"!\$",excelLetter(categoriescol),'$',datarow,":\$",excelLetter(categoriescol),'\$',datarow+rows-1)	
	
	#line1
	valuescol=weightCol
		nameref=string("=",dataSheet,"!\$",excelLetter(valuescol),'$',headerrow)
		valueref=string("=",dataSheet,"!\$",excelLetter(valuescol),'$',datarow,":\$",excelLetter(valuescol),'\$',datarow+rows-1)		
		series=Dict{Any,Any}("name"=>nameref,"categories"=>categoriesref,"values"=>valueref)	
		chartDict["series1"]=deepcopy(series)
	#line2 on Secondary Axis
	valuescol=scoreCol
		nameref=string("=",dataSheet,"!\$",excelLetter(valuescol),'$',headerrow)
		valueref=string("=",dataSheet,"!\$",excelLetter(valuescol),'$',datarow,":\$",excelLetter(valuescol),'\$',datarow+rows-1)				
		series=Dict{Any,Any}("name"=>nameref,"categories"=>categoriesref,"values"=>valueref,"y2_axis"=>1)
		chartDict["combChart"]=Dict{Any,Any}("type"=>"line")
		chartDict["series_comb1"]=deepcopy(series) #combined Charts (bar/line charts) use a specific notation of combChart instead of Series2
		
	#options
	chartDict["set_x_axis"]=Dict{Any,Any}("name"=>charttitle, "date_axis"=>false)
	chartDict["set_y_axis"]=Dict{Any,Any}("name"=>"Weight", "major_gridlines"=>Dict{Any,Any}("visible"=>true))
	chartDict["set_y2_axis"]=Dict{Any,Any}("name"=>"Average Score")
	chartDict["set_legend"]=Dict{Any,Any}("position"=>"bottom")
	chartDict["set_size"]=Dict{Any,Any}("x_scale"=>3.4, "y_scale"=>1.4)
	chartDict["set_title"]=Dict{Any,Any}("name"=>charttitle)
	
	resChart=Chart(convert(String,sheetWhereChartIsLocated),chartDict,convert(String,uppercase(location)))
	return resChart
end


function defineUnivariateChartWith2Lines(sheetWhereChartIsLocated::T,dataSheet::T,location::T,charttitle::String,rows::Int,categoriescol::Int,line1Col::Int,line2Col::Int,weightCol::Int,headerrow::Int) where {T <: AbstractString}
	#headercol1::Int,headerrow1::Int)
	chartDict=Dict{AbstractString,Dict{AbstractString,Any}}()
	chartDict["add_chart"]=Dict{Any,Any}("type"=>"column")
	#categoriescol=line1Col
	#headerrow=1
	datarow=headerrow+1		
	categoriesref=string("=",dataSheet,"!\$",excelLetter(categoriescol),'$',datarow,":\$",excelLetter(categoriescol),'\$',datarow+rows-1)	
	
	#line1 (which is a bar here)
	valuescol=weightCol
		nameref=string("=",dataSheet,"!\$",excelLetter(valuescol),'$',headerrow)
		valueref=string("=",dataSheet,"!\$",excelLetter(valuescol),'$',datarow,":\$",excelLetter(valuescol),'\$',datarow+rows-1)		
		series=Dict{Any,Any}("name"=>nameref,"categories"=>categoriesref,"values"=>valueref)	
		chartDict["series1"]=deepcopy(series)
	#line2 on Secondary Axis
	valuescol=line1Col
		nameref=string("=",dataSheet,"!\$",excelLetter(valuescol),'$',headerrow)
		valueref=string("=",dataSheet,"!\$",excelLetter(valuescol),'$',datarow,":\$",excelLetter(valuescol),'\$',datarow+rows-1)				
		series=Dict{Any,Any}("name"=>nameref,"categories"=>categoriesref,"values"=>valueref,"y2_axis"=>1)
		chartDict["combChart"]=Dict{Any,Any}("type"=>"line")
		chartDict["series_comb1"]=deepcopy(series) #combined Charts (bar/line charts) use a specific notation of combChart instead of Series2
	#line3 on Secondary Axis
	valuescol=line2Col
		nameref=string("=",dataSheet,"!\$",excelLetter(valuescol),'$',headerrow)
		valueref=string("=",dataSheet,"!\$",excelLetter(valuescol),'$',datarow,":\$",excelLetter(valuescol),'\$',datarow+rows-1)				
		series=Dict{Any,Any}("name"=>nameref,"categories"=>categoriesref,"values"=>valueref,"y2_axis"=>1)
		chartDict["combChart"]=Dict{Any,Any}("type"=>"line")
		chartDict["series_comb2"]=deepcopy(series) #combined Charts (bar/line charts) use a specific notation of combChart instead of Series2
		
	#options
	chartDict["set_x_axis"]=Dict{Any,Any}("name"=>charttitle, "date_axis"=>false)
	chartDict["set_y_axis"]=Dict{Any,Any}("name"=>"Weight", "major_gridlines"=>Dict{Any,Any}("visible"=>true))
	chartDict["set_y2_axis"]=Dict{Any,Any}("name"=>"Average Score")
	chartDict["set_legend"]=Dict{Any,Any}("position"=>"bottom")
	chartDict["set_size"]=Dict{Any,Any}("x_scale"=>3.4, "y_scale"=>1.4)
	chartDict["set_title"]=Dict{Any,Any}("name"=>charttitle)
	
	resChart=Chart(convert(String,sheetWhereChartIsLocated),chartDict,convert(String,uppercase(location)))
	return resChart
end

function defineTwoWayCharts(sheetWhereChartIsLocated::T,dataSheet::T,location::T,charttype::T,firstdatarow::Int,rows::Int,headercol::Int,nClasses::Int,descriptorcolumn::Int,xaxisname::T,yaxisname::T,
		title::T;legendpos="top",categoriescol=1,xaxisformat="",yaxisformat="",xscale=3,yscale=1.5) where {T <: AbstractString}
	#http://xlsxwriter.readthedocs.org/en/latest/chart.html
	chartDict=Dict{AbstractString,Dict{AbstractString,Any}}()
	chartDict["add_chart"]=Dict{Any,Any}("type"=>charttype)
	datarow=firstdatarow
	categoriesref=string("=",dataSheet,"!\$",excelLetter(categoriescol),'$',datarow,":\$",excelLetter(categoriescol),'\$',datarow+rows-1)		

	valuescol=headercol
	for count=1:nClasses # in headercols						
		chartDict[string("series",count)]=deepcopy(Dict{Any,Any}("name"=>string("=",dataSheet,"!\$",excelLetter(descriptorcolumn),'$',datarow),"categories"=>categoriesref,
				"values"=>string("=",dataSheet,"!\$",excelLetter(valuescol),'$',datarow,":\$",excelLetter(valuescol),'\$',datarow+rows-1)))
		datarow+=rows+2
	end
	#options
	chartDict["set_x_axis"]=Dict{Any,Any}("name"=>xaxisname, "date_axis"=>false)
	chartDict["set_y_axis"]=Dict{Any,Any}("name"=>yaxisname,"major_gridlines"=>Dict{Any,Any}("visible"=>true))
	chartDict["set_legend"]=Dict{Any,Any}("position"=>legendpos)
	chartDict["set_size"]=Dict{Any,Any}("x_scale"=>xscale, "y_scale"=>yscale)
	chartDict["set_title"]=Dict{Any,Any}("name"=>title)
	if length(xaxisformat)>0
		chartDict["set_x_axis"]["num_format"]=xaxisformat
	end
	if length(yaxisformat)>0
		chartDict["set_y_axis"]["num_format"]=yaxisformat
	end
	resChart=Chart(convert(String,sheetWhereChartIsLocated),chartDict,convert(String,uppercase(location)))
	return resChart
end


function defineChartWithNSeries(sheetWhereChartIsLocated::T,dataSheet::T,location::T,charttype::T,firstdatarow::Int,rows::Int,headercols::Array{Int,1},headerrow1::Int,xaxisname::T,yaxisname::T,title::T;legendpos="top",categoriescol=1,xaxisformat="",yaxisformat="",xscale=3,yscale=1.5) where {T <: AbstractString}
	#http://xlsxwriter.readthedocs.org/en/latest/chart.html
	chartDict=Dict{AbstractString,Dict{AbstractString,Any}}()
	chartDict["add_chart"]=Dict{Any,Any}("type"=>charttype)
	datarow=firstdatarow
	categoriesref=string("=",dataSheet,"!\$",excelLetter(categoriescol),'$',datarow,":\$",excelLetter(categoriescol),'\$',datarow+rows-1)		

	count=0
	for ii in headercols
		count+=1	
		valuescol=ii;headercol=ii #line1
		chartDict[string("series",count)]=deepcopy(Dict{Any,Any}("name"=>string("=",dataSheet,"!\$",excelLetter(headercol),'$',headerrow1),"categories"=>categoriesref,
				"values"=>string("=",dataSheet,"!\$",excelLetter(valuescol),'$',datarow,":\$",excelLetter(valuescol),'\$',datarow+rows-1)))
	end
	#options
	chartDict["set_x_axis"]=Dict{Any,Any}("name"=>xaxisname, "date_axis"=>false)
	chartDict["set_y_axis"]=Dict{Any,Any}("name"=>yaxisname,"major_gridlines"=>Dict{Any,Any}("visible"=>true))
	chartDict["set_legend"]=Dict{Any,Any}("position"=>legendpos)
	chartDict["set_size"]=Dict{Any,Any}("x_scale"=>xscale, "y_scale"=>yscale)
	chartDict["set_title"]=Dict{Any,Any}("name"=>title)
	if length(xaxisformat)>0
		chartDict["set_x_axis"]["num_format"]=xaxisformat
	end
	if length(yaxisformat)>0
		chartDict["set_y_axis"]["num_format"]=yaxisformat
	end
	resChart=Chart(convert(String,sheetWhereChartIsLocated),chartDict,convert(String,uppercase(location)))
	return resChart
end

function defineChartWith2Series(sheetWhereChartIsLocated::T,dataSheet::T,location::T,charttype::T,rows::Int,headercol1::Int,headercol2::Int,headerrow1::Int) where {T <: AbstractString}
	#http://xlsxwriter.readthedocs.org/en/latest/chart.html
	chartDict=Dict{AbstractString,Dict{AbstractString,Any}}()
	chartDict["add_chart"]=Dict{Any,Any}("type"=>charttype)
	categoriescol=1
	datarow=headerrow1+1
	
	categoriesref=string("=",dataSheet,"!\$",excelLetter(categoriescol),'$',datarow,":\$",excelLetter(categoriescol),'\$',datarow+rows-1)	
	
	#line1
	valuescol=headercol1
	headercol=headercol1
		nameref=string("=",dataSheet,"!\$",excelLetter(headercol),'$',headerrow1)
		valueref=string("=",dataSheet,"!\$",excelLetter(valuescol),'$',datarow,":\$",excelLetter(valuescol),'\$',datarow+rows-1)		
		series=Dict{Any,Any}("name"=>nameref,"categories"=>categoriesref,"values"=>valueref)
		chartDict["series1"]=deepcopy(series)
	#line2
		headercol=headercol2
		valuescol=headercol2
		nameref=string("=",dataSheet,"!\$",excelLetter(headercol),'$',headerrow1)
		valueref=string("=",dataSheet,"!\$",excelLetter(valuescol),'$',datarow,":\$",excelLetter(valuescol),'\$',datarow+rows-1)		
		series=Dict{Any,Any}("name"=>nameref,"categories"=>categoriesref,"values"=>valueref)
		chartDict["series2"]=deepcopy(series)
	#options
	chartDict["set_x_axis"]=Dict{Any,Any}("name"=>"Iteration", "date_axis"=>false)
	chartDict["set_y_axis"]=Dict{Any,Any}("major_gridlines"=>Dict{Any,Any}("visible"=>true))
	chartDict["set_legend"]=Dict{Any,Any}("position"=>"top")
	chartDict["set_size"]=Dict{Any,Any}("x_scale"=>3, "y_scale"=>1.5)

	resChart=Chart(convert(String,sheetWhereChartIsLocated),chartDict,convert(String,uppercase(location)))
	return resChart
end


function add_iteration_charts!(xlData::ExcelData,settings::ModelSettings,startrow)
	chartRowStep=22
	
	currentChartRow=startrow
	thischart=defineChartWith2Series("Overview","ModelStatistics",string("A",currentChartRow),"line",settings.niter,2,5,1)
	push!(xlData.charts,deepcopy(thischart))
	currentChartRow+=chartRowStep
	thischart=defineChartWith2Series("Overview","ModelStatistics",string("A",currentChartRow),"line",settings.niter,6,7,1)
	push!(xlData.charts,deepcopy(thischart))
	currentChartRow+=chartRowStep
	thischart=defineChartWith2Series("Overview","ModelStatistics",string("A",currentChartRow),"line",settings.niter,3,4,1)
	push!(xlData.charts,deepcopy(thischart))
	currentChartRow+=chartRowStep
	thischart=defineChartWith2Series("Overview","ModelStatistics",string("A",currentChartRow),"column",settings.niter,8,9,1)
	push!(xlData.charts,deepcopy(thischart))
	currentChartRow+=chartRowStep
	thischart=defineChartWithNSeries("Overview","ModelStatistics",string("A",currentChartRow),"line",2,settings.niter,Int[14,15,18,19],1,"Iteration","","")		
	push!(xlData.charts,deepcopy(thischart))
	currentChartRow+=chartRowStep
	thischart=defineChartWithNSeries("Overview","ModelStatistics",string("A",currentChartRow),"line",2,settings.niter,Int[12,13,16,17],1,"Iteration","","")		
	push!(xlData.charts,deepcopy(thischart))
	currentChartRow+=chartRowStep
	thischart=defineChartWith2Series("Overview","ModelStatistics",string("A",currentChartRow),"line",settings.niter,22,23,1)
	push!(xlData.charts,deepcopy(thischart))
	currentChartRow+=chartRowStep
	thischart=defineChartWith2Series("Overview","ModelStatistics",string("A",currentChartRow),"line",settings.niter,20,21,1)
	push!(xlData.charts,deepcopy(thischart))
	#Feb 2018 
	#add charts for RSquared of Numerator
	currentChartRow+=chartRowStep
	thischart=defineChartWith2Series("Overview","ModelStatistics",string("A",currentChartRow),"line",settings.niter,26,27,1)
	push!(xlData.charts,deepcopy(thischart))	
	#add charts for RSS of numerator
	currentChartRow+=chartRowStep
	thischart=defineChartWith2Series("Overview","ModelStatistics",string("A",currentChartRow),"line",settings.niter,30,31,1)
	push!(xlData.charts,deepcopy(thischart))		
return nothing
end
