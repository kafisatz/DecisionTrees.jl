export resample_trnvalidx!,removeBOM,sampleData,subset_pda_mod #for debugging purposes only

import Base: append!,mean,length,findin,isless,eltype,resize!,convert
#import PooledArrays.levels
#import MySQL: mysql_query,mysql_display_error,mysql_store_result,mysql_affected_rows,mysql_free_result,mysql_num_fields,mysql_fetch_fields,mysql_num_rows,mysql_get_julia_type,mysql_fetch_row,MYSQL_RES,MYSQL_TYPE,MYSQL_ROW #temporarily disabled
export createZipFile

#this could certainly be done nice
#tries to find the 'max' type of the refs (which is likely UInt8 or UInt16)
function find_max_type(x::DataFrame)
    for col in 1:size(x,2)
        if eltype(x[col].refs)==UInt32 
            return UInt32
        end
    end
    for col in 1:size(x,2)
        if eltype(x[col].refs)==UInt16
            return UInt16
        end
    end
    return UInt8
end

#it is critical to define levels properly for PooledArray
function Missings.levels(x::PooledArray)
    return deepcopy(x.pool)
end

function resample_trnvalidx!(x::DTMTable,trnprop::Float64)
	@assert trnprop>0
	@assert trnprop<1
	sz=length(x.weight)
	trnsz=round(Int,trnprop*sz)
	trnsz=max(1,min(sz,trnsz))
	t=sample(1:sz,trnsz,replace=false,ordered=true)
	v=setdiff(1:sz,t)
	@assert issorted(t)
	@assert issorted(v)
	x.trnidx=t
	x.validx=v
	   
	return nothing
end

function get_stats(model::Tree)
	#typeof(b)==Tree
	#fieldnames(b)
	statsdf=model.exceldata.sheets[2].data
	#find the last rows of the results sheet
	str="Correlation"
	rows=findin(statsdf[:,1],[str])
	if length(rows)!=1
		@show rows
		error("Multiple rows with 'name' Correlation found. This is not expected")
	end
	stats=statsdf[rows[1]:end,1:2]	
	desc=string.(stats[:,1])
	numbers=float.(stats[:,2])
	setti=model.exceldata.sheets[1].data[:,2]
	setti_desc=string.(model.exceldata.sheets[1].data[:,1])
	return desc,numbers,setti_desc,setti
end

function get_stats(model::Union{BoostedTree,BaggedTree})
	numbers=Array{Float64,2}(model.modelstats[end,:])[:]
	desc=names(model.modelstats)
	settdf=convert(DataFrame,writeAllFieldsToArray(model.settings))
	setti_desc=string.(settdf[:,1])
	setti=settdf[:,2]
	return desc,numbers,setti_desc,setti
end

function try_convert(T,x)
	local res::T
	try
		res=convert(T,x)
	catch e
		#@show e
		#dump(e)
		@assert false "Unable to convert x" # to Integer array in function quadratic_weighted_kappa"
	end
	return res
end

function write_estimatedNumerator_to_mat!(mat::Array{Float64,2},col::Int,idx::Vector{Int},ref::Vector{Float64},scores::Vector{Int})
	@inbounds for z in idx
		mat[z,col]=ref[scores[z]]
	end				
return nothing
end

function map_these_values(rawObservedRatioPerScore::Vector{Float64},scores::Vector{Int})
	#estimateUnsmoothed_x_some_test=map_these_values(rawObservedRatioPerScore,scores)
	vec=zeros(Float64,length(scores))
	for i=1:length(scores)
		@inbounds sc=scores[i]
		@inbounds vec[i]=rawObservedRatioPerScore[sc]
	end
	return vec 
end

function update_mapped_estFromScores!(estFromScores::Vector{Float64},referenceForNumEstimates::Vector{Float64},scores::Vector{Int})
	for i=1:length(scores)
		@inbounds sc=scores[i]
		@inbounds estFromScores[i]=referenceForNumEstimates[sc]
	end
	return nothing
end

function write_estimatedNumerator_to_mat!(mat::Array{Float64,2},col::Int,ref::Vector{Float64},scores::Vector{Int})
	@inbounds for z=1:size(mat,1)
		mat[z,col]=ref[scores[z]]
	end				
return nothing
end

function fill_some_elements!(vec::Vector{T},idx::Vector{Int},x::T) where T<:Number
	for i in idx
		@inbounds vec[i]=x
	end
	return nothing
end

function write_column!(idx,mat::Array{T,2},col::Int,vec::Vector{T}) where T<:Number
    for i in idx 
        @inbounds mat[i,col]=vec[i]
    end
    return nothing
end

function write_column!(mat::Array{T,2},col::Int,vec::Vector{T}) where T<:Number
    for i=1:size(mat,1) 
        @inbounds mat[i,col]=vec[i]
    end
    return nothing
end

function update_current_rels!(idx,vec,estratios,value::Float64)
	for i in idx 
		@inbounds vec[i]=estratios[i]/value
	end				
	return nothing
end

function update_current_rels!(vec,estratios,value::Float64)
	for i=1:length(vec)
		@inbounds vec[i]=estratios[i]/value
	end				
	return nothing
end

function update_current_rels!(idx,vec,estratios,mat::Array{Float64,1})
	for i in idx 
		@inbounds vec[i]=estratios[i]/mat[i,1]
	end				
	return nothing
end

function update_current_rels!(vec,estratios,mat::Array{Float64,1})
	for i=1:length(vec)
		@inbounds vec[i]=estratios[i]/mat[i,1]
	end				
	return nothing
end

function quadratic_weighted_kappa(actual::Array{T,1},estimate::Array{U,1}) where {T <: Number,U <: Number}
	estimate_int=try_convert(Array{Int,1},estimate)
	actual_int=try_convert(Array{Int,1},actual)
	return quadratic_weighted_kappa(actual_int,estimate_int)
end

function frequency_count(actual::Array{Int,1},estimate::Array{Int,1})
(lo, hi) = extrema(actual)
	(lo_e,hi_e)=extrema(estimate)
	#out of range estimates are not valid and should produce an error
	@assert lo<=lo_e
	@assert hi_e<=hi
	#NOTE: what do we do if the observed range of actual is limited, e.g. if there is no actual with value 1? shall we start at 2?
	#I assume that this is a rare case
	if lo!=1
		warn("Lo of actual is not equal to 1: lo=$(lo), hi=$(hi)")
	end
	ooo=one(lo)-lo
	vecsize=hi+ooo
	cnt = zeros(Int, vecsize,vecsize)
	hist_a = zeros(Int, vecsize)
	hist_e = zeros(Int, vecsize)
	for count=1:length(actual)
		@inbounds idxrow=actual[count] + ooo
		@inbounds idxcol=estimate[count] + ooo
		@inbounds cnt[idxrow,idxcol] += 1
		@inbounds hist_a[idxrow]+=1
		@inbounds hist_e[idxcol]+=1
	end
	return hist_a,hist_e,vecsize,hi,lo,hi_e,lo_e,cnt
end

function quadratic_weighted_kappa(actual::Array{Int,1},estimate::Array{Int,1})
	#this is temporarily disabled
	return -999.999

	@assert size(actual)==size(estimate)
	hist_a,hist_e,vecsize,hi,lo,hi_e,lo_e,o=frequency_count(actual,estimate)
	w=zeros(Float64,vecsize,vecsize)
	#nq=Float64((vecsize-1)^2)
	for i=1:vecsize,j=1:vecsize
		w[i,j]=((i-j)^2) #/nq
	end
	len=vecsize #size(actual,1)
	hist_a=hist_a./len
	hist_e=hist_e./len
	e=hist_a*reshape(hist_e,1,length(hist_e))
	e=e./sum(e)
	o=o./sum(o)
	kappa=1-sum(w.*o)/sum(w.*e)
	#hack return -kappa since some post processing methods in DTM minimize the error
	kappa=(-1.0)*kappa
	return kappa
end

#function tryQWKappa{U<:Number,T<:Number}(actual::Array{U,1},estimate::Array{T,1})
function tryQWKappa(actual,estimate)
	return NaN #currently disabled
	local estimate_int,actual_int
	try
		actual_int=try_convert(Array{Int,1},actual)
		estimate2=round.(estimate)
		estimate_int=try_convert(Array{Int,1},estimate2)
	catch
		return NaN
	end
	mi,ma=extrema(estimate_int)
	limit_estimate!(estimate_int,mi,ma)
	return quadratic_weighted_kappa(actual_int,estimate_int)
end

function limit_estimate!(est::Array{T,1},mi::Number,ma::Number) where {T <: Number}
	@assert mi<=ma
	for i=1:length(est)
		est[i]=min(max(est[i],mi),ma)
	end
end

function select_metric_col(allmodelStats::DataFrame,metric::AbstractString)
	@assert in(Symbol(metric),names(allmodelStats)) "DataFrame has no column named $(metric)"
	return allmodelStats[Symbol(metric)]
end

function best_models(allmodelStats::DataFrame,metric)
	metric_col=select_metric_col(allmodelStats,metric)
	res2=Array{Any}(1+length(unique(allmodelStats[1])),3)
		res2[1,1]="Model Number"
		res2[1,2]="Minimal Error"
		res2[1,3]=string(metric," at Final Iteration")
	itercol=convert(Array{Int,1},allmodelStats[2])
	#sort data
		srt=sortperm(metric_col)
		column=allmodelStats[1] #hard coded we could also select the column ModelNumber...
		column_srt=column[srt]
		metric_col_srt=metric_col[srt]
		itercol=itercol[srt]
	uq=unique(column_srt)
	k=2
	#keep only one entry (the best one) per model
	for i in uq
		res2[k,1]=i
		idx=findlast(column,i)
		res2[k,3]=metric_col[idx]

		idx2=findfirst(column_srt,i)
		res2[k,2]=metric_col_srt[idx2]
		k+=1
	end
	return res2
end

function best_model_per_iter(allmodelStats,metric)
#consider the best model at each iteration (the model which has lowest error at the final iteration is not necessarily the best model, if we stop already after 200 iterations!)
	metric_col=select_metric_col(allmodelStats,metric)
	narows=isna(metric_col)
	itercol=convert(Array{Int,1},allmodelStats[2])
	modelnrcol=allmodelStats[1]
	narows2=find(narows)
	sort!(narows2)
	if size(narows2,1)>0
		deleteat!(itercol,narows2)
		deleteat!(metric_col,narows2)
		deleteat!(modelnrcol,narows2)
	end
	itrs=unique(itercol)
	sort!(itrs)
	res=Array{Any}(1+length(itrs),3)
	res[1,1]="Iteration"
	res[1,2]="Best Model"
	res[1,3]=metric
	for i in itrs
		val,mdl=select_best_model(itercol,metric_col,modelnrcol,i)
		res[i+1,1]=i
		res[i+1,2]=mdl
		res[i+1,3]=val
	end
	return res
end

function select_best_model(itercol,metric_col,modelnrcol,i)
	idx=itercol.==i
	v=metric_col[idx]
	m=modelnrcol[idx]
	val,loc=findmin(v)
	return val,m[loc]
end

function eta_before(t0::DateTime,model_i_which_starts_now::Int,nmodels::Int)
	#here the model with number model_i_which_starts_now is just starting (thus it is not yet completed)
	elapsed=(now()-t0)
	elapsed_s=elapsed.value/1000
	elapsed_m=elapsed_s/60
	eta=elapsed_s/max(1,model_i_which_starts_now-1)*(nmodels-model_i_which_starts_now-1)
	return eta
end

function eta_after(t0::DateTime,model_i_which_is_completed::Int,nmodels::Int)
	#here the model with number model_i_which_is_completed has just starting completed
	elapsed=(now()-t0)
	elapsed_s=elapsed.value/1000
	elapsed_m=elapsed_s/60
	eta=elapsed_s/model_i_which_is_completed*(nmodels-model_i_which_is_completed)
	return eta
end

function readSettings(settingsFilename)
	local settingsArray,number_of_num_features
	try
		println("Importing Settings CSV: \r\n $(settingsFilename)")
		settingsArray=readcsv(settingsFilename,String) #it is read as String because otherwise some values are parsed as floats (which I want to suppress here)
	catch readerror
		@show readerror
		error("Unable to read file $(settingsFilename)")
	end
	try
		settingsArray=removeBOM(settingsArray)
		coln=findin([lowercase(x) for x in settingsArray[1,:]],["number_of_num_features"])
		@assert size(coln)==(1,) global_number_of_num_f_warning_mandatory_field
		@assert size(settingsArray,1)>1 "SettingsArray must have more than 1 row. Abort"
		coln2=coln[1]
		number_of_num_features=parse(Int,settingsArray[2,coln2])
	catch readerror
		@show readerror
		error("Unable to remove BOM from file $(settingsFilename)")
	end
	@assert size(settingsArray,1)>=2 "Error: Settings file: $(settingsFilename) has less than two row. Abort."
	return settingsArray,number_of_num_features
end

function convert(t::Type{Array},setti::ModelSettings)
	fn=fieldnames(setti)
	res=Array{Any}(2,length(fn))
	i=0
	for f in fn
		i+=1
		res[1,i]=string(f)
		res[2,i]=string(getfield(setti,f))
	end
	return res
end

function queryWithTypes(db,sql::AbstractString, values=[])
#this is a modification of query from SQLite.jl
#the function infers the type from the SQLite DB which is (intentionally) not done by the author of the package.
    stmt = SQLiteStmt(db,sql)
    bind(stmt, values)
    status = execute(stmt)
    ncols = SQLite.sqlite3_column_count(stmt.handle)
    if status == SQLite.SQLITE_DONE || ncols == 0
        return changes(db)
    end
    colnames = Array{AbstractString}(ncols)
    results = Array{Any}(ncols)
    for i = 1:ncols
        colnames[i] = bytestring(SQLite.sqlite3_column_name(stmt.handle,i-1))
        #results[i] = Any[]
        #if t == SQLITE_INTEGER
        #    results[i] = Int64[]
        #elseif t == SQLITE_FLOAT
		t = SQLite.sqlite3_column_type(stmt.handle,i-1)
		if t == SQLite.SQLITE_FLOAT
            results[i] = Float64[]
        elseif t == SQLite.SQLITE_TEXT
            results[i] = String[]
        else
            results[i] = Any[]
        end
    end
    while status == SQLite.SQLITE_ROW
        for i = 1:ncols
            t = SQLite.sqlite3_column_type(stmt.handle,i-1)
            if t == SQLite.SQLITE_INTEGER
                r = SQLite.sqlite3_column_int64(stmt.handle,i-1)
            elseif t == SQLite.SQLITE_FLOAT
                r = SQLite.sqlite3_column_double(stmt.handle,i-1)
            elseif t == SQLite.SQLITE_TEXT
                #TODO: have a way to return text16?
                r = bytestring(SQLite.sqlite3_column_text(stmt.handle,i-1) )
            elseif t == SQLITE_BLOB
                blob = SQLite.sqlite3_column_blob(stmt.handle,i-1)
                b = SQLite.sqlite3_column_bytes(stmt.handle,i-1)
                buf = zeros(UInt8,b)
                unsafe_copy!(pointer(buf), convert(Ptr{UInt8},blob), b)
                r = SQLite.sqldeserialize(buf)
            else
                r = NULL
            end
            push!(results[i],r)
        end
        status = SQLite.sqlite3_step(stmt.handle)
    end
    if status == SQLite.SQLITE_DONE
        return SQLite.ResultSet(colnames, results)
    else
        SQLite.sqliteerror(stmt.db)
    end
end

function write_any_array_to_excel_file(arr::Array,file::String)
	ascii(file)
	isfile(file)&&rm(file)
	sheets=[ExcelSheet("sheet1",convert(DataFrame,arr))]
	xld=ExcelData(sheets,Array{Chart}(0))
	write_statistics(xld,file,false,false)
	return nothing
end

function write_any_array_to_excel_file(arr::Array,arr2::Array,file::String)
	ascii(file)
	isfile(file)&&rm(file)
	sheets=[ExcelSheet("sheet1",convert(DataFrame,arr)),ExcelSheet("sheet2",convert(DataFrame,arr2))]
	xld=ExcelData(sheets,Array{Chart}(0))
	write_statistics(xld,file,false,false)
	return nothing
end

function write_any_array_to_excel_file(arr::Array,arr2::Array,arr3::Array,file::String)
	ascii(file)
	isfile(file)&&rm(file)
	sheets=[ExcelSheet("sheet1",convert(DataFrame,arr)),ExcelSheet("sheet2",convert(DataFrame,arr2)),ExcelSheet("sheet3",convert(DataFrame,arr3))]
	xld=ExcelData(sheets,Array{Chart}(0))
	write_statistics(xld,file,false,false)
	return nothing
end

function write_any_array_to_excel_file(arr::Array,arr2::Array,arr3::Array,arr4::Array,file::String)
	ascii(file)
	isfile(file)&&rm(file)
	sheets=[ExcelSheet("sheet1",convert(DataFrame,arr)),ExcelSheet("sheet2",convert(DataFrame,arr2)),ExcelSheet("sheet3",convert(DataFrame,arr3)),ExcelSheet("sheet4",convert(DataFrame,arr4))]
	xld=ExcelData(sheets,Array{Chart}(0))
	write_statistics(xld,file,false,false)
	return nothing
end

function write_any_array_to_excel_file(arr::Array{Array{T,2},1},file::String) where {T <: Any}
	ascii(file)
	isfile(file)&&rm(file)
	sheets=[ExcelSheet(string("sheet",i),convert(DataFrame,arr[i])) for i=1:length(arr)]
	xld=ExcelData(sheets,Array{Chart}(0))
	write_statistics(xld,file,false,false)
	return nothing
end

function gini_single_argument(x)
    n = length(x)
    xx = sort(x)
    2*(sum(collect(1:n).*xx))/(n*sum(xx))-1
end

function gini(actual,estimate,weight)
	#weight=ones(eltype(estimate),length(estimate))
	#@assert (size(actual)==size(estimate)==size(weight))&&(size(actual,2)==1) # I think it takes too much time to check this as we call this fn very often
	arrlen=length(actual)
	srt=zeros(Int,arrlen) #pre allocation seems to increase performance here (not sure why though....)
	sortperm!(srt,estimate,rev=true)
	#srt=sortperm(estimate,rev=true)
	sumactual=0.0
	@inbounds for i=1:arrlen
		sumactual+=actual[i]*weight[i]
	end
	cums=0.0
	res=0.0
	#the summation is not quite accurate here (accuracy is at e-7 or so) (the reason is the fact that we do not sort and sum the data)
	@inbounds for i in srt #1:arrlen
		cums+=(actual[i]*weight[i])
		res+=cums/sumactual-i/arrlen
	end
	return res
end

function gini_more_accurate_but_slower(a,estimate::Array{T,1},weight=ones(T,length(estimate))) where {T <: Number}
	#@assert (size(actual)==size(estimate)==size(weight))&&(size(actual,2)==1) # I think it takes too much time to check this as we call this fn very often
	actual=deepcopy(a)
	arrlen=length(actual)
	srt=zeros(Int,arrlen) #pre allocation seems to incrase performance here (not sure why though....)
	sortperm!(srt,estimate,rev=true)
	#srt=sortperm(estimate,rev=true)
	for i=1:arrlen
		@inbounds actual[i]=actual[i]*weight[i]
	end
	actual_sorted=actual[srt]
	incrementalToCumulative!(actual_sorted)
	#actual_sorted_float=convert(Array{Float64,1},actual_sorted)

	sumactual=float(sum(actual))
	res=0.0
	for i=1:arrlen
		@inbounds res+=actual_sorted[i]/sumactual-float(i/arrlen)
	end
	return res
end

function normalized_gini(actual,estimate::SubArray{T,1},weight=ones(T,length(estimate))) where {T <: Number}
	return gini(actual,estimate,weight)/gini(actual,actual,weight)
end

function position_in_sorted_array(arr,x)
	#arr needs to be sorted for this to work
	#returns the index i at which arr equals to x
	#returns length(arr)+1 if x does not ocurr within arr
	pos=searchsortedfirst(arr,x)
	arr_length_plus1=length(arr)+1
	if (pos<arr_length_plus1)&&(arr[pos]==x)
		return pos
	else
		return arr_length_plus1
	end
end

function merged_fitted_vectors(obs::Int,sampleVector,fitted_values_on_sampled_data,unusedSamplePart,ooBagEstimates)
#note, we could use inbounds below in the loop provided that we first check
# all_indices=append!(unique(sampleVector),unique(unusedSamplePart))
# countsort!(all_indices)
# @assert all_indices==collect(1:obs)
	#obs=length(denominator)
	fitted=Array{Float64}(obs)
	srt=sortperm(sampleVector)
	sortedSampleVectorWhichMayContainDuplicates=sampleVector[srt]
	#maxIndexOfSample=sortedSampleVectorWhichMayContainDuplicates[end]
	fitted_values_on_sampled_data_sorted=fitted_values_on_sampled_data[srt]	
	len_sample_vector=length(sortedSampleVectorWhichMayContainDuplicates)	
	for z=1:obs		
		pos=position_in_sorted_array(sortedSampleVectorWhichMayContainDuplicates,z)		
		if pos<=len_sample_vector
			fitted[z]=fitted_values_on_sampled_data_sorted[pos]			
		else			
			#it is important that unusedSamplePart is SORTED (which should hold by definition!)
			pos=position_in_sorted_array(unusedSamplePart,z)
			fitted[z]=ooBagEstimates[pos]			
		end
	end
	return fitted
end

function print_some_settings(s::ModelSettings,fields::Array{String,1})
		str=String("")
		for x in fields
			str=string(str,string(", ",x,": ",getfield(s,Symbol(x))))
		end
		println(str)
end
#=
temporarily disabled
		function mypopulate_row!(arr::Array, mysqlfield_types::Array{MYSQL_TYPE}, result::MYSQL_ROW, row)
		 for i = 1:length(mysqlfield_types)
				strval = mysql_load_string_from_resultptr(result, i)
				if strval == Nothing
					arr[i][row] = NA
				else
					arr[i][row] = mysql_interpret_field(strval,mysql_get_julia_type(mysqlfield_types[i]))
				end
			end
		end

		function mysql_result_to_columns(result::MYSQL_RES)
			nfields = mysql_num_fields(result)
			fields = mysql_fetch_fields(result)
			nrows = mysql_num_rows(result)

			jfield_types = Array{DataType}( nfields)
			mysqlfield_types = Array{Cuint}( nfields)
			field_headers = Array{Symbol}( nfields)
			resultingArray=Array{Any}(0)

			for i = 1:nfields
				mysql_field = unsafe_load(fields, i)
				jfield_types[i] = mysql_get_julia_type(mysql_field.field_type)
				field_headers[i] = Symbol(bytestring(mysql_field.name))
				mysqlfield_types[i] = mysql_field.field_type
				push!(resultingArray,Array{jfield_types[i]}(nrows))
			end

			for row = 1:nrows
				mypopulate_row!(resultingArray, mysqlfield_types, mysql_fetch_row(result), row)
			end

			return resultingArray,field_headers
		end

		function mysql_execute_query_DF_wo_NA(con,command)
			response = mysql_query(con, command)
			mysql_display_error(con, response,"Error occured while executing mysql_query on \"$command\"")

			result = mysql_store_result(con)
			if (result == C_NULL)
				affectedRows = mysql_affected_rows(con)
				return convert(Int, affectedRows)
			end

			#alternative implementation with Arrays instead of DataArrays
			columnarr,field_headers = mysql_result_to_columns(result)
			retval=DataFrame(columnarr,field_headers)

			mysql_free_result(result)
			return retval
		end
=#

function connectsqlite(file)
	local sqldb
	try
		sqldb=SQLiteDB(file)
	catch eri
		@show eri
		error("Unable to connect to SQLiteDB: $(file)")
		return nothing
	end
	return sqldb
end

function querysqlite(sqldb,;tblname="tmp",rows_to_read::Int=-1)
	local qres
	#@show tblname
	#@show sqldb.file
	try
		if rows_to_read>0
			qres=queryWithTypes(sqldb, string("select * from ",tblname," LIMIT ",rows_to_read));
		else
			qres=queryWithTypes(sqldb, string("select * from ",tblname));
		end
	catch eri
		#close db in any case
		close(sqldb)
		@show eri
		error("Unable to perform query on SQLiteDB.")
	end
	#close(sqldb) #close db in any case
	return qres
end

function queryToDF(qres)
		local df
		try
			df=DataFrame(qres.values, Symbol[Symbol(uppercase(x)) for x in qres.colnames]);
		catch eri
			@show eri
			error("Unable to convert SQLQuery to DataFrame")
			return nothing
		end
		return df
end

function readDFfromSQLDB()
	pw="as212.11asdfoope002319..,cvMMndx7"
	local con,df
	try
		con = mysql_connect("localhost", "root", pw, "julia")
		#the following line enforces UTF8 format for this connection
		#there may be nicer solutions for this (e.g. the my.ini file) but for now this is how we solved umlaute issues
		command="""SET NAMES 'utf8';"""
		mysql_execute_query(con,command);
		print("Connected to SQL DB...")
	catch err
		@show err
		error("Unable to connect to SQL DB.")
	end
	#command = """SELECT * from tmp LIMIT 12;"""
	command = """SELECT * from tmp;"""
	try
		@info "TODO: Read performace may be imrpved with: M:\\Non-Life\\Personal\\Bernhard\\Resources\\IT\\Julia\\PackageModifications\\ImproveReadPerformance_MySQL"
		df=mysql_execute_query_DF_wo_NA(con,command);
	catch err
		mysql_disconnect(con)
		@show err
		error("Unable to read from SQL DB.")
	end
	mysql_disconnect(con)

	#if the first column (irkey) is a numeric variable in SAS, then it will be read as Float64 here.
	#further it will have digits which will yield an output (in the end) of the form [1.0, 2.0, 3.0]
	#to suppress this we try to convert the first column to an Int column
	#if it fails there are either digits (which are nonzero) or it is indeed a string
		try
			if eltype(df[1])<:Number
				df[1]=convert(Array{Int,1},df[1])
			end
		catch
			#we don't care
		end
return df
end

function readDFfromSQLiteDB(file,const_shift_cols::Int,number_of_num_features::Int)
	#determine eltype by reading only a few rows
		rows_to_read=100
		sqldb=connectsqlite(file)
		print("Connected to SqliteDB at $(file)")
		local qres,tblname,df
		qres=querysqlite(sqldb,rows_to_read=rows_to_read)
		#close db in any case
		#close(sqldb)
		#create df
		df=queryToDF(qres)
		#modify_eltypes_by_guessing!(df)
		eltypev=define_eltypevector(df,const_shift_cols,number_of_num_features)
	#read full file
		#sqldb=connectsqlite(file)
		#print("Connected to SqliteDB at $(file)")
		qres=querysqlite(sqldb)
		print("Closing SqliteDB")
		close(sqldb) #close db in any case
		#create df
		print("Converting query result to  DataFrame")
		df=queryToDF(qres)
		#warn("this is obsolete now!	 change it")
		try
			modify_eltypes!(df,eltypev)
		catch eri
			@show eri
			error("Unable to convert datatypes of DataFrame after reading from SQLiteDB.")
		end
	return df
end

function variablesUsed(bt)
	utf8ListVarsUsed=Array{String}(0)
	niter=bt.settings.niter
	nmat=zeros(Int,bt.settings.number_of_num_features+bt.settings.number_of_char_features,niter)
	for i=1:niter
		#note: bt.intNumVarsUsed[i] is an array denoting which split threshold or id was used
		for j=1:length(bt.intVarsUsed[i])
			nmat[j,i]+=sum(bt.intVarsUsed[i][j])
		end
	end
	vars_used=sum(nmat,2).>0
	return deepcopy(bt.settings.df_name_vector[vars_used[:]]),deepcopy(vars_used)
end

function update_matched_vectors!(trnidx,validx,trnMatched,valMatched,feat,thresh)
	fill!(trnMatched,false)
	fill!(valMatched,false)
	for x in validx
		#@inbounds valMatched[x] = ifelse(feat[x].==thresh,true,false)
		@inbounds valMatched[x] = ifelse(feat[x].==thresh,true,false)
	end
	for x in trnidx
		@inbounds trnMatched[x] = ifelse(feat[x].==thresh,true,false)
	end
	return nothing
end

function createTwoWayValidationCharts(trnidx,validx,nameOfValidationSheet,scoreBandLabels::Array{String,1},mappings,candMatWOMaxValues,sett,scores,actualNumerator,denominator,weight,finalEstimateForCharts,features)
denominatortrn=view(denominator,trnidx)

trnMatched=BitVector(length(actualNumerator))
valMatched=BitVector(length(actualNumerator))

locrow=3
locincrement=22
locrow-=locincrement
#loop through vars
local loc,thismat,emptyline,resultingMatrix,statsThisIteration,singleRowWithKeyMetrics,columnOfRelativityTrn,strVarname,list,intVarnamePOSITIVE,feat,featVAL,outMatrix
emptyOutMatrix=repmat([],0,16) #note the 16 is hardcoded here which is very bad programming, it is the width of the ouput "statsThisIteration" of the fn createTrnValStatsForThisIteration
resultingMatrix=deepcopy(emptyOutMatrix)
integerVarlist=sett.statsByVariables
validationCharts=Array{Chart}(0)
if sett.statsRandomByVariable>0
	#random groups are requested
	integerVarlist=vcat(integerVarlist,[0])
end
firstdatarow=4
for intVarname in integerVarlist
	outMatrix=deepcopy(emptyOutMatrix)
	if intVarname==0 # 0 refers to the random group
		strVarname=string("Random Group")
		rng=1:sett.statsRandomByVariable
		list=UInt8[rng;]
		feat=convert(Vector{UInt8},rand(rng,length(weight)))
		@assert eltype(feat)==UInt8	
	elseif intVarname>0
		strVarname=sett.df_name_vector[intVarname]
		feat=features[intVarname]
		list=feat.pool #list=features[intVarname+sett.number_of_char_feautres].pool #candMatWOMaxValues[intVarname]
	end
	#create stats
	for thresh in list
		strthresh=string(thresh)
		#trn_matched=(view(feat,trnidx).==thresh)	
		#val_matched=(view(feat,validx).==thresh)
		update_matched_vectors!(trnidx,validx,trnMatched,valMatched,feat,thresh)
				
		#=
		v1=view(view(actualNumerator,trnidx),trn_matched)
		v2=view(view(denominator,trnidx),trn_matched)
		v3=view(view(weight,trnidx),trn_matched)
		v4=view(view(finalEstimateForCharts,trnidx),trn_matched)
		v5=view(view(scores,trnidx),trn_matched)
		v6=view(view(actualNumerator,validx),val_matched)
		v7=view(view(denominator,validx),val_matched)
		v8=view(view(weight,validx),val_matched)
		v9=view(view(finalEstimateForCharts,validx),val_matched)
		v10=view(view(scores,validx),val_matched)
		=#

		v1=view(actualNumerator,trnMatched)
		v2=view(denominator,trnMatched)
		v3=view(weight,trnMatched)
		v4=view(finalEstimateForCharts,trnMatched)
		v5=view(scores,trnMatched)
		v6=view(actualNumerator,valMatched)
		v7=view(denominator,valMatched)
		v8=view(weight,valMatched)
		v9=view(finalEstimateForCharts,valMatched)
		v10=view(scores,valMatched)
		
		#=
			@assert all(v1.==v1b)
			@assert all(v2.==v2b)
			@assert all(v3.==v3b)
			@assert all(v4.==v4b)
			@assert all(v5.==v5b)
			@assert all(v6.==v6b)
			@assert all(v7.==v7b)
			@assert all(v8.==v8b)
			@assert all(v9.==v9b)
			@assert all(v10.==v10b)
		=#
		statsThisIteration,singleRowWithKeyMetrics,columnOfRelativityTrn=createTrnValStatsForThisIteration(scoreBandLabels,-1,sett.scorebandsstartingpoints,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,sett.boolCalculateGini)
		statsThisIteration[1,1]="Scoreband" #this cell has the value "Cumulative Stats n=-1" by default which is not useful here.
		statsThisIteration=hcat(statsThisIteration[:,1],vcat([strVarname],repmat([strthresh],size(statsThisIteration,1)-1,1)),statsThisIteration[:,2:end])
		outMatrix=vcat(deepcopy(outMatrix),deepcopy(statsThisIteration))
	end
	#trn charts
		thischartTRN=defineTwoWayCharts(nameOfValidationSheet,nameOfValidationSheet,string("Q",locrow+=locincrement),"line",firstdatarow,length(scoreBandLabels)-1,9,length(list),2,"Scoreband","Observed Ratio",string("Observed Ratio Trn by ",convert(typeof(nameOfValidationSheet),strVarname)))
		thischartTRNexp=defineTwoWayCharts(nameOfValidationSheet,nameOfValidationSheet,string("Q",locrow+=locincrement),"line",firstdatarow,length(scoreBandLabels)-1,3,length(list),2,"Scoreband","Weight",string("Weight Trn by ",convert(typeof(nameOfValidationSheet),strVarname)))
	#val charts
		thischartVAL=defineTwoWayCharts(nameOfValidationSheet,nameOfValidationSheet,string("Q",locrow+=locincrement),"line",firstdatarow,length(scoreBandLabels)-1,10,length(list),2,"Scoreband","Observed Ratio",string("Observed Ratio Val by ",convert(typeof(nameOfValidationSheet),strVarname)))
		thischartVALexp=defineTwoWayCharts(nameOfValidationSheet,nameOfValidationSheet,string("Q",locrow+=locincrement),"line",firstdatarow,length(scoreBandLabels)-1,4,length(list),2,"Scoreband","Weight",string("Weight Val by ",convert(typeof(nameOfValidationSheet),strVarname)))
		locrow+=4 # Add some empty space between charts of different variables
	#add charts to list
		push!(validationCharts,deepcopy([thischartTRN thischartTRNexp thischartVAL thischartVALexp])...)
		firstdatarow+=length(list)*(1+length(scoreBandLabels))+4
	#add Header to data
		thismat=deepcopy(outMatrix)
		outMatrix=deepcopy(emptyOutMatrix) #reset the matrix
		hd=[string("Validation by ",strVarname)]
		header=hcat(hd,repmat([""],1,size(thismat,2)-1))
		isdefined(:emptyline) ? true : emptyline=repmat([""],1,size(thismat,2))
		thismat=vcat(header,emptyline,thismat,emptyline,emptyline)
		resultingMatrix=vcat(resultingMatrix,thismat)
end

return resultingMatrix,validationCharts
end

function set_leaf_numbers!(t::Tree)
	leaves_of_tree=create_leaves_array(t);
	fittedValuesPerLeaf::Vector{Float64}=map(x->x.fitted,leaves_of_tree)
	leafnrlist=sortperm(sortperm(fittedValuesPerLeaf));
	i=0
	for x in leaves_of_tree
		x.id=leafnrlist[i+=1]
	end
return nothing
end

function modify_eltypes!(df::DataFrame,elt)
	@assert size(df,2)==length(elt)
	for i=1:size(df,2)
		if eltype(df[i]) != elt[i]
			df[i]=myconvertDAType(df[i],elt[i])
			# df[i]=convert(DataArray{elt[i],1},df[i])
		end
	end
	return nothing
end

function myconvertDAType(d::D,e::T) where {D <: PooledArray,T <: DataType}
	@assert size(d,2)==1 #this can probably be done in a  nicer way d::DataArry{sometype,1}
	r= e<:AbstractString ? PooledArray(e[convert(e,string(x)) for x in d]) : convert(PooledArray{e,1},d)
	return r
end

function myconvertDAType(d::D,e::T) where {D <: Array,T <: DataType}
	@assert size(d,2)==1 #this can probably be done in a  nicer way d::DataArry{sometype,1}
	r= e<:AbstractString ? Array(e[convert(e,string(x)) for x in d]) : convert(Array{e,1},d)
	return r
end

function modify_eltypes_by_guessing!(df::DataFrame)
	#assumption is that the input df has columns which have eltype Any
	local col,canbefloat
	canbefloat=false
	for i=1:size(df,2)
		col=df[i]
		canbefloat=false
		try
			convert(Array{Float64,1},col)
			canbefloat=true
		catch
			#nothing
		end
		#@show i,canbefloat
		if canbefloat
			df[i]=convert(Array{Float64,1},col)
		else
			df[i]=convert(Array{String,1},col)
			#df[i]=String[string(x) for x in df[i]]
		end
	end
	return nothing
end

function define_eltypevector(df::DataFrame,const_shift_cols::Int64,number_of_num_features::Int64)
#Output is an eltypevector with elements Float64 and String
  #const_shift_cols=6
  cols=size(df,2)
  eltypevector=Array{DataType}(cols)
  @assert cols>6 #we have 6 default columns and need at least 1 independent variable
  #the first five columns are defined by the desing of the algorithm and thus always have the same type
  eltypevector[1]=String
  eltypevector[2]=Float64
  eltypevector[3]=Float64
  eltypevector[4]=Float64
  eltypevector[5]=Int64
  eltypevector[6]=String
  #tbd/todo if numvars need to come first: add a check that checks if only charvars appear towards the end
  i2=0
  for i=const_shift_cols+1:cols
	i2+=1
	if i2<=number_of_num_features
		eltypevector[i]=Float64
	else
		eltypevector[i]=String
	end
    #=
		if eltype(df[i])==Bool #todo/tbd improve how booleans are handled
			eltypevector[i]=String
		else
			if typeof(df[1,i]) <: Number
				eltypevector[i]=Float64
			else
				eltypevector[i]=String
			end
		end
	=#
  end
  return eltypevector
end

function convertDFbooleanColsToUTF8String!(d::DataFrame)
cols=size(d,2)
for i=1:cols
	@show i
	#@show eltype(d[i])
	#todo/tbd improve this
	#this is nasty and probably inefficient hack
	@show typeof(d[1,i])
	@show eltype(d[i])
	if typeof(d[1,i])==Bool
		@show 1
		dtmp=convert(Array,d[i])
		#@show typeof(dtmp)
		#@show eltype(dtmp)
		@show 2
		dtmp2=String[boolToUTF8String(x)  for x in dtmp]
		@show 3
		da=PooledArray(dtmp2)
		@show 4
		d[i]=deepcopy(da)
		@show 5
	end
end
return nothing
end

function boolToUTF8String(x::Bool)
	t=String(string("true"))
	f=String(string("false"))
	if x==true
		return t
	else
		return f
	end
end

function writeStatisticsFile!(statsfileExcel,xlData,filelistWithFilesToBeZipped)
	println("Exporting Stats to Excel file: \r\n $(statsfileExcel)")
	try 
		isfile(statsfileExcel)&&rm(statsfileExcel)
		@time write_statistics(xlData,statsfileExcel,false,false)
		push!(filelistWithFilesToBeZipped,statsfileExcel)
	catch e
        @show e
        dump(e)
		warn("DTM: Failed to create Excel Statistics file. \r\n $(statsfileExcel)")
	end
	return nothing
end

function removeBOM(s::Array{T,2}) where {T}
	#remove Byte Order Mark https://en.wikipedia.org/wiki/Byte_order_mark
	#if first(s[1])=='\ufeff' #parse(Int,first(s[1]))==65279
	if first(s[1])==global_byte_order_mark #'\ufeff' #parse(Int,first(s[1]))==65279
		tmpchart,tmpidx=next(s[1],1)
		copyOfS=deepcopy(s) #copying is likely not the most efficient way here...
		copyOfS[1]=s[1][tmpidx:end]
		return copyOfS
	end
	return s
end


function writeAllFields(f::IOStream,s)
	for x in fieldnames(s)
		write(f,x,"\r\n  ",string(getfield(s,x)),"\r\n")
	end
	return nothing
end

function writeAllFieldsToString(s)
	res=string()
	for x in fieldnames(s)
		res=string(res,x,"\r\n  ",string(getfield(s,x)),"\r\n")
	end
	return res
end

function writeAllFieldsToArray(s)
	nn=fieldnames(s)
	res=Array{AbstractString}(length(nn),2)
	i=0
	for x in nn
		i+=1
		res[i,:]=AbstractString[string(x) string(getfield(s,x))][:]
	end
	return res
end

function absrel(x,y)
	@assert length(x)==length(y)
	res=zeros(x)
	for i=1:length(x)
		@inbounds xi=x[i]
		@inbounds yi=y[i]
		r=abs(xi-yi)/xi
		@inbounds res[i]=r
	end
	return res
end

function histmod(v::AbstractVector,edg::AbstractVector)
	#last bucket (i.e. elements exceeding edge[end]) is not counted
	#a,b=hist(v,edg) #deprecated since julia 0.5
	hfit=StatsBase.fit(Histogram,v,edg,closed=:right)
	a=hfit.edges[1]
	b=hfit.weights
	push!(b,length(v)-sum(b))
	incrementalToCumulative!(b)
	return a,b
end

function errorhistogram(histbins,actualNumerator,estimateUnsmoothedTRN,estimateSmoothedTRN,estimateFromRelativitiesTRN,actualNumeratorVAL,estimateUnsmoothedVAL,estimateSmoothedVAL,estimateFromRelativitiesVAL)
	nTRN=length(actualNumerator)
	nVAL=length(actualNumeratorVAL)
	ed,count1TRN=histmod(absrel(actualNumerator,estimateUnsmoothedTRN),histbins)
	ed,count2TRN=histmod(absrel(actualNumerator,estimateSmoothedTRN),histbins)
	ed,count3TRN=histmod(absrel(actualNumerator,estimateFromRelativitiesTRN),histbins)
	ed,count1VAL=histmod(absrel(actualNumeratorVAL,estimateUnsmoothedVAL),histbins)
	ed,count2VAL=histmod(absrel(actualNumeratorVAL,estimateSmoothedVAL),histbins)
	ed,count3VAL=histmod(absrel(actualNumeratorVAL,estimateFromRelativitiesVAL),histbins)
	errorhist=hcat(histbins,count1TRN,count2TRN,count3TRN,count1VAL,count2VAL,count3VAL)

	#relative count
		#note: the explicit call of the convert funciton is likely obsolete, but it should not hurt either
	count1TRNf=convert(Array{Float64,1},count1TRN./nTRN)
	count2TRNf=convert(Array{Float64,1},count2TRN./nTRN)
	count3TRNf=convert(Array{Float64,1},count3TRN./nTRN)
	count1VALf=convert(Array{Float64,1},count1VAL./nVAL)
	count2VALf=convert(Array{Float64,1},count2VAL./nVAL)
	count3VALf=convert(Array{Float64,1},count3VAL./nVAL)

	errorhist=hcat(errorhist,count1TRNf,count2TRNf,count3TRNf,count1VALf,count2VALf,count3VALf)
return errorhist
end

#errors_num_estimates=addTariffEstimationStatsAndGraphs!(xlData,trnidx,validx,actualNumerator,estimateUnsmoothed,estimateSmoothed,estimateFromRelativities,est_matrix,est_matrixFromScores)
function addTariffEstimationStatsAndGraphs!(xlData,trnidx::Vector{Int},validx::Vector{Int},
	actualNumerator::Array{Float64,1},estimateUnsmoothed::Array{Float64,1},estimateSmoothed::Array{Float64,1},estimateFromRelativities::Array{Float64,1},est_matrix::Array{Float64,2},est_matrixFromScores::Array{Float64,2})
	niter=size(est_matrixFromScores,2)-1
	sumactual=sum(view(actualNumerator,trnidx))
    sumactualVAL=sum(view(actualNumerator,validx))	
	#[name i mse mrse mae mrae sumrae mseVAL mrseVAL maeVAL mraeVAL sumraeVAL][:]
	header=["Summation Type" "Iteration" "MSE TRN" "MRSE TRN" "MAE TRN" "MRAE TRN" "Sum MRAE TRN" "MSE VAL" "MRSE VAL" "MAE VAL" "MRAE VAL" "Sum MRAE VAL"]
	
    m1=tariffEstStats(sumactual,actualNumerator,est_matrix,sumactualVAL,trnidx,validx,"Raw Estimates")    
	m2=tariffEstStats(sumactual,actualNumerator,est_matrixFromScores,sumactualVAL,trnidx,validx,"Estimates from Scores")
	z1=tariffEstStats(sumactual,actualNumerator,estimateUnsmoothed,sumactualVAL,trnidx,validx,"Final Estimate Unsmoothed",0)
	z2=tariffEstStats(sumactual,actualNumerator,estimateSmoothed,sumactualVAL,trnidx,validx,"Final Estimate Smoothed",0)
	z3=tariffEstStats(sumactual,actualNumerator,estimateFromRelativities,sumactualVAL,trnidx,validx,"Final Estimate from Relativities",0)

	result=vcat(header,z1,z2,z3,m1,m2)
	#define errors_num_estimates (which will be attached to modelstatistics)
		#the first two columns of m1 and m2 are obsolete now
		header=header[:,3:end];m1=m1[:,3:end];;m2=m2[:,3:end];
		tmpvector_2_t=[string(x," (Raw Estimates)") for x in header]
		tmpvector_2=reshape(tmpvector_2_t,1,length(tmpvector_2_t))
		tmpvector_2_t=[string(x," (Estimates from Scores)") for x in header]
		tmpvector_3=reshape(tmpvector_2_t,1,length(tmpvector_2_t))
		errors_num_estimates=hcat(vcat(tmpvector_2,m1),vcat(tmpvector_3,m2)) # [:,3:end]
	#add two emtpy rows
		result=vcat(result,repmat([""],2,size(result,2)))
	#todo tbd potentially add error statistics by score band (maybe consider the average estimated premium per score band versus the actual)
	#Error Histogram
	errstartrow=1+size(result,1)
	errorhistheader=["Absolute Relative Error (Bins)" "Estimate Unsmoothed TRN" "Estimate Smoothed TRN" "Estimate from Relativities TRN" "Estimate Unsmoothed VAL" "Estimate Smoothed VAL" "Estimate from Relativities VAL" "Estimate Unsmoothed TRN Pct" "Estimate Smoothed TRN Pct" "Estimate from Relativities TRN Pct" "Estimate Unsmoothed VAL Pct" "Estimate Smoothed VAL Pct" "Estimate from Relativities VAL Pct" ]
	#granular up to .2
		histbins=collect(0.0:0.001:0.2) #steps of 1 per mille
		histbinsChartRows=length(histbins)
		errorhist=errorhistogram(histbins,view(actualNumerator,trnidx),view(estimateUnsmoothed,trnidx),view(estimateSmoothed,trnidx),view(estimateFromRelativities,trnidx),view(actualNumerator,validx),view(estimateUnsmoothed,validx),view(estimateSmoothed,validx),view(estimateFromRelativities,validx))
		errorhist=vcat(errorhistheader,errorhist)
	#granular up to 1
		histbins=collect(0.0:0.001:1.0) #steps of 1 per mille
		errorhist2=errorhistogram(histbins,view(actualNumerator,trnidx),view(estimateUnsmoothed,trnidx),view(estimateSmoothed,trnidx),view(estimateFromRelativities,trnidx),view(actualNumerator,validx),view(estimateUnsmoothed,validx),view(estimateSmoothed,validx),view(estimateFromRelativities,validx))
		errorhist2=vcat(errorhistheader,errorhist2)
	#add two emtpy rows
		errorhist=vcat(errorhist,repmat([""],2,size(errorhist,2)))
		errorhist=vcat(errorhist,errorhist2)
	#attach histogram below errors per iteration
	#increase breadth
		mx=max(size(result,2),size(errorhist,2))
		if mx>size(result,2)
			result=hcat(result,repmat([""],size(result,1),mx-size(result,2)))
		else
			errorhist=hcat(errorhist,repmat([""],size(errorhist,1),mx-size(errorhist,2)))
		end
	#attach histogram
		result=vcat(result,errorhist)
	#add sheet to Excel
		resultDF=convert(DataFrame,result)
		tariffestsheetname="NumeratorEstimationError" #NOTE: if we have a space in the sheetname then the functions which produce the graphs need to be adjusted slighty: the formula must have single quotes around the sheetname
		tariffEstErrorSheet=ExcelSheet(tariffestsheetname,resultDF)
		push!(xlData.sheets,tariffEstErrorSheet)
	#Add Charts
		c1=defineChartWithNSeries(tariffestsheetname,tariffestsheetname,"C5","line",5,niter,[4,4+5],1,"Iteration","Error","Based on Raw Estimates",categoriescol=2,yaxisformat="0.00%",xscale=3.5)
		push!(xlData.charts,deepcopy(c1))
		c1=defineChartWithNSeries(tariffestsheetname,tariffestsheetname,"C27","line",5,niter,[6,6+5],1,"Iteration","Error","Based on Raw Estimates",categoriescol=2,yaxisformat="0.00%",xscale=3.5)
		push!(xlData.charts,deepcopy(c1))
		c1=defineChartWithNSeries(tariffestsheetname,tariffestsheetname,"C49","line",5,niter,[7,7+5],1,"Iteration","Error","Based on Raw Estimates",categoriescol=2,yaxisformat="0.00%",xscale=3.5)
		push!(xlData.charts,deepcopy(c1))

		c1=defineChartWithNSeries(tariffestsheetname,tariffestsheetname,"C71","line",niter+5,niter,[4,4+5],1,"Iteration","Error","Based on Estimates by Score",categoriescol=2,yaxisformat="0.00%",xscale=3.5)
		push!(xlData.charts,deepcopy(c1))
		c1=defineChartWithNSeries(tariffestsheetname,tariffestsheetname,"C93","line",niter+5,niter,[6,6+5],1,"Iteration","Error","Based on Estimates by Score",categoriescol=2,yaxisformat="0.00%",xscale=3.5)
		push!(xlData.charts,deepcopy(c1))
		c1=defineChartWithNSeries(tariffestsheetname,tariffestsheetname,"C113","line",niter+5,niter,[7,7+5],1,"Iteration","Error","Based on Estimates by Score",categoriescol=2,yaxisformat="0.00%",xscale=3.5)
		push!(xlData.charts,deepcopy(c1))
	#error dist charts
		c1=defineChartWithNSeries(tariffestsheetname,tariffestsheetname,"C135","line",errstartrow+1,histbinsChartRows,[8,8+3],errstartrow,"Error","Percentage of Data","Cumulative Distribution Absolute Relative Errors",categoriescol=1,xaxisformat="0.00%",yaxisformat="0.00%",yscale=3.1,xscale=3.5)
		push!(xlData.charts,deepcopy(c1))
		c1=defineChartWithNSeries(tariffestsheetname,tariffestsheetname,"C181","line",errstartrow+1,histbinsChartRows,[9,9+3],errstartrow,"Error","Percentage of Data","Cumulative Distribution Absolute Relative Errors",categoriescol=1,xaxisformat="0.00%",yaxisformat="0.00%",yscale=3.1,xscale=3.5)
		push!(xlData.charts,deepcopy(c1))
		c1=defineChartWithNSeries(tariffestsheetname,tariffestsheetname,"C227","line",errstartrow+1,histbinsChartRows,[10,10+3],errstartrow,"Error","Percentage of Data","Cumulative Distribution Absolute Relative Errors",categoriescol=1,xaxisformat="0.00%",yaxisformat="0.00%",yscale=3.1,xscale=3.5)
		push!(xlData.charts,deepcopy(c1))
	#all in one
		c1=defineChartWithNSeries(tariffestsheetname,tariffestsheetname,"C273","line",errstartrow+1,histbinsChartRows,[8,8+3,9,9+3,10,10+3],errstartrow,"Error","Percentage of Data","Cumulative Distribution Absolute Relative Errors",categoriescol=1,xaxisformat="0.00%",yaxisformat="0.00%",yscale=3.1,xscale=3.5)
		push!(xlData.charts,deepcopy(c1))
return errors_num_estimates
end

function tariffEstStats(sumactual::Float64,actual::Array{Float64,1},estimateMat::Array{Float64,2},sumactualVAL::Float64,trnidx::Vector{Int},validx::Vector{Int},name::T) where {T <: AbstractString}
	result=Array{Any,2}(size(estimateMat,2)-1,12)
	for i=1:size(estimateMat,2)-1 #NOTE the first column of estimateMat is iteration 0 which is the mean fitted value for all observations (which is not interesting at all!)
		mse,mrse,mae,mrae,sumrae=calcErrors(sumactual,view(actual,trnidx),view(estimateMat,trnidx,i+1))        
		mseVAL,mrseVAL,maeVAL,mraeVAL,sumraeVAL=calcErrors(sumactualVAL,view(actual,validx),view(estimateMat,validx,i+1))
		result[i,:]=[name i mse mrse mae mrae sumrae mseVAL mrseVAL maeVAL mraeVAL sumraeVAL][:]
	end
	return result
end

function tariffEstStats(sumactual::Float64,actual::Array{Float64,1},estimate::Array{Float64,1},sumactualVAL::Float64,trnidx::Vector{Int},validx::Vector{Int},name::T,count::Int) where {T <: AbstractString}
	mse,mrse,mae,mrae,sumrae=calcErrors(sumactual,view(actual,trnidx),view(estimate,trnidx))
	mseVAL,mrseVAL,maeVAL,mraeVAL,sumraeVAL=calcErrors(sumactualVAL,view(actual,validx),view(estimate,validx))
	res=[name count mse mrse mae mrae sumrae mseVAL mrseVAL maeVAL mraeVAL sumraeVAL]	
end

function calcErrors(sumactual::Float64,actual,estimate) #::Array{Float64,1},estimate::Array{Float64,1})
	@assert length(estimate)==length(actual)
	mae=mrse=mrae=mse=0.0
	sz=length(actual)
	for i=1:sz
			@inbounds act=actual[i]
			@inbounds est=estimate[i]
			diff=est-act
			diffRel=diff/est
			mse+=abs2(diff)
			mae+=abs(diff)
			mrse+=abs2(diffRel)
			mrae+=abs(diffRel)
		end
		sumrae=mae/sumactual
		mse/=sz
		mae/=sz
		mrse/=sz
		mrae/=sz
	return mse,mrse,mae,mrae,sumrae
end

function initBootstrapSample(sampleSizeCanBeNEGATIVE::Int) #,trnidx::Vector{Int},validx::Vector{Int})
		abssampleSize=abs(sampleSizeCanBeNEGATIVE)
		@assert abssampleSize>0
		sampleVector=Vector{Int}(abssampleSize)
	return abssampleSize,sampleVector
end
    
function defineWeightPerTree(vectorOfAllTreesWithErrorStats::Array{TreeWithErrorStats,1},baggingWeightTreesError::String)
	local res
	if baggingWeightTreesError=="uniform"
		res=ones(Float64,length(vectorOfAllTreesWithErrorStats))
	else
		res=[getfield(x.err,Symbol(baggingWeightTreesError)) for x in vectorOfAllTreesWithErrorStats]
	end
	#["mse" "pearsoncorrelation" "mrae" "mrse" "mae" "uniform"]
	#weightPerTree=[x.err.mae for x in vectorOfAllTreesWithErrorStats]
	customInverseNormalized!(res)
	if sum(isfinite.(res))==0
		#if all errors are set to 0 (or NaN, Inf)
		return ones(Float64,length(vectorOfAllTreesWithErrorStats))*1/length(vectorOfAllTreesWithErrorStats)
	else
		return res
	end
end

function customInverseNormalized!(weightPerTree::Array{T,1}) where {T <: Number}
	#todo, this can be done more efficiently
	for i=1:length(weightPerTree)
		@inbounds weightPerTree[i]=1/weightPerTree[i]
	end
	mx=-Inf
	for i=1:length(weightPerTree)
		@inbounds w=weightPerTree[i]
		if isfinite(w)
			if w>mx
				mx=w
			end
		end
	end
	for i=1:length(weightPerTree)
		@inbounds w=weightPerTree[i]
		if !isfinite(w)
			@inbounds weightPerTree[i]=mx
		end
	end
	sm=sum(weightPerTree)
	for i=1:length(weightPerTree)
		@inbounds weightPerTree[i]/=sm
	end
	return nothing
end

function aggregateByLeafNr(leafNrVector::Array{Int,1},numeratorEST::Array{Float64,1},numerator::Array{Float64,1},denominator::Array{Float64,1},weight::Array{Float64,1})
	(lo, hi) = extrema(leafNrVector) #this should by definition always be 1,length(leafNrVector), but the evaluation of this function should be nearly instant as normally leafNrVector has at most a few hundered entries
	ooo=one(lo)-lo
	vecsize=hi+ooo
    cnt = zeros(Int, vecsize)
    sumnumeratorEST = zeros(Float64, vecsize)
	sumnumerator = zeros(Float64, vecsize)
	sumdenominator = zeros(Float64, vecsize)
	sumweight= zeros(Float64, vecsize)
    #note: inbounds increases efficiency here (about a factor of 2), however if the bounds are violated something nasty might happen (quote: If the subscripts are ever out of bounds, you may suffer crashes or silent corruption.)
	for count=1:length(leafNrVector)
        #i=a[count] #i in a
		#idx=i - lo + 1
        @inbounds idx=leafNrVector[count] + ooo
		@inbounds cnt[idx] += 1
        @inbounds sumnumerator[idx] += numerator[count]
		@inbounds sumnumeratorEST[idx] += numeratorEST[count]
		@inbounds sumdenominator[idx] += denominator[count]
		@inbounds sumweight[idx] += weight[count]
    end
    return cnt,sumnumeratorEST,sumnumerator,sumdenominator,sumweight
end


function corw(x::Array{T,1},y::Array{U,1},wvec::W) where {T <: Number,U <: Number,W <: AbstractWeights}
  @assert length(x)==length(y)==length(wvec)
  #todo this can be done more efficiently
  #soon it will probably be available in StatsBase or in some other package
  w=wvec.values
  meanwx=mean(x,wvec)
  meanwy=mean(y,wvec)
  #xbar=x.-meanwx
  #ybar=y.-meanwy
  res=corwInternal(w,x,y,meanwx,meanwy)
  return res
end

function corwInternal(w,x,y,meanwx,meanwy)
	n = length(w)
    xi=x[1]-meanwx
	yi=y[1]-meanwy
	xx = w[1]*abs2(xi)
	yy = w[1]*abs2(yi)
    xy = w[1]*xi*yi
    i = 1
    while i < n
        i += 1
        @inbounds xi = x[i]-meanwx
        @inbounds yi = y[i]-meanwy
		@inbounds wi = w[i]
        xx += abs2(xi)*wi
        yy += abs2(yi)*wi
        xy += wi*xi*yi
    end
    return xy / (sqrt(xx*yy))
end

function calcErrorStats(fitted::Array{Float64,1},actual::Array{Float64,1},weight::Array{Float64,1})
#need to change/improve this, the measures are probably not meaningful for ratios IF the data is not aggreated in some manner
	@assert length(weight)==length(actual)==length(fitted)
	#this function can probably be optimized quited a bit as many of the quantities could be derived in one go (see also calcErrors where we calculate UNWEIGHTED errors for TariffEstimation Statistics)
	#note stats below could be calculated more efficiently, however if we only have 20 groups runtime should not be any issue
		weightVec=FrequencyWeights(weight)
		diff=fitted.-actual
	#least squares devaiation, mean square error
		msevec=abs2(diff)
		mse=mean(msevec,weightVec)
	#average deviation, mean absolute error
		maevec=abs(diff)
		mae=mean(maevec,weightVec)
	#relative squared error,mean realtive squared error
		diffRel=diff./fitted
		mrsevec=abs2(diffRel)
		mrse=mean(mrsevec,weightVec)
	#relative absolute deviation, mean relative absolute error
		mraevec=abs(diffRel)
		mrae=mean(mraevec,weightVec)
	#Pearson product moment correlation / correlation coefficient
		pearsoncorrelation=corw(fitted,actual,weightVec)
		res=ErrorStats(mse,mae,mrse,mrae,pearsoncorrelation)
return res
end

function goodnessOfFit(idx::Vector{Int},numESTIMATE,numACTUAL,w,d)
	#todo/tbd	note: this needs additional work
	#currently the data is aggregated into (up to) 20 random groups here
	obs=length(idx)
	ngroups_max=20
	f_randgroups=PooledArray(rand(1.0:float(ngroups_max),obs))
	not_meaningful_scores=zeros(Int,obs)
	
	numeratorEst=view(numESTIMATE,idx)
	numeratorActual=view(numACTUAL,idx)
	weight=view(w,idx)
	denominator=view(d,idx)
	error("this has not yet been updated")
	valuelist,cnt,NOT_MEANINGFUL_UNUSEDsumscores,sumnumeratorEst,sumnumerator,sumdenominator,sumweight=aggregate_data(f_randgroups,not_meaningful_scores,numeratorEst,numeratorActual,denominator,weight)
	#valuelist,cnt,NOT_MEANINGFUL_UNUSEDsumscores,sumnumeratorEst,sumnumerator,sumdenominator,sumweight=aggregate_data(randgroups_pda,not_meaningful_scores,numESTIMATE,numACTUAL,denominator,weight)
	
	return calcErrorStats(sumnumeratorEst./sumdenominator,sumnumerator./sumdenominator,sumweight)
end

function calcIterationsPerCore(ntot::Int,ncores::Int)
	@assert ncores>0
	@assert ntot>0
	if ncores==1
		return [ntot]
	else
		nPerCore=fld(ntot,ncores)
		res=ones(Int,ncores).*nPerCore
		res[end]=ntot-sum(res[1:end-1])
		excess=res[end,1]-res[1] #by definition this must be smaller than ncores
		add_this=vcat(ones(Int,excess),zeros(Int,ncores-excess))
		res=res.+add_this
		res[end]=ntot-sum(res[1:end-1])
		@assert sum(res)==ntot
		return res
	end
end

function calcIterationsPerCore2(ntot::Int,ncores::Int)
#this function returns a cumulative counter in the second column of res
#the counter states the starting point of the iteration for the respective core
	@assert ncores>0
	@assert ntot>0
	if ncores==1
		return Int[ntot 1]
	else
		nPerCore=fld(ntot,ncores)
		res=ones(Int,ncores).*nPerCore
		res[end]=ntot-sum(res[1:end-1])
		excess=res[end,1]-res[1] #by definition this must be smaller than ncores
		add_this=vcat(ones(Int,excess),zeros(Int,ncores-excess))
		res=res.+add_this
		res[end]=ntot-sum(res[1:end-1])

		col2=deepcopy(res)
		incrementalToCumulative!(col2)
		unshift!(col2,0)
		pop!(col2)
		@assert sum(res)==ntot
		res2=hcat(res,col2.+1)
		return res2
	end
end

function sampleData!(trnidx::Vector{Int},samplesize::Int,smp::Vector{Int}) #,numerator::Array{Float64,1},denominator::Array{Float64,1},weight::Array{Float64,1}, numfeatures::Array{pdaMod,1},charfeatures::Array{pdaMod,1},smp::Array{Int,1},n::Array{Float64,1},d::Array{Float64,1},w::Array{Float64,1},numf::Array{pdaMod,1},charf::Array{pdaMod,1},un::Array{Float64,1},ud::Array{Float64,1},uw::Array{Float64,1},unumf::Array{pdaMod,1},ucharf::Array{pdaMod,1})
#this version should (does it?) allocate less memory
		repl=(samplesize<0)
		@assert repl==false "DTM: Sampling without replacement is currently not supported." #bk -> tbd/todo check how we can make this work
		#e.g. fitted_values_all_data_this_vector_is_modified_by_build_tree must possibly look different if replace=true 
		#currently a CRITICAL assumption is that validx and trnidx are disjoint and unique(validx)==validx . The same must holde for the sample and oobag sample
		@assert length(smp)==abs(samplesize)
		sample!(trnidx,smp,replace=repl)
	#define unused part of data
		repl==false ? smpUniqueSorted=smp : smpUniqueSorted=unique(smp)
		sort!(smpUniqueSorted) 

		smpUnusedPart=Vector{eltype(smp)}(0)
		sizehint!(smpUnusedPart,length(trnidx)-length(smpUniqueSorted))
		for x in trnidx
			tmp=searchsorted(smpUniqueSorted,x)
			if length(tmp)<=0
				push!(smpUnusedPart,x)
			end
		end
		@assert size(smpUnusedPart,1)+size(smp,1)==size(trnidx,1)
return smpUnusedPart
end

function sampleData(samplesize::Int,numerator::Array{Float64,1},denominator::Array{Float64,1},weight::Array{Float64,1}, numfeatures::Array{pdaMod,1},charfeatures::Array{pdaMod,1})
		repl=!(samplesize<0)
		smp=sample(1:size(numerator,1),abs(samplesize),replace=repl)
	#draw sample from data, this currently creates a COPY of the data
		n=numerator[smp]
		d=denominator[smp]
		w=weight[smp]
		numf=subset_pda_mod(numfeatures,smp)
		charf=subset_pda_mod(charfeatures,smp)
return smp,n,d,w,numf,charf
end

function sampleData(prop::Float64,numerator::Array{Float64,1},denominator::Array{Float64,1},weight::Array{Float64,1}, numfeatures::Array{pdaMod,1},charfeatures::Array{pdaMod,1})
	samplesize=convert(Int,round(prop*size(numerator,1)))
	return sampleData(samplesize,numerator,denominator,weight,numfeatures,charfeatures)
end

function initSettingsWhichAreTheSameForBoostingAndBagging(trnidx::Vector{Int},validx::Vector{Int},actualNumerator::Array{Float64,1},denominator::Array{Float64,1},sett::ModelSettings)
	empty_rows_after_iteration_stats=2
	showTimeUsedByEachIteration=sett.showTimeUsedByEachIteration
	chosen_apply_tree_fn=sett.chosen_apply_tree_fn
	BoolStartAtMean=sett.BoolStartAtMean
	adaptiveLearningRate=sett.adaptiveLearningRate
	moderationvector=sett.moderationvector
	iterations=sett.niter

	obstrn = length(trnidx) #originally obs
	obsval = length(validx) 
	trn_numtot=sum(view(actualNumerator,trnidx))
	trn_denomtot=sum(view(denominator,trnidx))
	trn_meanobservedvalue=trn_numtot/trn_denomtot
	val_numtot=sum(view(actualNumerator,validx))
	
	val_denomtot=sum(view(denominator,validx))
	val_meanobservedvalue=val_numtot/val_denomtot
	nameOflistofpredictorsSheet,nameOfpredictorsSheet,nameOfModelStatisticsSheet,nameOfScoresSheet,nameOfOverviewSheet,nameOfSettingsSheet,nameOfValidationSheet,xlData=initExcelData()
	showProgressBar_time=sett.showProgressBar_time
	return obstrn,obsval,trn_meanobservedvalue,val_meanobservedvalue,trn_numtot,val_numtot,trn_denomtot,val_denomtot,empty_rows_after_iteration_stats,showTimeUsedByEachIteration,chosen_apply_tree_fn,BoolStartAtMean,adaptiveLearningRate,moderationvector,iterations,nameOflistofpredictorsSheet,nameOfpredictorsSheet,nameOfModelStatisticsSheet,nameOfScoresSheet,nameOfOverviewSheet,nameOfSettingsSheet,nameOfValidationSheet,xlData,showProgressBar_time,sett.boolProduceEstAndLeafMatrices
end

function initExcelData()
#Define ExcelSheet
	xlData=ExcelData(Array{ExcelSheet}(0),Array{Chart}(0))
	#nameOfModelStatisticsSheet=global_nameOfModelStatisticsSheet
	#nameOfSettingsSheet=global_nameOfSettingsSheet
	nameOfScoresSheet="Scores"
	nameOfOverviewSheet="Overview"
	nameOfpredictorsSheet="Predictors"
	nameOflistofpredictorsSheet="PredictorList"
return nameOflistofpredictorsSheet,nameOfpredictorsSheet,global_nameOfModelStatisticsSheet,nameOfScoresSheet,nameOfOverviewSheet,global_nameOfSettingsSheet,global_nameOfTWOWayValidationSheet,xlData
end

macro nogc(ex)
  quote
    try
      gc_disable()
      local val = $(esc(ex))
    finally
      gc_enable()
    end
    val
  end
end

macro timeConditional(bool,ex)
	return :(if $(esc(bool));@time ($(esc(ex)));else;$(esc(ex));end)
end

#=
    function leaf_numbers(leaves::Array{Leaf,1},trnsize::Int)
        leafnr=zeros(Int,trnsize)
        for i=1:length(leaves)
            thisid=leaves[i].id
            for j in leaves[i].idx
                leafnr[j]=thisid
            end
        end
        return leafnr
    end
=#

"""
returns l and r indices (l corresponds to all elements of trnidx that 'match' subset)
This verison of the function is for String / Categorical variables
"""
function lrIndices(trnidx::Vector{Int},f::PooledArray{String,T,1},subset::Array) where T<:Unsigned
	l=Vector{Int}(0)
	r=Vector{Int}(0)
	sizehint!(l,length(trnidx))
	sizehint!(r,length(trnidx))
	for i in trnidx
		@inbounds thisref=f.refs[i]
		if in(thisref,subset)
			push!(l,i)
		else
			push!(r,i)
		end
	end
	return l,r
end

"""
returns l and r indices (l corresponds to all elements of trnidx that 'match' subset)
This verison of the function is for numerical variables
"""
#todo tbd, we could restrict the input to U<:Number and T<:Unsigned for the pda f
function lrIndices(idx::Vector{Int},f,subset::Array) 
	l=Vector{Int}(0)
	r=Vector{Int}(0)
	sizehint!(l,length(idx))
	sizehint!(r,length(idx))
	for i in idx
		@inbounds thisref=f.refs[i]
		#warn("this can be done much more efficiently for  numerical variables! <=")
        #if !(in(thisref,subset)==(thisref<=subset[end]))
        #    error("I did not expect this....")
        #end
		if  thisref<=subset[end] # this should be equivalent to 'in(thisref,subset)' by construction
			push!(l,i)
		else
			push!(r,i)
		end
	end
	return l,r
end

function update_part_of_this_vector!(idx::Vector{Int},read::Vector{T},write::Vector{T}) where T<:Number
	for k in idx
		@inbounds write[k]=read[k]
	end
	return nothing
end

function createIntegerIndices(trn_val_idx)
	l=length(trn_val_idx)
	t=Vector{Int}(0)
	v=Vector{Int}(0)
	sizehint!(v,l)
	sizehint!(t,l)
	for i=1:l 
		if trn_val_idx[i]==1
			push!(t,i)
		else 
			push!(v,i)
		end
	end
	return t,v
end

function randint(n::Int)
	return convert(Int,ceil(rand()*abs(n)))
end

function createZipFile(zipFilename::T,filelist::Array{U,1},silent::Bool=false) where {T <: AbstractString,U <: AbstractString}
#use 7zip to zip the result files
  if isfile(zipFilename)
    try
      rm(zipFilename)
    catch err
      @show err
      warn("Unable to delete file $(zipFilename)!")
    end
  end

    #flist=join(filelist," ")
    #mx5 is compression level 5 http://sevenzip.osdn.jp/chm/cmdline/ http://sevenzip.osdn.jp/chm/cmdline/switches/method.htm#LZMA2
	cmd=`7z a -mx5 -m0=LZMA2:d=25 $(zipFilename) $(filelist)`
	if silent
		tmptxt="tmplog.log"
		tmperr="tmperr.log"
		pip=pipeline(cmd, stdout=tmptxt, stderr=tmperr)
		rc=run(pip)
		isfile(tmperr)&&rm(tmperr)
		isfile(tmptxt)&&rm(tmptxt)
	else
		rc=run(cmd)
	end
  return rc
end
createZipFile(zipFilename::T,filelist::U) where {T <: AbstractString,U <: AbstractString}=createZipFile(zipFilename,[filelist])

function mergeLeftAndRightFittedValues(fitted_valuesl::Array{T,1},fitted_valuesr::Array{T,1},left::BitVector) where {T <: Number}
	  fitted_values=Array{eltype(fitted_valuesl)}(length(left))
	  loopl=loopr=1
	  @inbounds for jj=1:length(left)
		if left[jj]
			fitted_values[jj]=fitted_valuesl[loopl]
			loopl+=1
		else
			fitted_values[jj]=fitted_valuesr[loopr]
			loopr+=1
			end
		end
		return fitted_values
end

	function coded_numdata(numfeatures::Array{Any,2},trn::BitVector,val::BitVector,candidates::Array{Float64,2})
error("The version where a candidate matrix is supplied does not yet work (is not yet implented).")
return "miro"
end


"""adds numerical explanatory vars as PDAs to the DF f"""
function add_coded_numdata!(wholeDF::DataFrame,sett::ModelSettings,trn_val_idx::Vector{UInt8},max_splitting_points_num::Int64,features::DataFrame) #Version where no candidates are supplied -> Julia chooses the candidates
  nobs=length(trn_val_idx)
  cols=sett.number_of_num_features 
  
  local this_column #need to initialize the variable here otherwise it will only exist locally in the try loop
  candMatWOMaxValues=Array{Array{Float64,1}}(0)
  for i=1:cols
	thisname=Symbol(sett.df_name_vector[i])	
	original_data_vector=wholeDF[thisname] 	
	#thiscol_as_utf8=view(dfIndata,:,colnumber)
    header=deepcopy(sett.df_name_vector) #[string(x) for x in names(numfeatures)]
    elt=eltype(original_data_vector)
	@assert (elt<:Number) "Datatype of numeric predictor $(i):$(header[i]) is not numeric in the CSV. Please check the order of the columns in the data, their type and the value of number_of_num_features in the settings!" #the assert is probably not needed here, as we try to convert afterwards anyway
	cond2=(typeof(original_data_vector)<:AbstractVector{Float64})	 
    if !cond2
		try
		  this_column=convert(Array{Float64,1},original_data_vector)
		catch thiserr
		  @show thiserr
		  dump(thiserr)
		  error("Unable to convert column $(i):$(header[i]) to a floating point vector.")
		end
    else
		this_column=original_data_vector
	end
	#Generate candidate list
	candlist=define_candidates(this_column,max_splitting_points_num)
	if max_splitting_points_num>500
		warn("DTM: max_splitting_points_num=$(max_splitting_points_num) is larger than 500. That may not be a good idea")
		if max_splitting_points_num>2000
			error("DTM: max_splitting_points_num too large max_splitting_points_num=$(max_splitting_points_num)")
		end
	end
	(length(candlist)==1)&&(@info "DTM: Numeric column $(i):$(string(thisname)) (numeric) has zero splitting points.") #throw an error when not splitting point exists, was this intended? does this ever happen?
    sort!(candlist,alg=QuickSort) #although it should already be sorted
    mini,maxi=extrema(this_column) #it IS CRUCIAL that the maximum is taken from the whole vector here! (and not only the training part)
	push!(candMatWOMaxValues,deepcopy(candlist))
    candlist[end]=max(candlist[end],maxi) #set the last value of candlist to the largest observed training (and validation) element -> all datapoints which are larger than the former maximum(candlist) will be identified with the maxi
    @assert length(candlist)<255 "Currently at most 254 splitting points are supported. Please choose less splitting points for numeric column $(i):$(header[i])" #see note below
	#prepare trn data	
	compressed_vector=map_numdata_to_candidates(this_column,candlist)
	features[thisname]=PooledArray(compressed_vector)
end

return candMatWOMaxValues
end


function map_numdata_to_candidates(this_column,candlist)
	mapped_vec=zeros(this_column)
	for ii=1:length(this_column)
		#x=this_column[1]
		mapped_idx=searchsortedfirst(candlist,this_column[ii])
		mapped_vec[ii]=candlist[mapped_idx]
	end
	return mapped_vec
end

function define_candidates(feature_column,max_splitting_points_num::Int64)
	#step 1: get quantiles of data
	domain_i = unique(quantile(feature_column, linspace(0.5/Float64(max_splitting_points_num), 1.0-0.5/Float64(max_splitting_points_num),max_splitting_points_num)))
	@assert size(domain_i,1)<=max_splitting_points_num #that should not happen

	#step 2 ensure each element of domain_i corresponds to an observation
	uq_obs=unique(feature_column);
	sort!(uq_obs);
	result=zeros(domain_i)
	for tmploopvar=1:length(domain_i)
		thisval=domain_i[tmploopvar]
		idx2=searchsortedfirst(uq_obs,thisval)
		if idx2>1
			if abs(uq_obs[idx2-1]-thisval)<abs(uq_obs[idx2]-thisval)
				idx2-=1
			end
		end
		result[tmploopvar]=uq_obs[idx2]
	end
	result_unique=unique(result)
	return result_unique #NOTE! the length of domain_i may be smaller than max_splitting_points_num
end

function deleteAllOutputFiles(datafolder::AbstractString,outfilename::AbstractString)
	#deleteIfExists(string(datafolder,outfilename,".log"))
	deleteIfExists(string(datafolder,outfilename,".zip"))
	deleteIfExists(string(datafolder,outfilename,".jld2"))
	deleteIfExists(string(datafolder,outfilename,".7z"))
	deleteIfExists(string(datafolder,outfilename,".txt"))
	deleteIfExists(string(datafolder,outfilename,".sas"))
	deleteIfExists(string(datafolder,outfilename,".cs"))
	deleteIfExists(string(datafolder,outfilename,".csv"))
	deleteIfExists(string(datafolder,outfilename,".xlsx"))
	deleteIfExists(string(datafolder,outfilename,".xlsm"))
	deleteIfExists(string(datafolder,outfilename,".xlsb"))
	deleteIfExists(string(datafolder,outfilename,".stats.csv"))
	deleteIfExists(string(datafolder,outfilename,".stats.xlsx"))
	deleteIfExists(string(datafolder,outfilename,".roptstats.csv"))
	deleteIfExists(string(datafolder,outfilename,".leafnumbers.csv"))
	deleteIfExists(string(datafolder,outfilename,"_iteration_matrix.csv.jld2"))
	deleteIfExists(string(datafolder,outfilename,"_iteration_matrix.csv"))
	deleteIfExists(string(datafolder,outfilename,"_iteration_matrixFromScores.csv"))
return nothing
end

function my_iswriteable(f)
	if isfile(f)
		return '1'==bits(uperm(f))[end-1]
	end
	return false
end

function deleteIfExists(f::AbstractString)
	!isfile(f)&&(return nothing)
	if !my_iswriteable(f)
		warn("No write access to file $(f)")
	else
	try
		rm(f)
	catch err
		warn("Unable to delete file $(f) It might be open by another application.\r\n\r\n\r\n")
		@show err
	end
	end
end

#see also julia base functions basename, dirname, splitext, splitdir, splitdrive http://docs.julialang.org/en/release-0.4/stdlib/file/
function filename(x::AbstractString)
	idx=search(reverse(x),'\\')
	res= idx>0 ? copy(x)[end-idx+2:end] : copy(x)
	return res
end

function pathname(x::AbstractString)
	idx=search(reverse(x),'\\')
	res= idx>0 ? copy(x)[1:end-idx] : copy(x)
	return res
end

function fileroot(f)
	a,fe=splitext(f)
	return a	
end

fileending(f)=splitext(f)[2]	

function printover(s::AbstractString)
	printover(STDOUT,s,:green)
end

function printover(io::IO, s::AbstractString, color::Symbol = :green)
    if isdefined(Main, :IJulia)
        print(io, "\r" * s)
    else
        print(io, "\u1b[1G")   # go to first column
        print_with_color(color, io, s)
        print(io, "\u1b[K")    # clear the rest of the line
    end
end

function add_backslashes!(f::Array)
	f[:]=[x[end]=='\\' ? x : string(x,"\\") for x in f]
	return nothing
end

function add_backslashes(x)
	if x[end]!='\\'
		return string(x,"\\")
	else
		return x
	end
end

function someRankOptStats(mat::Array{Float64,2},numerator::Array{Float64,1},PremiumOfNextRank::Array{Float64,1})
	@assert size(mat,1)==length(numerator)
	nobs,iters=size(mat)
	iters=iters-1 #iteration 0 is not of interest
	premsum=sum(numerator)
	rklostcountlist=zeros(Int64,iters)
	premincrementlist=zeros(Float64,iters)
	for i=1:iters
		rklostcountlist[i]=sum(mat[:,i+1].>PremiumOfNextRank)
		premincrementlist[i]=sum(mat[:,i+1].-numerator)
	end
	rklpct=rklostcountlist./nobs
	prempct=premincrementlist./premsum
	res=[[1:iters] rklostcountlist premincrementlist rklpct prempct rklpct./prempct]
	header=["Iteration" "Ranks Lost" "Additional Premium" "Ranks Lost Pct" "Additional Premium Pct" "Ratio"]
	res=[header,res,repmat([""],2,length(header))]
	return res
end

function firstmatch(x::Array{T,1},y::T) where {T}
  res=0
  for i=1:length(x)
    if y==x[i]
      res=i
      break
    end
  end
  res
end

#function aggregateScores(num::Array{Float64,1},denom::Array{Float64,1},weight::Array{Float64,1},est::Array{Float64,1},scores::Array{Int64,1},scorebandsstartingpoints::Array{Int64,1})
function aggregateScores(num,denom,weight,est,scores,scorebandsstartingpoints)
	nClasses=length(scorebandsstartingpoints)
	sumnumer = zeros(Float64, nClasses)
	sumdenom = zeros(Float64, nClasses)
	sumweight = zeros(Float64, nClasses)
	sumnumeratorEst = zeros(Float64, nClasses)
	@assert length(num)==length(denom)==length(weight)==length(est)==length(scores)	
	@inbounds for i=1:length(weight)
	#todo, check if @inbounds makes this faster! It should be the case
		class=searchsortedlast(scorebandsstartingpoints,scores[i])
		sumnumer[class]+=num[i]
		sumdenom[class]+=denom[i]
		sumweight[class]+=weight[i]
		sumnumeratorEst[class]+=denom[i]*est[i]
	end
	return sumnumer,sumdenom,sumweight,sumnumeratorEst
end

function createScorebandsUTF8List(scorebandsstartingpoints,nscores::Int;addtotal::Bool=false)
	scorebandsUTF8List=Array{String}(0)
	nClasses=length(scorebandsstartingpoints)
	for i=1:nClasses
		i<nClasses ? last=string(scorebandsstartingpoints[i+1]-1) : last=string(nscores)
		#push!(scorebandsUTF8List,string('\'',scorebandsstartingpoints[i],"-",last)) #we add an apostrophe here to avoid autoformatting in Excel (it will interpret 1-99 as a date otherwise)
		push!(scorebandsUTF8List,string("=",DoubleQuote,scorebandsstartingpoints[i],"-",last,DoubleQuote)) #we add an apostrophe here to avoid autoformatting in Excel (it will interpret 1-99 as a date otherwise)
	end
	if addtotal
	 scorebandsUTF8List=vcat(scorebandsUTF8List,["Total"])
	end
	return scorebandsUTF8List
end

function createTrnValStatsForThisIteration(description,iter::Int,scorebandsstartingpoints::Array{Int64,1}
,numeratortrn,denominatortrn,weighttrn,numeratorEsttrn,scorestrn
,numeratorval,denominatorval,weightval,numeratorEstval,scoresval,boolCalculateGini::Bool)
		nClasses=length(scorebandsstartingpoints)		
		sumnumeratortrn,sumdenominatortrn,sumweighttrn,sumnumeratorEstimatetrn=aggregateScores(numeratortrn,denominatortrn,weighttrn,numeratorEsttrn,scorestrn,scorebandsstartingpoints)
		sumnumeratorval,sumdenominatorval,sumweightval,sumnumeratorEstimateval=aggregateScores(numeratorval,denominatorval,weightval,numeratorEstval,scoresval,scorebandsstartingpoints)		
		if boolCalculateGini
			#these statistics are optinally disabled as the gini takes a LOT of time especially due to sorting			
			#NOTE: we currently have two different gini values which we calculate, one approach takes only one argument (the estimator) while the other appraoch considers the truth and the estimator
			#it seems that the single_agrument_gini is more common (it should correspond to the gini in library(reldist) in R).
			gini_single_argtrn=gini_single_argument(numeratorEsttrn.*denominatortrn)
			gini_single_argval=gini_single_argument(numeratorEstval.*denominatorval)			
			normalized_unweighted_gini_numeratorTrn=normalized_gini(numeratortrn,numeratorEsttrn) 
			normalized_unweighted_gini_numeratorVal=normalized_gini(numeratorval,numeratorEstval)
			total_sum_squares_of_numeratorTrn,residual_sum_squares_of_numeratorTrn,total_sum_squares_of_ratioTrn,residual_sum_squares_of_ratioTrn=calc_sum_squares(numeratortrn,numeratorEsttrn,denominatortrn)
			total_sum_squares_of_numeratorVal,residual_sum_squares_of_numeratorVal,total_sum_squares_of_ratioVal,residual_sum_squares_of_ratioVal=calc_sum_squares(numeratorval,numeratorEstval,denominatorval)		
		else
			gini_single_argtrn=1.0
			gini_single_argval=1.0
			normalized_unweighted_gini_numeratorTrn=1.0 
			normalized_unweighted_gini_numeratorVal=1.0 	
			total_sum_squares_of_numeratorTrn=residual_sum_squares_of_numeratorTrn=total_sum_squares_of_ratioTrn=residual_sum_squares_of_ratioTrn=1.0
			total_sum_squares_of_numeratorVal=residual_sum_squares_of_numeratorVal=total_sum_squares_of_ratioVal=residual_sum_squares_of_ratioVal=1.0
		end
		r2_of_numeratortrn=1.0-residual_sum_squares_of_numeratorTrn/total_sum_squares_of_numeratorTrn
		r2_of_ratiotrn=1.0-residual_sum_squares_of_ratioTrn/total_sum_squares_of_ratioTrn
		r2_of_numeratorval=1.0-residual_sum_squares_of_numeratorVal/total_sum_squares_of_numeratorVal
		r2_of_ratioval=1.0-residual_sum_squares_of_ratioVal/total_sum_squares_of_ratioVal

		qwkappaTrn=tryQWKappa(numeratortrn,numeratorEsttrn)
		qwkappaVal=tryQWKappa(numeratorval,numeratorEstval)
		observedratiotrn,observedratioval,fittedratiotrn,fittedratioval,relativitytrn,relativityval,lifttrn,liftval,reversalstrn,reversalsval,relDiffTrnObsVSValObs,relDiffValObsVSTrnFitted, stddevtrn,stddevval,stddevratio,correlation=buildStatisticsInternal(sumnumeratortrn,sumdenominatortrn,sumweighttrn,sumnumeratorEstimatetrn,sumnumeratorval,sumdenominatorval,sumweightval,sumnumeratorEstimateval)
		cumulativeStatsMatrix=[description sumweighttrn sumweightval sumnumeratortrn sumnumeratorval sumdenominatortrn sumdenominatorval observedratiotrn observedratioval sumnumeratorEstimatetrn sumnumeratorEstimateval fittedratiotrn fittedratioval relativitytrn relativityval]
		reltrntitle="Relativity Trn (Observed)" #this is used for the charts
		header=["Cumulative Stats n=$(iter)" "Weight Trn" "Weight Val" "Numerator Trn" "Numerator Val" "Denominator Trn" "Denominator Val" "Observed Ratio Trn" "Observed Ratio Val" "Numerator Est Trn" "Numerator Est Val" "Fitted Ratio Trn" "Fitted Ratio Val" reltrntitle "Relativity Val (Observed)"]
		columnOfRelativityTrn=findin(header,[reltrntitle])[1]
		cumulativeStatsMatrix=vcat(header,cumulativeStatsMatrix)
		singleRowWithKeyMetrics=[iter correlation stddevtrn stddevval stddevratio lifttrn liftval reversalstrn reversalsval relDiffTrnObsVSValObs relDiffValObsVSTrnFitted minimum(relativitytrn) maximum(relativitytrn) minimum(observedratiotrn) maximum(observedratiotrn) minimum(relativityval) maximum(relativityval) minimum(observedratioval) maximum(observedratioval) normalized_unweighted_gini_numeratorTrn normalized_unweighted_gini_numeratorVal gini_single_argtrn gini_single_argval r2_of_ratiotrn r2_of_ratioval r2_of_numeratortrn r2_of_numeratorval residual_sum_squares_of_ratioTrn residual_sum_squares_of_ratioVal residual_sum_squares_of_numeratorTrn residual_sum_squares_of_numeratorVal]
return cumulativeStatsMatrix,singleRowWithKeyMetrics,columnOfRelativityTrn
end

function calc_sum_squares(num_actual,num_estimate,denom)
	@assert length(num_estimate)==length(num_actual)==length(denom)
	sstot=0.0;ssres=0.0;sstot_ratio=0.0;ssres_ratio=0.0
	sz=length(num_actual)
	mean_act=mean(num_actual)
	mean_ratio_pointwise=0.0	
	for i=1:sz
		@inbounds mean_ratio_pointwise += ifelse(denom[i]==0.0,0.0,num_actual[i]/denom[i])
	end 	
	
	for i=1:sz
		@inbounds act=num_actual[i]
		@inbounds den=denom[i]
		#(for boosting) WE NEED TO MULTIPLY THE NUMERATOR WITH THE DENOMINATOR HERE \ TBD, TODO: SOMEONE NEEDS TO REVIEW THIS
		@inbounds est=num_estimate[i]*den
		ssres+=abs2(est-act)
		sstot+=abs2(act-mean_act) #THIS SHOULD BE CACHED (be careful as the training data set may change over time)! todo, tbd

		ssres_ratio+=ifelse(den==0,0.0,abs2((est-act)/den))
		sstot_ratio+=ifelse(den==0,0.0,abs2(act/den-mean_ratio_pointwise))		
		#if i<3
		#	@show mean(num_actual),mean(num_estimate),mean(denom),sum(num_actual),sum(num_estimate),sum(denom)
		#	@show i,act,est,abs2(est-act),abs2(act-mean_act),abs2((est-act)/den),abs2(act/den-mean_ratio_pointwise)
		#end
	end
		
	#@show sstot,ssres,sstot_ratio,ssres_ratio
	return sstot,ssres,sstot_ratio,ssres_ratio
end

#createSingleTreeExcel(trnidx,validx,sett,tree,leaves_of_tree,fitted_values_tree,leaf_number_vector,numerator,denominator,weight)
function createSingleTreeExcel(trnidx::Vector{Int},validx::Vector{Int},sett::ModelSettings,tree::Tree,leaves_of_tree::Vector{Leaf},est::Vector{Float64},leafNumbers::Array{Int64,1},numerator::Array{Float64,1},denominator::Array{Float64,1},weight::Array{Float64,1})
	nameOflistofpredictorsSheet,nameOfpredictorsSheet,nameOfModelStatisticsSheet,nameOfScoresSheet,nameOfOverviewSheet,nameOfSettingsSheet,nameOfValidationSheet,xlData=initExcelData()
	overallstats=createTrnValStats!(trnidx,validx,sett,nameOfModelStatisticsSheet,nameOfSettingsSheet,xlData,tree,leaves_of_tree,est,leafNumbers,numerator,denominator,weight)
return xlData,overallstats
end

#,leaves_of_tree::Vector{Leaf},est::Vector{Float64},leafNumbers::Array{Int64,1},numerator::Array{Float64,1},denominator::Array{Float64,1},weight::Array{Float64,1})
function createTrnValStats!(trnidx::Vector{Int},validx::Vector{Int},sett::ModelSettings,nameOfModelStatisticsSheet::T,nameOfSettingsSheet::T,xlData,tree::Tree,leaves_of_tree::Vector{Leaf},est::Vector{Float64},leafNumbers::Array{Int64,1},numerator::Array{Float64,1},denominator::Array{Float64,1},weight::Array{Float64,1}) where {T <: AbstractString}
t=tree.rootnode
#function for SINGLE Tree
	nClasses=length(leaves_of_tree)
	#srt=sortperm(fittedValuesPerLeaf) #srt[1] is the leaf number (canonical leaf numbering according to concatenation of leaves -> see create_leaves_of_tree) with the lowest fitted value
	leafnrlist=Int[x.id for x in leaves_of_tree]
	#@info "it is important that this list is constructed consistently with the definition of the leaf numbers"
	srt=sortperm(leafnrlist)

	#trn data
	sumnumeratortrn = zeros(Float64, nClasses)
	sumdenominatortrn = zeros(Float64, nClasses)
	sumweighttrn = zeros(Float64, nClasses)
	sumnumeratorEstimatetrn = zeros(Float64, nClasses)
	for i=1:nClasses
		thisleaf=leaves_of_tree[srt[i]]
		sumnumeratortrn[i]+=thisleaf.sumnumerator
		sumdenominatortrn[i]+=thisleaf.sumdenominator
		sumweighttrn[i]+=thisleaf.size
		sumnumeratorEstimatetrn[i]+=thisleaf.fitted*thisleaf.sumdenominator #should be denominator..... nodesize(thisleaf)	#this was rowcount instead of nodesize in an earlier version
	end

	#val data
	sumnumeratorval = zeros(Float64, nClasses)
	sumdenominatorval = zeros(Float64, nClasses)
	sumweightval = zeros(Float64, nClasses)
	sumnumeratorEstimateval = zeros(Float64, nClasses)
	#for i=1:nClasses
	@inbounds for i in validx
		this_is_the_leaf_no=leafNumbers[i]
		sumnumeratorval[this_is_the_leaf_no]+=numerator[i] 
		sumdenominatorval[this_is_the_leaf_no]+=denominator[i] 
		sumweightval[this_is_the_leaf_no]+=weight[i]
		sumnumeratorEstimateval[this_is_the_leaf_no]+=denominator[i]*est[i]
	end
	segmentlist=collect(1:length(leaves_of_tree))
	segmentlist=[segmentlist; "Total"]

	header=["Segment" "Weight" "Relative Weight" "Numerator" "Denominator" "Observed Ratio" "Numerator Estimate" "Fitted Ratio" "Relativity"]
	thisres,overallstats=buildStatistics(header,segmentlist,sumnumeratortrn,sumdenominatortrn,sumweighttrn,sumnumeratorEstimatetrn,sumnumeratorval,sumdenominatorval,sumweightval,sumnumeratorEstimateval)
	if sett.boolCalculateGini
        nest=denominator.*est
        giniTrn=gini_single_argument(view(nest,trnidx))			
        giniVal=gini_single_argument(view(nest,validx))        
        #normalized_unweighted_gini_numeratorTrn=normalized_gini(numeratortrn,numeratorEsttrn) 
        #normalized_unweighted_gini_numeratorVal=normalized_gini(numeratorval,numeratorEstval)
        total_sum_squares_of_numeratorTrn,residual_sum_squares_of_numeratorTrn,total_sum_squares_of_ratioTrn,residual_sum_squares_of_ratioTrn=calc_sum_squares(view(numerator,trnidx),view(est,trnidx),view(denominator,trnidx))
        total_sum_squares_of_numeratorVal,residual_sum_squares_of_numeratorVal,total_sum_squares_of_ratioVal,residual_sum_squares_of_ratioVal=calc_sum_squares(view(numerator,validx),view(est,validx),view(denominator,validx))
    else 
    #gini
		giniTrn=1.0 #normalized_gini(view(numerator,trnidx),view(est,trnidx))  #currently disabled gini takes a LOT of time especially due to sorting
		giniVal=1.0 #normalized_gini(view(numerator,validx),view(est,validx))        	
        total_sum_squares_of_numeratorTrn=residual_sum_squares_of_numeratorTrn=total_sum_squares_of_ratioTrn=residual_sum_squares_of_ratioTrn=1.0
        total_sum_squares_of_numeratorVal=residual_sum_squares_of_numeratorVal=total_sum_squares_of_ratioVal=residual_sum_squares_of_ratioVal=1.0
    end
    warn("need to review/amend the rsquared and rss which are not correct!")
    r2_of_numeratortrn=1.0-residual_sum_squares_of_numeratorTrn/total_sum_squares_of_numeratorTrn
    r2_of_ratiotrn=1.0-residual_sum_squares_of_ratioTrn/total_sum_squares_of_ratioTrn
    r2_of_numeratorval=1.0-residual_sum_squares_of_numeratorVal/total_sum_squares_of_numeratorVal
    r2_of_ratioval=1.0-residual_sum_squares_of_ratioVal/total_sum_squares_of_ratioVal    
    #attach Gini
		ginirowtrn=["Gini Numerator Trn" giniTrn repmat([""],1,size(thisres,2)-2)]
		thisres=[thisres;ginirowtrn];overallstats=[overallstats;ginirowtrn];
		ginirowval=["Gini Numerator Val" giniVal repmat([""],1,size(thisres,2)-2)]
		thisres=[thisres;ginirowval];overallstats=[overallstats;ginirowval];
	#quadratic weighted kappa
		qwkappaTrn=tryQWKappa(view(numerator,trnidx),view(est,trnidx))
		qwkappaVal=tryQWKappa(view(numerator,validx),view(est,validx))
	#attach qwkappa
		qwkapparowtrn=["Quadratic Weighted Kappa (Numerator) Trn" qwkappaTrn repmat([""],1,size(thisres,2)-2)]
		thisres=[thisres;qwkapparowtrn];overallstats=[overallstats;qwkapparowtrn];
		qwkapparowval=["Quadratic Weighted Kappa (Numerator) Val" qwkappaVal repmat([""],1,size(thisres,2)-2)]
		thisres=[thisres;qwkapparowval];overallstats=[overallstats;qwkapparowval];
    #attach r squared
        r2rowtrn=["R Squared of Numerator Trn" r2_of_numeratortrn repmat([""],1,size(thisres,2)-2)]
		thisres=[thisres;r2rowtrn];overallstats=[overallstats;r2rowtrn];
		r2rowval=["R Squared of Numerator Val" r2_of_numeratorval repmat([""],1,size(thisres,2)-2)]
		thisres=[thisres;r2rowval];overallstats=[overallstats;r2rowval];
    #attach rss
        rss_rowtrn=["RSS of Numerator Trn" residual_sum_squares_of_numeratorTrn repmat([""],1,size(thisres,2)-2)]
		thisres=[thisres;rss_rowtrn];overallstats=[overallstats;rss_rowtrn];
		rss_rowval=["RSS of Numerator Val" residual_sum_squares_of_numeratorVal repmat([""],1,size(thisres,2)-2)]
		thisres=[thisres;rss_rowval];overallstats=[overallstats;rss_rowval];
    
	modelStatisticsSheet=ExcelSheet(nameOfModelStatisticsSheet,thisres)
	modelsettingsSheet=ExcelSheet(nameOfSettingsSheet,convert(DataFrame,writeAllFieldsToArray(sett)))
	xlData.sheets=[modelsettingsSheet,modelStatisticsSheet]
	#add charts
		thischart=defineRelativityChart(nameOfModelStatisticsSheet,nameOfModelStatisticsSheet,"K2",size(leaves_of_tree,1),1,1,headerrow2=5+size(leaves_of_tree,1)+1,headercol2=1,datarow2=5+size(leaves_of_tree,1)+2+1,xtitle="Segment",ytitle="Observed Ratio",xscale=2.4,title="Observed Ratio per Segment",datarow=3,valuescol=6,categoriescol=1,valuescol2=6,yscale=1.5)
		push!(xlData.charts,deepcopy(thischart))
		thischart=defineRelativityChart(nameOfModelStatisticsSheet,nameOfModelStatisticsSheet,string("K",5+size(leaves_of_tree,1)+2),size(leaves_of_tree,1),1,1,headerrow2=5+size(leaves_of_tree,1)+1,headercol2=1,datarow2=5+size(leaves_of_tree,1)+2+1,xtitle="Segment",ytitle="Relative Weight per Segment",xscale=2.4,title="Relative Weight per Segment",datarow=3,valuescol=3,categoriescol=1,valuescol2=3,yscale=1.5)
		push!(xlData.charts,deepcopy(thischart))
	return overallstats[:,1:2]
	#return thisres
end

function pushsum!(x::Array)
	push!(x,sum(x))
end

function wAverage(vec,w)
	@assert size(vec)==size(w)
	res=zero(eltype(vec))
	wsum=zero(eltype(vec))
	@inbounds for i=1:size(vec,1)
		res+=vec[i]*w[i]
		wsum+=w[i]
	end
	return res/wsum
end
function wAverage(vec::Array{T,1},w::Array{U,1}) where {T,U}
	@assert size(vec)==size(w)
	res=zero(eltype(vec))
	wsum=zero(eltype(vec))
	@inbounds for i=1:size(vec,1)
		res+=vec[i]*w[i]
		wsum+=w[i]
	end
	return res/wsum
end

function buildStatisticsInternal(sumnumeratortrn,sumdenominatortrn,sumweighttrn,sumnumeratorEstimatetrn,sumnumeratorval,sumdenominatorval,sumweightval,sumnumeratorEstimateval)
	relativitytrn=(sumnumeratortrn./sumdenominatortrn)/(sum(sumnumeratortrn)/sum(sumdenominatortrn))
	relativityval=(sumnumeratorval./sumdenominatorval)/(sum(sumnumeratorval)/sum(sumdenominatorval)) #todo/tbd this can be improved: sum(sumnumeratorval)/sum(sumdenominatorval) will never change through the full boosting model, thus we do not need to calculate it all the time
	correlation=cor(relativitytrn,relativityval)
	#this is not meaningful if we have reversals!
    #mi,ma=extrema(relativitytrn);lifttrn=ma/mi
    #mi,ma=extrema(relativityval);liftval=ma/mi
    lifttrn=relativitytrn[end]/relativitytrn[1]
    liftval=relativityval[end]/relativityval[1]
	reversalstrn=sum(reversals(relativitytrn))	
	reversalsval=sum(reversals(relativityval))

	#Attach relativity of one for Total
	relativitytrn=push!(relativitytrn,1.0)
	relativityval=push!(relativityval,1.0)
	#Calculate Totals
	pushsum!(sumweighttrn)
	pushsum!(sumnumeratortrn)
	pushsum!(sumdenominatortrn)
	pushsum!(sumnumeratorEstimatetrn)
	pushsum!(sumweightval)
	pushsum!(sumnumeratorval)
	pushsum!(sumdenominatorval)
	pushsum!(sumnumeratorEstimateval)
	#calculate ratios
	observedratiotrn=sumnumeratortrn./sumdenominatortrn
	fittedratiotrn=sumnumeratorEstimatetrn./sumdenominatortrn
	observedratioval=sumnumeratorval./sumdenominatorval
	fittedratioval=sumnumeratorEstimateval./sumdenominatorval

	#Calculate metrics
	reldiff=relativitytrn.-relativityval
	sz_1=length(sumweighttrn)-1
	@assert length(sumweighttrn)==length(sumweightval)==length(fittedratioval)==length(fittedratiotrn)==length(reldiff)==length(sumnumeratorEstimatetrn)
	stddevtrn=sum(view(sumweighttrn,1:sz_1).*(view(fittedratiotrn,1:sz_1).-fittedratiotrn[end]).^2)/sumweighttrn[end]
    stddevval=sum(view(sumweightval,1:sz_1).*(view(fittedratioval,1:sz_1).-fittedratioval[end]).^2)/sumweightval[end]
	stddevratio=stddevval/stddevtrn
    
    relDiffTrnObsVSValObs=wAverage(abs.(view(reldiff,1:sz_1)),view(sumnumeratorEstimatetrn,1:sz_1))
	reldiff2=relativityval.-fittedratiotrn
	relDiffValObsVSTrnFitted=wAverage(abs.(view(reldiff2,1:sz_1)),view(sumnumeratorEstimatetrn,1:sz_1))    
return observedratiotrn,observedratioval,fittedratiotrn,fittedratioval,relativitytrn,relativityval,lifttrn,liftval,reversalstrn,reversalsval,relDiffTrnObsVSValObs,relDiffValObsVSTrnFitted, stddevtrn,stddevval,stddevratio,correlation
end

function buildStatistics(header,description,sumnumeratortrn,sumdenominatortrn,sumweighttrn,sumnumeratorEstimatetrn,sumnumeratorval,sumdenominatorval,sumweightval,sumnumeratorEstimateval)
	#Correlation of relativities
    observedratiotrn,observedratioval,fittedratiotrn,fittedratioval,relativitytrn,relativityval,lifttrn,liftval,reversalstrn,reversalsval,relDiffTrnObsVSValObs,relDiffValObsVSTrnFitted,stddevtrn,stddevval,stddevratio,correlation=buildStatisticsInternal(sumnumeratortrn,sumdenominatortrn,sumweighttrn,sumnumeratorEstimatetrn,sumnumeratorval,sumdenominatorval,sumweightval,sumnumeratorEstimateval)

	statstrn=[description sumweighttrn float(sumweighttrn)/sumweighttrn[end] sumnumeratortrn sumdenominatortrn observedratiotrn sumnumeratorEstimatetrn fittedratiotrn relativitytrn]
	statsval=[description sumweightval float(sumweightval)/sumweightval[end] sumnumeratorval sumdenominatorval observedratioval sumnumeratorEstimateval fittedratioval relativityval]

	#add a header row and two empty rows in between
	headerval=["Validation Data" repmat([""],1,length(header)-1)]
	headertrn=deepcopy(headerval)
	headertrn[1]="Training Data"
	empty=repmat([""],1,length(header))

	corrow=["Correlation" correlation repmat([""],1,length(header)-2)]
	stddevrows=[["Std Dev Trn" stddevtrn];["Std Dev Val" stddevval];["Std Dev Ratio" stddevratio]]
	stddevrows=[stddevrows repmat([""],3,length(header)-2)]
	#lifttrn,liftval,reversalstrn,reversalsval,relDiffTrnObsVSValObs,relDiffValObsVSTrnFitted
	lifttrn=["Lift Trn" lifttrn repmat([""],1,length(header)-2)]
	liftval=["Lift Val" liftval repmat([""],1,length(header)-2)]
	reversalstrn=["Reversals Trn" reversalstrn repmat([""],1,length(header)-2)]
	reversalsval=["Reversals Val" reversalsval repmat([""],1,length(header)-2)]
	relDiffTrnObsVSValObs=["RelDiff TrnObs vs ValObs" relDiffTrnObsVSValObs repmat([""],1,length(header)-2)]
	relDiffValObsVSTrnFitted=["RelDiff ValObs vs TrnFitted" relDiffValObsVSTrnFitted repmat([""],1,length(header)-2)]

	stats=[headertrn;header;statstrn;empty;empty;headerval;header;statsval;empty;empty]
	overallstats=[corrow;stddevrows;lifttrn;liftval;reversalstrn;reversalsval;relDiffTrnObsVSValObs;relDiffValObsVSTrnFitted]
	stats=vcat(stats,overallstats)
	return stats,overallstats
end

#todo/tbd find out how to parametrize these functions properly....
#this syntax does not seem to work when it is called with a labellist::Array{UInt8,1}
#function sortlists!(catSortBy::SortByMean,labellist::Array{T,1},sumnumerator::Array{Float64,1},sumdenominator::Array{Float64,1},sumweight::Array{Float64,1},countlistfloat::Array{Float64,1}) where T<:UInt
#ERROR: MethodError: no method matching sortlists!(::DecisionTrees.SortByMean, ::Array{UInt8,1}, ::Array{Float64,1}, ::Array{Float64,1}, ::Array{Float64,1}, ::Array{Float64,1})
#Closest candidates are:
#sortlists!(::DecisionTrees.SortByMean, ::Array{T<:UInt64,1}, ::Array{Float64,1}, ::Array{Float64,1}, ::Array{Float64,1}, ::Array{Float64,1}) where T<:UInt64 at C:\Users\bernhard.konig\Documents\ASync\home\Code\Julia\DecisionTrees.jl\src\helper_functions.jl:2060
#sortlists!(::DecisionTrees.SortByMean, ::Array{T<:UInt64,1}, ::Array{Float64,1}, ::Array{Float64,1}, ::Array{Float64,1}, ::Array{Float64,1}, ::Array{Float64,1}) where T<:UInt64 at C:\Users\bernhard.konig\Documents\ASync\home\Code\Julia\DecisionTrees.jl\src\helper_functions.jl:2085
#sortlists!(::DecisionTrees.SortByMean, ::Array{T<:UInt64,1}, ::Array{Int64,1}, ::Array{Float64,1}, ::Array{Float64,1}) where T<:UInt64 at C:\Users\bernhard.konig\Documents\ASync\home\Code\Julia\DecisionTrees.jl\src\helper_functions.jl:2073
function sortlists!(catSortBy::SortByMean,labellist::Array{T,1},sumnumerator::Array{Float64,1},sumdenominator::Array{Float64,1},sumweight::Array{Float64,1},countlistfloat::Array{Float64,1}) where T<:Unsigned
	#todo/tbd generalize this function for an arbitrary number of input vectors
		#function for sorting all the lists by the mean
		meanratioperclass=sumnumerator./sumdenominator
		srt=sortperm(meanratioperclass,alg=QuickSort)
		myresort!(labellist,srt)
		myresort!(sumnumerator,srt)
		myresort!(sumdenominator,srt)
		myresort!(sumweight,srt)
		myresort!(countlistfloat,srt)
		return nothing
end

function sortlists!(catSortBy::SortByMean,labellist::Array{T,1},intSumrklost::Array{Int64,1},sumweight::Array{Float64,1},countlistfloat::Array{Float64,1}) where T<:Unsigned
	#function for sorting all the lists by the mean
#	meanratioperclass=sumnumerator./sumdenominator
	#srt=sortperm(meanratioperclass,alg=QuickSort)
	srt=sortperm(intSumrklost,alg=QuickSort)	 #todo/tbd check if this makes sense sorting by #rklost....?
	myresort!(labellist,srt)
	myresort!(intSumrklost,srt)
	#myresort!(sumdenominator,srt)
	myresort!(sumweight,srt)
	myresort!(countlistfloat,srt)
	return nothing
end

function sortlists!(catSortBy::SortByMean,labellist::Array{T,1},sumnumerator::Array{Float64,1},sumdenominator::Array{Float64,1},sumweight::Array{Float64,1},countlistfloat::Array{Float64,1},moments_per_pdaclass::Array{Float64,1}) where T<:Unsigned
	#todo/tbd generalize this function for an arbitrary number of input vectors
		#function for sorting all the lists by the mean
		meanratioperclass=sumnumerator./sumdenominator
		srt=sortperm(meanratioperclass,alg=QuickSort)
		myresort!(labellist,srt)
		myresort!(sumnumerator,srt)
		myresort!(sumdenominator,srt)
		myresort!(sumweight,srt)
		myresort!(countlistfloat,srt)
		myresort!(moments_per_pdaclass,srt)
		return nothing
end

function myresort!(x,srt)
  x[:]=x[srt]
  #z=copy(x)[srt]
  #for i=1:length(x)
  #   x[i]=z[i]
  #end
return nothing
end

function print_graycode(m)
  #this functions print the first m gray codes (binary numbers such that with each step exactly one bit flips)
  n=uint8(m)
  current=zero(n)
  o=one(n)
  for i=zero(n):n
    @show Int64(i),bits(current)
    t=trailing_ones(i)
    x=o<<t
    current=current$x
  end
#usage: print_graycode(15)
end

function crit_string_to_type(str)
  if str=="_difference"
    return DifferenceSplit()
  end
  if str=="_mse"
    return MSESplit()
  end
  if str=="_RankOpt"
    return RankOptSplit()
  end
  @assert false "Unknown split criterion $(str). ABORT."
end

function countsort(a::Array{T,1}) where {T <: Integer}
    (lo, hi) = extrema(a)
    b = zeros(T, length(a))
    cnt = zeros(Int, hi - lo + 1)
    for i in a
        cnt[i - lo + 1] += 1
    end
    z = one(Int)
    for i in lo:hi
        while cnt[i - lo + 1] > 0
            b[z] = i
            z += 1
            cnt[i - lo + 1] -= 1
        end
    end
    return b
end

function countsort!(a::Array{T,1}) where {T <: Integer}
    (lo, hi) = extrema(a)
#     b = zeros(T, length(a))
    cnt = zeros(Int, hi - lo + 1)
    for i in a
     cnt[i - lo + 1] += 1
    end

    z = one(Int)
    for i in lo:hi
        while cnt[i - lo + 1] > 0
            a[z] = i
            z += 1
            cnt[i - lo + 1] -= 1
        end
    end
end

function get_firstpos(v::Array{T,1}) where {T}
  #v=volume vector (v[i] = number of occurrences of element k in the feature vector)
  #return value is the a vector containing the first positions of element k in an ordered set
  firstpos=zeros(Int,length(v))
  firstpos[1]=1
  #@show "gaht nöd...."
  for i=2:length(v)
      firstpos[i]=firstpos[i-1]+v[i-1]
  end
  return firstpos
end

function aggregate_data(f::pdaMod,scores::Array{Int64,1},numeratorEst::Array{Float64,1},numerator::Array{Float64,1},denominator::Array{Float64,1},weight::Array{Float64,1})
    #this is the core function of the modelling process
	#besides copying of the data, the vast majority of time is spent in here!
	#most of the time is spent here, if we can improve the for loop below, that would improve performance greatly!
	#one possibility would be to introduce parallelization here (which is not straightforward, I think....)
	a= f.pda.refs
	#this should hold by definition/construction
	(lo, hi) = extrema(levels(f))
	ooo=one(lo)-lo
	vecsize=hi+ooo
	@assert length(f.ids)<=vecsize #need to investigate if it can happen that this does not hold true (tod/tbd)
	if length(f.values)!=vecsize #need to investigate if it can happen that this does not hold true (tod/tbd)
		rows_to_be_deleted=myfindinInt(collect(lo:hi),f.ids)
	end
	valuelist=deepcopy(f.values)
	#@assert vecsize==length(f.pda.
    cnt = zeros(Int, vecsize)
    sumnumerator = zeros(Float64, vecsize)
	sumnumeratorEst = zeros(Float64, vecsize)
	sumdenominator = zeros(Float64, vecsize)
	sumweight= zeros(Float64, vecsize)
	sumscores=zeros(eltype(scores),vecsize)
	#note: inbounds increases efficiency here (about a factor of 2), however if the bounds are violated something nasty might happen (quote: If the subscripts are ever out of bounds, you may suffer crashes or silent corruption.)
	for count=1:length(a)
        #i=a[count] #i in a
		#idx=i - lo + 1
		@inbounds idx=a[count] + ooo
		@inbounds cnt[idx] += 1
		@inbounds sumnumeratorEst[idx] += numeratorEst[count]
		@inbounds sumnumerator[idx] += numerator[count]
		@inbounds sumdenominator[idx] += denominator[count]
		@inbounds sumweight[idx] += weight[count]
		@inbounds sumscores[idx] += scores[count]
    end
	if length(f.values)!=vecsize #delete some rows
		sort!(rows_to_be_deleted,rev=true)
		for i in rows_to_be_deleted
			deleteat!(cnt,i)
			deleteat!(sumscores,i)
			deleteat!(sumnumeratorEst,i)
			deleteat!(sumnumerator,i)
			deleteat!(sumdenominator,i)
			deleteat!(sumweight,i)
		end
	end
  return valuelist,cnt,sumscores,sumnumeratorEst,sumnumerator,sumdenominator,sumweight
end



function aggregate_data(f::PooledArray,scores,numeratorEst,numerator,denominator,weight)
	lo=one(eltype(f.refs))
	hi=convert(eltype(f.refs),length(f.pool)) 	
	ooo=one(lo)-lo
	vecsize=hi+ooo

	valuelist=levels(f) #we should not use levels here, it is the wrong function
    cnt = zeros(Int, vecsize)
    sumnumerator = zeros(Float64, vecsize)
	sumnumeratorEst = zeros(Float64, vecsize)
	sumdenominator = zeros(Float64, vecsize)
	sumweight= zeros(Float64, vecsize)
	sumscores=zeros(eltype(scores),vecsize)		
			  
	@inbounds for count in 1:length(f) #f.indices[1]
		idx=f.refs[count] + ooo
		cnt[idx] += 1    			
		sumnumeratorEst[idx] += numeratorEst[count]
		sumnumerator[idx] += numerator[count]
		sumdenominator[idx] += denominator[count]
		sumweight[idx] += weight[count]
		sumscores[idx] += scores[count]
	end
	
	#is the following faster, if we do rows_to_be_deleted=findin(cnt,0);deleteat!(cnt,rows_to_be_deleted);....
	for i=vecsize:-1:1
		@inbounds if cnt[i]==0
	#if length(f.values)!=vecsize #delete some rows
			#warn("BK, I am not sure if this works as anticipated")
		#for i in rows_to_be_deleted
			deleteat!(cnt,i)
			deleteat!(sumscores,i)
			deleteat!(sumnumeratorEst,i)
			deleteat!(sumnumerator,i)
			deleteat!(sumdenominator,i)
			deleteat!(sumweight,i)
			deleteat!(valuelist,i)
		end		
	end
  return valuelist,cnt,sumscores,sumnumeratorEst,sumnumerator,sumdenominator,sumweight
end

function aggregate_data_diff(f::T,numerator::Array{Float64,1},denominator::Array{Float64,1},weight::Array{Float64,1})  where T<:SubArray
	#! preliminary testing showed that we can save 10% if we explicitly type f as in
	# aggregate_data_diff(f::SubArray{String,1,PooledArray{String,UInt8,1},Tuple{Array{Int64,1}},false},denominator....)	
    #if we never subset the pda and always work on a view of the original data, then lo, hi are always 1: size(pool)!
    lo=one(eltype(f.parent.refs)) #UInt8(1)    
    hi=convert(eltype(f.parent.refs),length(f.parent.pool)) #hi=convert(UInt8,length(f.parent.pool))

	ooo=one(lo)-lo
	vecsize=hi+ooo
    cnt = zeros(Int, vecsize)
    sumnumerator = zeros(Float64, vecsize)
	sumdenominator = zeros(Float64, vecsize)
	sumweight= zeros(Float64, vecsize)

	#warn("BK: need to add inbounds here as soon as things work.... ")
	 @inbounds for count in f.indices[1]
	#for count=1:length(a)
		#@inbounds idx=a[count] + ooo
		idx=f.parent.refs[count] + ooo
		cnt[idx] += 1
		sumnumerator[idx] += numerator[count]
		sumdenominator[idx] += denominator[count]
		sumweight[idx] += weight[count]
	end
  return cnt,sumnumerator,sumdenominator,sumweight
end

function myfindin(big::Array{T,1},small::Array{T,1}) where {T}
	#returns boolean array which indicates which elements are not found in the smaller array
	res=Array{Bool}(length(big))
	for i=1:length(big)
		z=big[i]
		res[i]=in(z,small)
	end
	return res
end

function myfindinInt(big::Array{T,1},small::Array{T,1}) where {T}
	#returns integer array which indicates which elements are not found in the smaller array
	res=Array{Int}(0)
	for i=1:length(big)
		z=big[i]
		if !in(z,small)
			push!(res,i)
		end
	end
	return res
end

function aggregate_data_rklostcount(f::pdaMod,increasedPremiumVector::Array{Float64,1},intRanksLostVector::Array{Int64,1},numerator::Array{Float64,1},denominator::Array{Float64,1},weight::Array{Float64,1})
    a=f.pda.refs
    (lo, hi) = extrema(levels(f))
    lengtha=length(a)
    cnt = zeros(Int, hi - lo + 1)
	intSumrklost = zeros(Int64, hi - lo + 1)
	sumweight= zeros(Float64, hi - lo + 1)
	ooo=one(lo)-lo
    #note: inbounds increases efficiency here (about a factor of 2), however if the bounds are violated something nasty might happen (quote: If the subscripts are ever out of bounds, you may suffer crashes or silent corruption.)
	for count=1:lengtha
        @inbounds idx=a[count] + ooo
        @inbounds cnt[idx] += 1
		@inbounds intSumrklost[idx] += intRanksLostVector[count]
		@inbounds sumweight[idx] += weight[count]
    end
  return cnt,intSumrklost,sumweight
end

function build_listOfMeanResponse(crit::MaxValueSplit,numerator::Array{Float64,1},denominator::Array{Float64,1},weight::Array{Float64,1},features::pdaMod,feature_levels::Array{UInt8,1},minweight::Float64)
  #todo/tbd if we sort the features anyway here, then, we can determine "feature_levels" more efficiently after the sorting (without using the levels function)
  ncategories=length(feature_levels)
  #countlist,sumobservedlist=aggregate_data_diff(features,numerator,denominator,weight)
  countlist,sumnumerator,sumdenominator,sumweight=aggregate_data_diff(features,numerator,denominator,weight)
  if length(countlist)!=ncategories
  #there are gaps with no data between min(features.ids),max(features.ids)
     zeroidx=findin(countlist,0)	 #@show length(zeroidx)
     deleteat!(countlist,zeroidx)
	 deleteat!(sumnumerator,zeroidx)
	 deleteat!(sumdenominator,zeroidx)
	 deleteat!(sumweight,zeroidx)
  end
  countlistfloat=Float64(countlist)
  #sorting is not done here! if we want the data sorted, is suggest to do this in another function which uses the output from this fn
return feature_levels,sumnumerator,sumdenominator,sumweight,countlistfloat
end

function build_listOfMeanResponse(crit::MaxMinusValueSplit,numerator::Array{Float64,1},denominator::Array{Float64,1},weight::Array{Float64,1},features::pdaMod,feature_levels::Array{UInt8,1},minweight::Float64)
  #todo/tbd if we sort the features anyway here, then, we can determine "feature_levels" more efficiently after the sorting (without using the levels function)
  ncategories=length(feature_levels)
  #countlist,sumobservedlist=aggregate_data_diff(features,numerator,denominator,weight)
  countlist,sumnumerator,sumdenominator,sumweight=aggregate_data_diff(features,numerator,denominator,weight)
  if length(countlist)!=ncategories
  #there are gaps with no data between min(features.ids),max(features.ids)
     zeroidx=findin(countlist,0)	 #@show length(zeroidx)
     deleteat!(countlist,zeroidx)
	 deleteat!(sumnumerator,zeroidx)
	 deleteat!(sumdenominator,zeroidx)
	 deleteat!(sumweight,zeroidx)
  end
  countlistfloat=Float64(countlist)
  #sorting is not done here! if we want the data sorted, is suggest to do this in another function which uses the output from this fn
return feature_levels,sumnumerator,sumdenominator,sumweight,countlistfloat
end


function build_listOfMeanResponse(crit::NormalDevianceSplit,numerator::Array{Float64,1},denominator::Array{Float64,1},weight::Array{Float64,1},features::pdaMod,feature_levels::Array{UInt8,1},minweight::Float64)
  ncategories=length(feature_levels)
  #calc variance for each "isolated" category
  countlist,sumnumerator,sumdenominator,sumweight,moments_per_pdaclass=aggregate_data_normal_deviance(features,numerator,denominator,weight)
  if length(countlist)!=ncategories
  #there are gaps with no data between min(features.ids),max(features.ids)
     zeroidx=findin(countlist,0)	 #@show length(zeroidx)
     deleteat!(countlist,zeroidx)
	 deleteat!(sumnumerator,zeroidx)
	 deleteat!(sumdenominator,zeroidx)
	 deleteat!(sumweight,zeroidx)
	 deleteat!(moments_per_pdaclass,zeroidx) #I think the delelat! fn should be able to handle any kind of vector, thus this should work just fine
  end
  countlistfloat=convert(Vector{Float64},countlist) #Float64[convert(Float64,x) for x in countlist]
  #sorting is not done here! if we want the data sorted, is suggest to do this in another function which uses the output from this fn
return feature_levels,sumnumerator,sumdenominator,sumweight,countlistfloat,moments_per_pdaclass
end

function build_listOfMeanResponse(crit::PoissonDevianceSplit,trnidx::Vector{Int},validx::Vector{Int},numerator::Array{Float64,1},denominator::Array{Float64,1},weight::Array{Float64,1},features,feature_levels,minweight::Float64)
    return build_listOfMeanResponse(DifferenceSplit(),trnidx,validx,numerator,denominator,weight,features,feature_levels,minweight)
end 

function build_listOfMeanResponse(crit::DifferenceSplit,trnidx::Vector{Int},validx::Vector{Int},numerator::Array{Float64,1},denominator::Array{Float64,1},weight::Array{Float64,1},features,feature_levels,minweight::Float64)  
  #todo/tbd if we sort the features anyway here, then, we can determine "feature_levels" more efficiently after the sorting (without using the levels function)
  ncategories=length(feature_levels)  #currently also contains the val feature levels!
  countlist,sumnumerator,sumdenominator,sumweight=aggregate_data_diff(features,numerator,denominator,weight)
  if length(countlist)!=ncategories
  #there are gaps with no data between min(features.ids),max(features.ids)
     zeroidx=findin(countlist,0)	 #@show length(zeroidx)
     deleteat!(countlist,zeroidx)
	 deleteat!(sumnumerator,zeroidx)
	 deleteat!(sumdenominator,zeroidx)
	 deleteat!(sumweight,zeroidx)
  end
  countlistfloat=convert(Vector{Float64},countlist) #loat64[convert(Float64,x) for x in countlist]
  #sorting is not done here! if we want the data sorted, is suggest to do this in another function which uses the output from this fn
return feature_levels,sumnumerator,sumdenominator,sumweight,countlistfloat
end

function build_listOfMeanResponse_mse(labels::Array{Float64,1},features::pdaMod,feature_levels::Array{UInt8,1},minweight::Float64,lcwiowatt::Array{Float64,1})
  #todo/tbd if we sort the features anyway here, then, we can determine "feature_levels" more efficiently after the sorting (without using the levels function)
  error("this is currently not working properly")
  warn("ensure that this does not modify features!!")
  ncategories=length(feature_levels)
  features_sorted=deepcopy(features)#   vol,labels_sorted=custom_countsort!(features_sorted,labels,lcwiowatt)
  countlist,firstindices=custom_countsort!(features_sorted,labels,labels_sorted)
  sum_of_sum_of_squareslist=zeros(Float64,ncategories)
  meanobservedlist=zeros(Float64,ncategories)
#   features_sorted=deepcopy(features) 
  #this determines the order of all lists (output of this function)
#   idx_sortfpda=sortperm(features_sorted.pda,alg=QuickSort)
#   features_sorted.pda.refs=features_sorted.pda.refs[idx_sortfpda]
  labellist_sorted=sort(feature_levels)
#   labels_sorted=deepcopy(labels)[idx_sortfpda]
#   firstindices=Array{Int}(ncategories)
#   for i=1:ncategories
#     firstindices[i]=searchsortedfirst(features_sorted.pda.refs,labellist_sorted[i])
#   end

  kmin=1
  for i=1:ncategories
    kmax= i<ncategories ? firstindices[i+1]-1 : length(labels_sorted)
    viewtmp=view(labels_sorted,kmin:kmax)
	meanobservedlist[i]+=sum(viewtmp)
	tmpsumcount=sum(viewtmp.*viewtmp)
    #countlist[i]=kmax-kmin+1
    sum_of_sum_of_squareslist[i]=tmpsumcount
    kmin=kmax+1
  end
      sumobservedlist=copy(sum_of_sum_of_squareslist)
      countlistfloat=Float64(countlist)
      sum_of_sum_of_squareslist=sum_of_sum_of_squareslist./countlistfloat
      #sorting is not done here! if we want the data sorted, is suggest to do this in another function which uses the output from this fn
    return labellist_sorted,meanobservedlist./countlistfloat,countlistfloat,sumobservedlist
end

function rand_bool_idx(n::Int,p::Float64)
  res=BitArray(n)
  for i=1:n
    rand()<p ? res[i]=true : res[i]=false
  end
  return res
end

function dynamic_increment(residuals::Array{Float64,1},current_mdf::Float64,ranks_lost_pct::Float64,variable_mdf_pct::Float64)
    incr=quantile(residuals,ranks_lost_pct)
    incr*=variable_mdf_pct
  return incr
end

function full_indices(parentidx,left)
  obs=size(parentidx,1)
  resl=copy(parentidx)
  resr=copy(parentidx)
  k=1
  for i=1:obs
    if parentidx[i]
      left[k] ? resr[i]=false : resl[i]=false
      k+=1
    end
  end
  return resl,resr
end

function my_write(f::IOStream,arr::AbstractArray;sep::Char=',')
  write(f,"\r\n")
  for i=1:size(arr,1)
    for j=1:size(arr,2)
      write(f,arr[i,j])
      if !(j==size(arr,2))
        write(f,sep)
      end
    end
    write(f,"\r\n")
  end
  write(f,"\r\n")
  nothing
end

function variable_importance_internal(leaves_array::Array{Leaf,1},namevec::Array{String,1},number_of_num_features::Int)
#number_of_num_features=length(intNumVarsUsed)
	sz=length(namevec)
	nleaves=length(leaves_array)
	onedimcount=zeros(Int64,sz)
	twodimcount=zeros(Int64,div(sz*(sz-1),2))
	seen1dim=Array{Int}(0)
	seen2dim=Array{Int}(0)
	for i=1:nleaves
		resize!(seen1dim,0) #reset the array to the empty array
		resize!(seen2dim,0) #reset the array to the empty array
		depth_of_this_leaf=length(leaves_array[i].rule_path)
		for depth=1:depth_of_this_leaf
            #@show leaves_array[i].rule_path
			this_featid = leaves_array[i].rule_path[depth].featid
			if this_featid>number_of_num_features
				dump(leaves_array[i].rule_path)
				@show i;@show depth;@assert false "This should not have happened: A feature id was invalid...."
			end
			variable1 = this_featid < 0 ? -this_featid+number_of_num_features : this_featid
			if !in(variable1,seen1dim)
				onedimcount[variable1]+=1
				push!(seen1dim,variable1)
			end
			for loopvar=depth+1:depth_of_this_leaf
				this_featid2 = leaves_array[i].rule_path[loopvar].featid
				if this_featid2!=this_featid
					variable2 = this_featid2 < 0 ? -this_featid2+number_of_num_features : this_featid2
					#@show variable1,variable2 
                    if variable1==variable2 #this should not happen here!
                        @show variable1,variable2
                        error("this should not have happend!")
                    end
					#add this combination to the list
					rownumber=twodimIndex(variable1,variable2,sz)
					if !in(rownumber,seen2dim) #for each leaf, a combination i,j can only count once. (this is a somewhat arbitrary choice; one could argue how 2dim variable importance is defined though)
						twodimcount[rownumber]+=1
						push!(seen2dim,rownumber)
					end
				end
			end
		end
	end
	return onedimcount,twodimcount
end

function variable_importance(leaves_array::Array{Leaf,1},namevec::Array{String,1},number_of_num_features::Int)
  sz=length(namevec)
  nleaves=length(leaves_array)
  onedimcount,twodimcount=variable_importance_internal(leaves_array,namevec,number_of_num_features)
  twodimcountfloat=convert(Vector{Float64},twodimcount./sum(twodimcount)) #Float64[convert(Float64,x) for x in twodimcount/sum(twodimcount)]
  res1dim=Array{String}(sz,2)
  res2dim=Array{String}(div(sz*(sz-1),2),3)
  count=0
  for i=1:sz, j=i+1:sz
    count+=1
    res2dim[count,1]=namevec[i]
    res2dim[count,2]=namevec[j]
    res2dim[count,3]=string(twodimcountfloat[count])
  end

  res1dim[:,1]=deepcopy(namevec)
  onedimcountfloat=float(onedimcount/sum(onedimcount))
  for j=1:sz
    res1dim[j,2]=string(Float64(onedimcountfloat[j]))
  end
  dropzero=onedimcountfloat.>0
  res1dim=res1dim[dropzero,:]
  onedimcountfloat=deepcopy(onedimcountfloat[dropzero])
  srt=sortperm(onedimcountfloat,rev=true,alg=QuickSort)

  dropzero=twodimcountfloat.>0
  twodimcountfloatdropz=deepcopy(twodimcountfloat[dropzero])
  res2dim=res2dim[dropzero,:]
  srt2=sortperm(twodimcountfloatdropz,rev=true,alg=QuickSort)
  return res1dim[srt,:],res2dim[srt2,:],onedimcount,twodimcount
end

function twodimIndex(variable1::Int,variable2::Int,nvars::Int)
	#calculates the (canonical) index in a two dimensional table i,j (both going from 1 to nvars) -> rownumber of combination i,j
	if variable1<variable2
		v1=variable1
		v2=variable2
	else
		if variable1==variable2
			return nothing
		end
		v1=variable2
		v2=variable1
	end
	idx=0
	for i=1:v1-1
		idx+=(nvars-i)
	end
	idx+=v2-v1
	#@show v1,v2,idx
	#twodimcount[idx]+=1
	return idx
end

function variable_importance_OLD(leaves_array::Array{Leaf,1},namevec::Array{String,1},number_of_num_features::Int64)
  #tbd/todo improve the performance of this code; better yet: construct the variable_importance during the tree building process!
  sz=length(namevec)
  nleaves=length(leaves_array)
  depth=0
  for i=1:nleaves
    depth=max(depth,length(leaves_array[i].rule_path))
  end
  list_of_featureids_on_path_to_leaves=Array{Int64}(nleaves,depth)
  fill!(list_of_featureids_on_path_to_leaves,0)
    for i=1:nleaves
      for k=1:length(leaves_array[i].rule_path)
        this_featid=leaves_array[i].rule_path[k].featid
        if this_featid<0
          this_featid=-this_featid+number_of_num_features
        end
        list_of_featureids_on_path_to_leaves[i,k]=this_featid
      end
  end

  onedimcount=zeros(Int64,sz)
  for j=1:sz, k=1:nleaves
   if in(j,view(list_of_featureids_on_path_to_leaves,k,:))
      onedimcount[j]+=1
   end
  end

  twodimcount=Array{Int64}(div(sz*(sz-1),2))
  fill!(twodimcount,0)
  #todo/tbd twodimcount and onedimcount can be created much more efficiently
  #we do not need to loop through all the variables, but only through the ones which are actually used by the tree!
  #"for i=1:used_by_tree, j=1:used_by_tree; if i!=j; i&j used in conjunction? ;end;end;"
  count=0
  #@info "check if this works as intended!"
  for i=1:sz, j=i+1:sz
      count+=1
	  if min(onedimcount[j],onedimcount[i])>0 #both variables need to be used in the tree
		  for k=1:nleaves
		   this_tmp=view(list_of_featureids_on_path_to_leaves,k,:) #this line and the next take up quite some time
		   if in(i,this_tmp) && in(j,this_tmp) #both indices appear on the path to leaf k
			twodimcount[count]+=1
		   end
	   end
    end
  end
  twodimcountfloat=Float64(twodimcount/sum(twodimcount))

  res1dim=Array{String}(sz,2)
  res2dim=Array{String}(div(sz*(sz-1),2),3)
  count=0
  for i=1:sz, j=i+1:sz
    count+=1
    res2dim[count,1]=namevec[i]
    res2dim[count,2]=namevec[j]
    res2dim[count,3]=string(twodimcountfloat[count])
  end

  res1dim[:,1]=deepcopy(namevec)
  onedimcountfloat=float(onedimcount/sum(onedimcount))
  for j=1:sz
    res1dim[j,2]=string(Float64(onedimcountfloat[j]))
  end
  dropzero=onedimcountfloat.>0
  res1dim=res1dim[dropzero,:]
  onedimcountfloat=deepcopy(onedimcountfloat[dropzero])
  srt=sortperm(onedimcountfloat,rev=true,alg=QuickSort)

  dropzero=twodimcountfloat.>0
  twodimcountfloatdropz=deepcopy(twodimcountfloat[dropzero])
  res2dim=res2dim[dropzero,:]
  srt2=sortperm(twodimcountfloatdropz,rev=true,alg=QuickSort)
  return res1dim[srt,:],res2dim[srt2,:],onedimcount,twodimcount
end

function some_tree_settings(trnidx,validx,fixedinds::Array{Int,1},candMatWOMaxValues::Array{Array{Float64,1},1},mappings::Array{Array{String,1},1},mw::Float64,weight::Array{Float64,1},subsampling_features_prop::Float64,n::Int)
	if abs(subsampling_features_prop)<1.0
		@assert length(fixedinds)==0 "subsampling_features_prop must be 1.0 (i.e. disabled) if fixed subsets of features are provided!"
	end

	intVarsUsed=Array{Array{Int,1}}(0)	
	for i=1:length(candMatWOMaxValues)
		push!(intVarsUsed,copy(zeros(Int,length(candMatWOMaxValues[i]))))
	end	
	for i=1:length(mappings)
		push!(intVarsUsed,copy(zeros(Int,length(mappings[i]))))
	end	

	minweight=copy(mw)
	totalweight=sum(view(weight,trnidx))
	if minweight < 0.0
		minweight=max(-1.0,minweight)
		minweight= totalweight*(-minweight)
	end
    @assert !(minweight<0)
	nLeavesEstimate=Int64(round(max(1,0.5*(totalweight/minweight+0.5*totalweight/minweight))))
	nDepthEstimate=ceil(log(2,nLeavesEstimate))
	#this Depthestimate is about 3 times too low as we create very granular trees!
	nDepthToStartParallelization=Int64(round(max(1,min(nDepthEstimate-1,log(2,nprocs())))))

	if length(fixedinds) >0 #features are pre set
		inds=deepcopy(fixedinds);
	else
		inds=randomFeatureSelection(n,subsampling_features_prop)
	end
	#end
	return intVarsUsed,inds,minweight,nDepthToStartParallelization
end

function randomFeatureSelection(n_features::Int,subsampling_features_prop::Float64)
	# subsampling_features_prop>0.0 -> chose a fixed number m=ceil(subsampling_features_prop*nfeatures) of explanatory factors
	# subsampling_features_prop<0.0 -> each explanatory factor has a probability of subsampling_features_prop to be chosen.
	# subsampling_features_prop<0.0 -> here, it may be that only 1 variable is chosen, or that all variables are chosen
	@assert abs(subsampling_features_prop)>0.0
	@assert n_features>=0
	inds=Vector{Int64}(0);
	if abs(subsampling_features_prop)>=1.0
		inds=collect(1:n_features)
	else
		if subsampling_features_prop>0.0
			candidates=collect(1:n_features)
			m=Int(cld(size(candidates,1),1/subsampling_features_prop)) #m should always be >0 with this definition
			inds=sample(candidates,m,replace=false)			
		elseif (subsampling_features_prop<0.0)&&(subsampling_features_prop>-1.0)
			for i=1:n_features
			  if rand()<subsampling_features_prop;push!(inds,copy(i));end;  #copy statment is obsolete here, as we work with a value in this case
			end			
		else
			#this should not happen
			@show subsampling_features_prop
			@show subsampling_features_prop=0.0
			@show subsampling_features_prop=1.0
			@show subsampling_features_prop=-1.0
			@show abs(subsampling_features_prop)
			@show abs(subsampling_features_prop)==0.0
			@show abs(subsampling_features_prop)==1.0
			error("this should not have happend. Check the value of subsampling_features_prop")
		end
	end
#if inds is empty, pick a random feature variable (otherwise we have 0 variables which makes no sense)
if (length(inds)==0)
	push!(inds,randint(n_features)) 
end
return inds
end

function get_left_index_vector(numfeatures::Array{Float64,2},id::Int64,thresh::Float64)
  sz=size(numfeatures,1)
  res=BitArray(sz)
  for i=1:sz
      res[i]=(numfeatures[i,id]<thresh)
  end
  return res
end

function myunique(C::Array{UInt8,1},maximalNumberOfDistinctValues::Int)
#output is a vector of unique values
#this function is faster than unique in some cases as it stops when "seen" is full (it uses apriori knowledge of the maximal number of distinct values which can occur)
    out = Array{eltype(C)}(0)
    seen = Set{eltype(C)}()
	sizehint!(seen,maximalNumberOfDistinctValues)
	sizehint!(out,maximalNumberOfDistinctValues)
    count=0
    for x in C
        if !in(x, seen)
            push!(seen, x)
            push!(out, x)
            count+=1
        end
        if count>=maximalNumberOfDistinctValues;break;end;
    end
    out
end

function my_findin(small::Array{UInt8,1},big::Array{UInt8,1})
res=Array{Int}(0)
sizehint!(res,length(small))
for i=1:size(small,1)
    val=small[i]
    for j=1:length(big)
      if big[j]==val
        push!(res,i)
        break
      end
    end
end
return res
end

function complement_int_index(idx::Array{Int64,1},sz::Int64)
    set=[1:sz]
    complement=Array{Int64}(length(set)-length(idx))
    k=1
    for i=1:size(set,1)
        if in(set[i],idx)
        else
            complement[k]=i
            k+=1
        end
    end
return complement
end

function int_idx_to_bitarray(idx::Array{Int64,1},sz::Int64)
  res=BitArray(sz)
  fill!(res,false)
  for i=1:length(idx)
    res[idx[i]]=true
  end
return res
end

function BitArrayIndexMatchingSubset(a::Array{String,1},leftsubset::Array{String,1})
  #this function is incredibly inefficient but the alternative below is not yet working
  #tbd / todo: improve splitting of PDAs
  res=BitArray(length(a))
  fill!(res,false)
  match = findin(a, leftsubset)
  for i in 1:length(match)
    res[match[i]]=true
  end
  return res
end

  function BitArrayIndexMatchingSubset(a::PooledArray{String,UInt8,1},leftsubset::Array{String,1})
  #this function is incredibly inefficient but the alternative below is not yet working
  #tbd / todo: improve splitting of PDAs
  res=BitArray(length(a))
  match = findin(a.pool, leftsubset)
  for i in 1:length(a)
    res[i]=(a.refs[i] in match)
  end
  return res
end


function int_idx_to_bool_idx(leftidx::Array{Int64,1},rightidx::Array{Int64,1})
  #this can probably be done more efficiently
  maxl=maximum(leftidx)
  maxr=maximum(rightidx)
  sz=max(maxr,maxl)
  loop_over_left_idx=(maxl<maxr)
  ressmaller=Vector(Bool,sz)
  fill!(ressmaller,false)
  if loop_over_left_idx
    for i=1:maxl
      ressmaller[leftidx[i]]=true
    end
      return ressmaller,!ressmaller
   else
      for i=1:maxr
        ressmaller[rightidx[i]]=true
      end
      return !ressmaller,ressmaller
  end
end

#print Tree Fn
function print_tree(tree::Leaf,sett::ModelSettings,candMatWOMaxValues::Array{Array{Float64,1},1},mappings::Array{Array{String,1},1},indent::Int64)
    println("L|S$(round(nodesize(tree),1)),Fit$(round(tree.fitted,2)),D$(tree.depth)")
end

function print_tree(tree::Node{T},sett::ModelSettings,candMatWOMaxValues::Array{Array{Float64,1},1},mappings::Array{Array{String,1},1},indent::Int64=0) where T<:Unsigned 
	orig_id=tree.featid
	number_of_num_features=sett.number_of_num_features
	df_name_vector=sett.df_name_vector
	orig_id<0 ? this_id=number_of_num_features-orig_id : this_id=orig_id
	orig_id>0 ? println("N|F$(this_id) $(df_name_vector[this_id]),T$(signif(candMatWOMaxValues[this_id][tree.subset[end]],5)),Fit$(round(fittedratio(tree),2)),S$(round(nodesize(tree),1)))") : println("N|F$(this_id)  $(df_name_vector[this_id]),[$(join(mappings[-orig_id][[convert(Int,z) for z in tree.subset]],","))],Fit$(round(fittedratio(tree),2)),S$(round(nodesize(tree),1)))")
	#orig_id>0 ? println("N|F$(this_id),T$(round(need_to_use_candidate_matrix_here,2)),Fit$(round(fittedratio(tree),2)),S$(round(nodesize(tree),1)))") : println("N|F$(this_id),[$(join(mappings[-orig_id][[convert(Int,z) for z in tree.subset]],","))],Fit$(fittedratio(tree)),S$(round(nodesize(tree),1)))")
	print(" " ^ indent * "L=")
	print_tree(tree.left,sett,candMatWOMaxValues,mappings,indent + 1)
	print(" " ^ indent * "R=")
	print_tree(tree.right,sett,candMatWOMaxValues,mappings,indent + 1)
end

#write Tree Fn
function write_tree(candMatWOMaxValues::Array{Array{Float64,1},1},ensemble::E,number_of_num_features::Int64,indent::Int64=0,fileloc="c:\\temp\\myt.txt",df_name_vector::Array{String,1}=Array{String}(1),mappings::Array{Array{String,1},1}=Array{Array{String,1}}(0)) where {E <: Ensemble}
  iterations=size(ensemble.trees,1)
  sz=length(df_name_vector)
  onedimintvec_sum=zeros(Float64,sz)
  twodimintvec_sum=zeros(Float64,div(sz*(sz-1),2))
  p = Progress(iterations, 2, "Calculating variable importance and writing tree structure...") # minimum update interval: 2 second
  f=open(fileloc,"a+")
  for i=1:iterations
	leaves_array=create_leaves_array(ensemble.trees[i])
    write(f,"\r\n\r\nIteration number $(i)\r\n")
    #1-dim Variable Importance
    #@show size(leaves_array)
	  var_imp1d_str_arr,var_imp2d_str_arr,onedimintvec,twodimintvec=variable_importance(leaves_array,df_name_vector,number_of_num_features)
	  onedimintvec_sum+=onedimintvec
	  twodimintvec_sum+=twodimintvec
	write_tree(candMatWOMaxValues,ensemble.trees[i],number_of_num_features,var_imp1d_str_arr,var_imp2d_str_arr,0,f,df_name_vector,mappings)
	next!(p)
  end

  #close(f) #todo/tbd check if we can comment this line and the line below which opens the file again. ideally the file is only opened once and closed once.
  onedimintvec_avg=onedimintvec_sum/sum(onedimintvec_sum)
  twodimintvec_avg=twodimintvec_sum/sum(twodimintvec_sum)

  res1dim=Array{String}(sz,2)
  res2dim=Array{String}(div(sz*(sz-1),2),3)
  count=0
  for i=1:sz, j=i+1:sz
    count+=1
    res2dim[count,1]=df_name_vector[i]
    res2dim[count,2]=df_name_vector[j]
    res2dim[count,3]=string(twodimintvec_avg[count])
  end

  res1dim[:,1]=deepcopy(df_name_vector)
  for j=1:sz
    res1dim[j,2]=string(Float64(onedimintvec_avg[j]))
  end
  dropzero=onedimintvec_avg.>0
  res1dim=res1dim[dropzero,:]
  onedimintvec_avg=deepcopy(onedimintvec_avg[dropzero])
  srt=sortperm(onedimintvec_avg,rev=true,alg=QuickSort)

  dropzero=twodimintvec_avg.>0
  twodimintvec_avg=deepcopy(twodimintvec_avg[dropzero])
  res2dim=res2dim[dropzero,:]
  srt2=sortperm(twodimintvec_avg,rev=true,alg=QuickSort)

  #f=open(fileloc,"a+")
    write(f,"\r\n\r\nOverall predictor importance of boosted tree: \r\n")
    my_write(f,res1dim[srt,:])
    my_write(f,res2dim[srt2,:])
  #bagging only
  if typeof(ensemble)==BaggedTree
	#wei=hcat([1:iterations],weightPerTree)
	write(f,"\r\n\r\nWeights per Iteration")
	for i=1:iterations
		write(f,i," ",ensemble.weightPerTree[i])
	end
  end
  close(f)

  return nothing
end

function write_tree(candMatWOMaxValues::Array{Array{Float64,1},1},tree::Tree,number_of_num_features::Int64,var_imp1d_str_arr::Array{String,2},var_imp2d_str_arr::Array{String,2},indent::Int64,fileloc::String,df_name_vector::Array{String,1}=Array{String}(1),mappings::Array{Array{String,1},1}=Array{Array{String,1}}(0))
	#version which opens fileloc
	f=open(fileloc,"a+")
		write_tree(candMatWOMaxValues,tree.rootnode,number_of_num_features,var_imp1d_str_arr,var_imp2d_str_arr,indent,f,df_name_vector,mappings)
    close(f)
end

function write_tree(candMatWOMaxValues::Array{Array{Float64,1},1},tree::Leaf,number_of_num_features::Int64,var_imp1d_str_arr::Array{String,2},var_imp2d_str_arr::Array{String,2},indent::Int64,f::IOStream,df_name_vector::Array{String,1}=Array{String}(1),mappings::Array{Array{String,1},1}=Array{Array{String,1}}(0))
#function in the special case where the whole tree is in fact a single leaf
	#version which uses f::IOStream as input
	write(f,"\r\n1 dimensional predictor importance:")
	my_write(f,var_imp1d_str_arr);
	write(f,"\r\n2 dimensional predictor importance:")
	my_write(f,var_imp2d_str_arr);
	write(f,"Number of observations in this tree: $(nodesize(tree))");write(f,"\r\n");
    write(f,"Number of leaves: $(size(create_leaves_array(tree),1))");write(f,"\r\n");
    write(f,"Number of Nodes: $(number_of_nodes(tree))");write(f,"\r\n");
    write(f,"Maximal depth of the tree: $(maxdepth(tree))");write(f,"\r\n");
    write(f,"\r\n")
    write(f,"This tree is equal to a single leaf.\r\n")

    write(f,"\r\n")
    write(f," " ^ indent * "Leaf=")
    write_tree(candMatWOMaxValues,tree,number_of_num_features, indent + 1,f,df_name_vector,mappings)
end

function write_tree(candMatWOMaxValues::Array{Array{Float64,1},1},tree::Node{T},number_of_num_features::Int64,var_imp1d_str_arr::Array{String,2},var_imp2d_str_arr::Array{String,2},indent::Int64,f::IOStream,df_name_vector::Array{String,1}=Array{String}(1),mappings::Array{Array{String,1},1}=Array{Array{String,1}}(0)) where T<:Unsigned 
	#version which uses f::IOStream as input
	write(f,"\r\n1 dimensional predictor importance:")
	my_write(f,var_imp1d_str_arr);
	write(f,"\r\n2 dimensional predictor importance:")
	my_write(f,var_imp2d_str_arr);
	orig_id=tree.featid
	orig_id<0 ? this_id=number_of_num_features-orig_id : this_id=orig_id
    write(f,"Number of observations in this tree: $(nodesize(tree))");write(f,"\r\n");
    write(f,"Number of leaves: $(size(create_leaves_array(tree),1))");write(f,"\r\n");
    write(f,"Number of Nodes: $(number_of_nodes(tree))");write(f,"\r\n");
    write(f,"Maximal depth of the tree: $(maxdepth(tree))");write(f,"\r\n");
    write(f,"\r\n")
    write(f,"\r\n")

    orig_id>0 ? write(f,"N|F$(this_id) $(df_name_vector[this_id]),T$(signif(candMatWOMaxValues[this_id][tree.subset[end]],5)),Fit$(round(fittedratio(tree),2)),S$(round(nodesize(tree),1)))") : write(f,"N|F$(this_id)  $(df_name_vector[this_id]),[$(join(mappings[-orig_id][[convert(Int,z) for z in tree.subset]],","))],Fit$(round(fittedratio(tree),2)),S$(round(nodesize(tree),1)))")
    write(f,"\r\n")
    write(f," " ^ indent * "L=")
    write_tree(candMatWOMaxValues,tree.left,number_of_num_features, indent + 1,f,df_name_vector,mappings)
    write(f," " ^ indent * "R=")
    write_tree(candMatWOMaxValues,tree.right,number_of_num_features, indent + 1,f,df_name_vector,mappings)
end

function write_tree(candMatWOMaxValues::Array{Array{Float64,1},1},tree::Leaf,number_of_num_features::Int64, indent::Int64=0,fileloc::String="c:\\temp\\myt.txt",df_name_vector::Array{String,1}=Array{String}(1),mappings::Array{Array{String,1},1}=Array{Array{String,1}}(0))
    warn("This should not happen. Opening file $(fileloc) to write Leaf information.")
	#This should only happen when the whole tree is a leaf.
	#We should by any means avoid to open and close the file for each leaf.
	f=open(fileloc,"a+")
	write_tree(candMatWOMaxValues,tree,number_of_num_features, indent ,f,df_name_vector,mappings)
    close(f)
end

function write_tree(candMatWOMaxValues::Array{Array{Float64,1},1},tree::Leaf,number_of_num_features::Int64,indent::Int64,f::IOStream,df_name_vector::Array{String,1}=Array{String}(1),mappings::Array{Array{String,1},1}=Array{Array{String,1}}(0))
    write(f,"L|S$(round(nodesize(tree),1)),Fit$(round(tree.fitted,2)),D$(tree.depth)\r\n")
end

function write_tree(candMatWOMaxValues::Array{Array{Float64,1},1},tree::Node{T},number_of_num_features::Int64, indent::Int64,f::IOStream,df_name_vector::Array{String,1}=Array{String}(1),mappings::Array{Array{String,1},1}=Array{Array{String,1}}(0)) where T<:Unsigned 
	orig_id=tree.featid
	orig_id<0 ? this_id=number_of_num_features-orig_id : this_id=orig_id
    orig_id>0 ? write(f,"N|F$(this_id)  $(df_name_vector[this_id]),T$(signif(candMatWOMaxValues[this_id][tree.subset[end]],5)),Fit$(round(fittedratio(tree),2)),S$(round(nodesize(tree),1)))") : write(f,"N|F$(this_id)  $(df_name_vector[this_id]),[$(join(mappings[-orig_id][[convert(Int,z) for z in tree.subset]],","))],Fit$(round(fittedratio(tree),2)),S$(round(nodesize(tree),1)))")
    write(f,"\r\n")
    write(f," " ^ indent * "L=")
    write_tree(candMatWOMaxValues,tree.left, number_of_num_features,indent + 1,f,df_name_vector,mappings)
    write(f," " ^ indent * "R=")
    write_tree(candMatWOMaxValues,tree.right,number_of_num_features, indent + 1,f,df_name_vector,mappings)
end

#write SAS Code fn
function write_sas_code(estimatesPerScore,candMatWOMaxValues::Array{Array{Float64,1},1},bt::BoostedTree,number_of_num_features::Int64,fileloc::String,df_name_vector::Array{String,1},settings::String,mappings::Array{Array{String,1},1}=Array{Array{String,1}}(0);leafvarname::String=convert(String,"leaf"),indent::Int64=0)	
	
iterations=size(bt.trees,1)
local vname
notsign=convert(String,"not")
insign=convert(String," in ")
quotestr=convert(String,"'")
listbeg=convert(String,"(")
listend=convert(String,")")
utf8ListVarsUsed,boolListVarsUsed=variablesUsed(bt)

fiostream=open(fileloc,"w")
#write BOM
write(fiostream,global_byte_order_mark)
write(fiostream,"/* \r\nSettings: \r\n",settings,"\r\n*/ \r\n")
write(fiostream,"/* \r\n Variables used by model: \r\n ",join(utf8ListVarsUsed,' ')," \r\nTODO:\r\n Replace the value for missing character values (default: '#') EVERYWHERE IN THE CODE with \'\' \r\n Remove the value for \"other\" character values (e.g.~) from the if conditions at the top (do not replace them in the actual model code further down). \r\n Define the lists of known character variables for each character variable. \r\n You can use the SAS macro %define_known_char_varlist(data=runmodel); for this. \r\n Define the macro variable RARE_CHARACTER_VALUE \r\n %let RARE_CHARACTER_VALUE=\'~\'; \r\n*/\r\n")
txtw="""
%put WARNING: (bk) All character variables need to be in uppercase in the data otherwise SAS will produce wrong results!;
"""
write(fiostream,"\r\n",txtw,"\r\n")

write(fiostream,"option nosource; \r\n\r\ndata vboosting/view=vboosting;\r\r\n format fitted_value fitted_value_best_guess best32. variables_have_unknown_values 8. variables_with_unknown_values \$1000. variables_have_rare_values 8. variables_with_rare_values \$1000.; format ")
for i=1:iterations
      write(fiostream,leafvarname,"_iter$(i) ")
end
write(fiostream,"; \r\nset runmodel;\r\nvariables_have_unknown_values=0;\r\nvariables_have_rare_values=0;\r\nvariables_with_unknown_values=\'\';\r\nvariables_with_rare_values=\'\';")
indexOfLastCharVarUsed=0
#data prep and check
	for i in 1:length(mappings)
		#upcase für char vars
		vname=df_name_vector[number_of_num_features+i]
		if boolListVarsUsed[number_of_num_features+i]
			indexOfLastCharVarUsed=max(i,indexOfLastCharVarUsed)
			write(fiostream,vname,"=upcase(",vname,");\r\n")
		end
	end
	for i in 1:length(mappings)
		vid=number_of_num_features+i
		vname=df_name_vector[vid]
		if boolListVarsUsed[number_of_num_features+i]
			write(fiostream,"/*Variable $(vname)*/\r\n")
			write(fiostream,"if ",vname," not in (&",uppercase(vname),"_VALS.) then do; _u_",vname,"=1;variables_have_unknown_values=1;end;else _u_",vname,"=0;\r\n")
			write(fiostream,"if ",vname," not in (",quotestr,join(mappings[i],string(quotestr," ",quotestr)),quotestr,") then ",vname,"=&RARE_CHARACTER_VALUE.;\r\n")
		end
	end

if length(mappings)>0&&(sum(boolListVarsUsed[number_of_num_features+1:end])>0)
	write(fiostream,"\r\nvariables_with_unknown_values=catx(\',\',\r\n")
	for i in 1:length(mappings)
		vid=number_of_num_features+i
		vname=df_name_vector[vid]
		if boolListVarsUsed[number_of_num_features+i]
			write(fiostream,"\tifc(_u_",vname,",\'",vname,"\',\'\')")
			if i<indexOfLastCharVarUsed
				write(fiostream,",\r\n")
			end
		end
	end
	write(fiostream,"\r\n);\r\n")
#variables_with_rare_values
	write(fiostream,"\r\nvariables_with_rare_values=catx(\',\',\r\n")
	for i in 1:length(mappings)
		vid=number_of_num_features+i
		vname=df_name_vector[vid]
		if boolListVarsUsed[number_of_num_features+i]
			write(fiostream,"\tifc(",vname,"=&RARE_CHARACTER_VALUE.,\'",vname,"\',\'\')")
			if i<indexOfLastCharVarUsed
				write(fiostream,",\r\n")
			end
		end
	end
	write(fiostream,"\r\n);\r\n")
end

#write tree code
    for i=1:iterations
		if typeof(bt.trees[i])!=Leaf
			orig_id=bt.trees[i].featid
			orig_id<0 ? this_id=number_of_num_features-orig_id : this_id=orig_id
			write(fiostream,"\r\n\r\n*Iteration $(i);\r\n")
			if orig_id>0
				write(fiostream,"if ",df_name_vector[this_id],"<=$(candMatWOMaxValues[this_id][bt.trees[i].subset[end]]) then do;\r\n")
					write_tree_at_each_node!(candMatWOMaxValues,bt.trees[i].left,number_of_num_features, indent+1,fiostream,df_name_vector,mappings,string(leafvarname,"_iter$(i)"),bt.moderationvector[i])
				write(fiostream," end;else do;\r\n")
					write_tree_at_each_node!(candMatWOMaxValues,bt.trees[i].right,number_of_num_features, indent+1,fiostream,df_name_vector,mappings,string(leafvarname,"_iter$(i)"),bt.moderationvector[i])
				write(fiostream," end;\r\n")
			else
				write(fiostream,"if ",df_name_vector[this_id],insign,listbeg,quotestr,join(mappings[-orig_id][[convert(Int,z) for z in bt.trees[i].subset]],string(quotestr," ",quotestr)),quotestr,listend," then do;\r\n")
					write_tree_at_each_node!(candMatWOMaxValues,bt.trees[i].left,number_of_num_features, indent+1,fiostream,df_name_vector,mappings,string(leafvarname,"_iter$(i)"),bt.moderationvector[i])
				write(fiostream," end;else do;\r\n")
					write_tree_at_each_node!(candMatWOMaxValues,bt.trees[i].right,number_of_num_features, indent+1,fiostream,df_name_vector,mappings,string(leafvarname,"_iter$(i)"),bt.moderationvector[i])
				write(fiostream," end;\r\n")
			end
		else
			#"Tree is not a node (hence it must be a Leaf). Code to write SAS code for a single leaf is not yet implemented"
			write_tree_at_each_node!(candMatWOMaxValues,bt.trees[i],number_of_num_features,indent,fiostream,df_name_vector,mappings,string(leafvarname,"_iter$(i)"),bt.moderationvector[i])
		end
	end

	#todo tbd, take note of this somewhere and do not forget it!
	#for the code which is generated here the last entry of maxRawRelativityPerScoreSorted is modified to be 50% (50% is arbitrary here) higher than the value which was derived by Julia
	#this means that any policy will get the maximal score if it is "worse" than any training policy
	#but if it is much worse, it will get the default score of -1
	modified_maxRawRelativityPerScoreSorted=deepcopy(bt.maxRawRelativityPerScoreSorted)
	modified_maxRawRelativityPerScoreSorted[end]=modified_maxRawRelativityPerScoreSorted[end]*1.5
	sas_write_ScoreMap(fiostream,modified_maxRawRelativityPerScoreSorted,estimatesPerScore)

	# old approach write(fiostream,"\r\nfitted_value=$(bt.meanobserved); \r\narray tmpnum{*} rel_mod_leaf_iter:; \r\ndo i_loop=1 to dim(tmpnum); \r\nfitted_value=fitted_value*tmpnum[i_loop]; \r\nend;\r\nif variables_have_unknown_values then do;\r\nfitted_value_best_guess=fitted_value;fitted_value=.;\r\nend;\r\nelse do;fitted_value_best_guess=fitted_value;\r\nend;\r\nrun; \r\n\r\noption source; \r\n")
	
	endtxt="""

	data fitted2;
		format irkey;
		format score raw_relativity;
		set vboosting(drop=i_loop leaf_iter: rel_mod_leaf_iter:);			
	run;

	proc sort data=res out=res_sorted;by irkey;quit;
	proc sort data=fitted out=fitted_sorted;by irkey;quit;
	proc compare criterion=1E-14 out=compared data=res_sorted compare=fitted_sorted(rename=fitted_value=est_from_trn);id irkey;var est_from_trn;quit;
	proc compare criterion=1E-14 out=compared data=res_sorted compare=fitted_sorted(rename=fitted_value=rawEstimatePerObs);id irkey;var rawEstimatePerObs;quit;

/*

	proc compare criterion=1E-14 out=compared data=res_sorted(rename=(est_from_trn=est_from_trn_other rawEstimatePerObs=est_from_trn)) compare=fitted_sorted(rename=fitted_value=est_from_trn);id irkey;var est_from_trn;quit;
	data sas_leafnrs;set fitted_sorted(keep=irkey leaf_iter:);run;
	proc sort data=julia_leafnrs;by irkey;quit;
	proc compare criterion=1E-14 out=cmp2 data=sas_leafnrs compare=julia_leafnrs;id irkey;quit;
	data test;set fitted;keep fitted_value fitted_value_best_guess &var_dep.;run;

*/

	"""
	write(fiostream,endtxt)
close(fiostream)
end

"""write_sas_code version for single trees"""
function write_sas_code(candMatWOMaxValues::Array{Array{Float64,1},1},inTree::Tree,number_of_num_features::Int64,fileloc::String,df_name_vector::Array{String,1},settings::String,mappings::Array{Array{String,1},1}=Array{Array{String,1}}(0),mdf::Float64=1.0;leafvarname::String=convert(String,"leaf"),indent::Int64=0)
tree=inTree.rootnode
notsign=convert(String,"not")
insign=convert(String," in ")
quotestr=convert(String,"'")
listbeg=convert(String,"(")
listend=convert(String,")")
fiostream=open(fileloc,"w")

  #write beginning
  #write BOM
	write(fiostream,global_byte_order_mark)
	write(fiostream,"/* \r\nSettings: \r\n WARNING: common issues with the SAS code: the data in sas may have leading zeros (categorical columns) which can be lost when exporting and importing!\r\n ")
	write(fiostream,"\r\n WARNING: common issues with the SAS code: variables are CUT off (at the end) if modified within SAS.\r\n ")
    write(fiostream,settings)
	write(fiostream,"\r\n*/ \r\n\r\n")	
    write(fiostream,"data sas_tree;format ",leafvarname,";set runmodel;;")
    write(fiostream,"\r\n")

if (typeof(tree)!=Leaf)  #in an earlier version we had ==Node ; however now node has a paramter -> (x==Node{UInt8}||x==Node{UInt16})
	orig_id=tree.featid
	orig_id<0 ? this_id=number_of_num_features-orig_id : this_id=orig_id

	if orig_id>0
		write(fiostream,"if ",df_name_vector[this_id],"<=$(candMatWOMaxValues[this_id][tree.subset[end]]) then do;\r\n\r\n")
			write_tree_at_each_node!(candMatWOMaxValues,tree.left,number_of_num_features, indent+1,fiostream,df_name_vector,mappings,leafvarname,mdf)
		write(fiostream," end;else do;\r\n")
			write_tree_at_each_node!(candMatWOMaxValues,tree.right,number_of_num_features, indent+1,fiostream,df_name_vector,mappings,leafvarname,mdf)
		write(fiostream," end;\r\n")
	else
		write(fiostream,"if ",df_name_vector[this_id],insign,listbeg,quotestr,join(mappings[-orig_id][[convert(Int,z) for z in tree.subset]],string(quotestr," ",quotestr)),quotestr,listend," then do;\r\n")
			write_tree_at_each_node!(candMatWOMaxValues,tree.left,number_of_num_features, indent+1,fiostream,df_name_vector,mappings,leafvarname,mdf)
		write(fiostream," end;else do;\r\n")
			write_tree_at_each_node!(candMatWOMaxValues,tree.right,number_of_num_features, indent+1,fiostream,df_name_vector,mappings,leafvarname,mdf)
		write(fiostream," end;\r\n")
	end
else
	#"Tree is not a node (hence it must be a Leaf). Code to write SAS code for a single leaf is not yet implemented"
	#@assert typeof(tree)==Leaf
	@info "DTM: Tree was a single leaf. No proper SAS Code was produced."
	@assert false
	write(fiostream," /*ERROR; the tree was a single leaf*/")
	#write_tree_at_each_node!(candMatWOMaxValues,tree,number_of_num_features,indent,fiostream,df_name_vector,mappings,leafvarname,mdf)
end
  #write rest of SAS code
  write(fiostream,"run; \r\nproc summary data=sas_tree missing;class ",leafvarname,";var &var_dep.;types ",leafvarname,";output out=_summary_( rename=(_freq_=n)) sum=;");
  write(fiostream,"where training_validation_bool eq 1;run; \r\ndata summary;set _summary_;est_sas=&var_dep./n;run; \r\nproc sql noprint;create table est_sas as select s.est_sas as est_from_trn,i.* from sas_tree as i ")
  write(fiostream,"left join summary as s on (s.",leafvarname,"=i.",leafvarname,") order by irkey;quit; \r\n")
  write(fiostream,"proc sort data=runmodel;by irkey;quit; \r\nproc sort data=res;by irkey;quit; \r\nproc compare out=compare base=res compare=est_sas;var &var_dep.;id irkey;quit; \r\n")
 #close file
close(fiostream)
end

function write_tree_at_each_node!(candMatWOMaxValues::Array{Array{Float64,1},1},tree::Node{T},number_of_num_features::Int64, indent::Int64,fiostream::IOStream,df_name_vector::Array{String,1}=Array{String}(1),mappings::Array{Array{String,1},1}=Array{Array{String,1}}(0),leafvarname::String=convert(String,"leaf"),mdf::Float64=1.0) where T<:Unsigned 
	orig_id=tree.featid
	orig_id<0 ? this_id=number_of_num_features-orig_id : this_id=orig_id
	notsign=convert(String,"not")
	insign=convert(String," in ")
	quotestr=convert(String,"'")
	listbeg=convert(String,"(")
	listend=convert(String,")")
	write(fiostream," " ^ indent)
	if orig_id>0
		write(fiostream,"if ",df_name_vector[this_id],"<=$(candMatWOMaxValues[this_id][tree.subset[end]]) then do;\r\n")
			write_tree_at_each_node!(candMatWOMaxValues,tree.left,number_of_num_features, indent+1,fiostream,df_name_vector,mappings,leafvarname,mdf)
		write(fiostream," " ^ (indent-1))
		write(fiostream," end;else do;\r\n")
			write_tree_at_each_node!(candMatWOMaxValues,tree.right,number_of_num_features, indent+1,fiostream,df_name_vector,mappings,leafvarname,mdf)
		write(fiostream," " ^ (indent-1))
		write(fiostream," end;\r\n")
	else
		write(fiostream,"if ",df_name_vector[this_id],insign,listbeg,quotestr,join(mappings[-orig_id][[convert(Int,z) for z in tree.subset]],string(quotestr," ",quotestr)),quotestr,listend," then do;\r\n")
			write_tree_at_each_node!(candMatWOMaxValues,tree.left,number_of_num_features, indent+1,fiostream,df_name_vector,mappings,leafvarname,mdf)
		write(fiostream," " ^ (indent-1))
		write(fiostream," end;else do;\r\n")
			write_tree_at_each_node!(candMatWOMaxValues,tree.right,number_of_num_features, indent+1,fiostream,df_name_vector,mappings,leafvarname,mdf)
		write(fiostream," " ^ (indent-1))
		write(fiostream," end;\r\n")
	end
end

function write_tree_at_each_node!(candMatWOMaxValues::Array{Array{Float64,1},1},tree::Leaf,number_of_num_features::Int64,indent::Int64,fiostream::IOStream,df_name_vector::Array{String,1}=Array{String}(1),mappings::Array{Array{String,1},1}=Array{Array{String,1}}(0),leafvarname::String=convert(String,"leaf"),mdf::Float64=1.0)
	writecomments=false
  if writecomments
    add=" *Nodesize: $(round(nodesize(tree),1)), Fitted Value: $(round(tree.fitted,2)), Depth: $(tree.depth);\r\n"
  else
    add="\r\n"
  end
	write(fiostream," " ^ indent,leafvarname,"=$(tree.id); $(add)")
    #this row is only meaningful for boosted multiplicative trees
	#warn("check if this definition of val is accurate for a boosted tree (e.g. lr model)"
	val=_moderate(tree.fitted,mdf)
	#val=tree.fitted
	write(fiostream," " ^ indent,"rel_mod_",leafvarname,"=$(val);\r\n")
end

function write_sas_code(leaves::Array{Leaf,1},number_of_num_features::Int64,fileloc::String,namevec::Array{String,1},settings::String,mappings::Array{Array{String,1},1}=Array{Array{String,1}}(0);leafvarname::String=convert(String,"leaf"))
	error("BK: This function should not be used anymore. (right?)")	
	#open file
fiostream=open(fileloc,"w")
  #write beginning
  #write BOM
	write(fiostream,global_byte_order_mark)
  write(fiostream,"/* \r\nSettings: \r\n")
  write(fiostream,settings)
  write(fiostream,"\r\n*/ \r\n")
    write(fiostream,"data sas_tree;format ",leafvarname,";set runmodel;")
    write(fiostream,"\r\n")
  #write leaf definitions
    write_rules_to_file(leaves,fiostream,namevec,number_of_num_features,mappings)
  #write rest of SAS code
   write(fiostream,"run; \r\nproc summary data=sas_tree missing;class ",leafvarname,";var &var_dep.;types ",leafvarname,";output out=_summary_( rename=(_freq_=n)) sum=;");
  write(fiostream,"where training_validation_bool eq 1;run; \r\ndata summary;set _summary_;est_sas=&var_dep./n;run; \r\nproc sql noprint;create table est_sas as select s.est_sas as est_from_trn,i.* from sas_tree as i ")
  write(fiostream,"left join summary as s on (s.",leafvarname,"=i.",leafvarname,") order by irkey;quit; \r\n")
  write(fiostream,"proc sort data=runmodel;by irkey;quit; \r\nproc compare out=compare base=res compare=est_sas;var &var_dep.;id irkey;quit; \r\n")
 #close file
close(fiostream)
end

function write_rules_to_file(leaves::Array{Leaf,1},f::IOStream,namevec::Array{String,1},number_of_num_features,mappings::Array{Array{String,1},1}=Array{Array{String,1}}(0);leafvarname::String=convert(String,"leaf"))
  for i=1:size(leaves,1)
    if i==1
      write(f,"if ")
    else
      write(f,"else if ")
    end
    write(f,rulpath_vector_to_str(leaves[i].rule_path,namevec,number_of_num_features,mappings))
    write(f," then $(leafvarname)=$(i);")
    write(f,"\r\n")
  end
end

function rulpath_vector_to_str(rp::Array{Rulepath,1},namevec::Array{String,1},number_of_num_features::Int64,mappings::Array{Array{String,1},1}=Array{Array{String,1}}(0))
  result=convert(String,"")
  for i=1:size(rp,1)
    if i>1
      result=convert(String,string(result," and ",reconstruct_rule(rp[i],namevec,number_of_num_features,mappings)))
    else
      result=convert(String,reconstruct_rule(rp[i],namevec,number_of_num_features,mappings))
    end
end
  result=convert(String,string("(",result,")"))
return result
end

function reconstruct_rule(r::Rulepath,namevec::Array{String,1},number_of_num_features,mappings::Array{Array{String,1},1}=Array{Array{String,1}}(0);notsign=convert(String,"not")::String,insign=convert(String," in ")::String,quotestr=convert(String,"'")::String,listbeg=convert(String,"(")::String,listend=convert(String,")")::String)
  orig_id=r.featid
  orig_id<0 ? this_id=number_of_num_features-orig_id : this_id=orig_id
  #[$(join(mappings[-orig_id][[convert(Int,z) for z in bt.trees[i].subset]],","))]
  if orig_id<0
    result=convert(String,string("(",namevec[this_id],insign,listbeg,quotestr, join(mappings[-orig_id][[convert(Int,z) for z in r.subset]],string(quotestr," ",quotestr)),quotestr,listend,")"))
  else
    error("this is not working anymore, to be modified")
    result=convert(String,string("(",namevec[this_id],"<",r.featthresholdvalue,")"))
  end
  if !(r.isLeftChild)
    result=convert(String,string(notsign,"(",result,")"))
  end
  #note "(" is not the same as '('
  if !(result[1]=='(')
    result=convert(String,string("(",result,")"))
  end
return result
end

#Statistics
mean(tree::Tree)=mean(tree.rootnode)
function mean(tree::Node{T}) where T<:Unsigned  
return ((nodesize(tree.left)*mean(tree.left)+nodesize(tree.right)*mean(tree.right))/nodesize(tree))
end
mean(tree::Leaf)=tree.mean

fittedratio(tree::Tree)=fittedratio(tree.rootnode)
function fittedratio(tree::Node{T}) where T<:Unsigned 
return numeratorsum(tree)/denominatorsum(tree)
end #((nodesize(tree.left)*fittedratio(tree.left)+nodesize(tree.right)*fittedratio(tree.right))/nodesize(tree))
fittedratio(tree::Leaf)=numeratorsum(tree)/denominatorsum(tree)

numeratorsum(tree::Tree)=numeratorsum(tree.rootnode)
function numeratorsum(tree::Node{T}) where T<:Unsigned 
return numeratorsum(tree.left)+numeratorsum(tree.right)
end
numeratorsum(x::Leaf)=x.sumnumerator

denominatorsum(tree::Tree)=denominatorsum(tree.rootnode)
function denominatorsum(tree::Node{T}) where T<:Unsigned  
return denominatorsum(tree.left)+denominatorsum(tree.right)
end
denominatorsum(x::Leaf)=x.sumdenominator

nodesize(tree::Tree)=nodesize(tree.rootnode)
	function nodesize(tree::Node{T}) where T<:Unsigned 
return nodesize(tree.left)+nodesize(tree.right)
end
nodesize(tree::Leaf)=tree.size

rowcount(tree::Tree)=rowcount(tree.rootnode)
function rowcount(tree::Node{T}) where T<:Unsigned  
return rowcount(tree.left)+rowcount(tree.right)
end
rowcount(tree::Leaf)=tree.rowcount::Int

maxdepth(tree::Tree)=maxdepth(tree.rootnode)
function maxdepth(tree::Node{T}) where T<:Unsigned  
return (1+max(maxdepth(tree.left),maxdepth(tree.right)))
end
maxdepth(tree::Leaf)=0

#helperfunctions

function ismonotonic(x)
    if endof(x)<2
        return true
    else
    for i=1:endof(x)-1
        if x[i]>x[i+1]
            return false
        end
    end
    end
return true
end

function _moderate!(idx_to_be_moderated::Vector{Int},estimatedRatio::Array{Float64,1},indicatedRelativity::Array{Float64,1},current_mdf::Float64)
	@assert size(estimatedRatio)==size(indicatedRelativity)
	for i in idx_to_be_moderated
		@inbounds estimatedRatio[i]=((indicatedRelativity[i]-1.0)*current_mdf+1.0)*estimatedRatio[i]
	end
	return nothing
end

function _moderate!(estimatedRatio::Array{Float64,1},indicatedRelativity::Array{Float64,1},current_mdf::Float64)
	@assert size(estimatedRatio)==size(indicatedRelativity)
	for i=1:length(estimatedRatio)
		@inbounds estimatedRatio[i]=((indicatedRelativity[i]-1.0)*current_mdf+1.0)*estimatedRatio[i]
	end
	return nothing
end

function _moderate(rel::Array{Float64,1},moderationfactor::Float64)
    l=size(rel,1)
    res=Array{Float64}(l)
    for i=1:l
		res[i]=_moderate(rel[i],moderationfactor)
    end
    return res
end

function _moderate(rel::Float64,moderationfactor::Float64)
    if rel<1.0 #this would probably also work without if (rel-1.0)*moderationfactor+1.0
	    return (1.0-(1.0-rel)*moderationfactor)
	  else
		  return ((rel-1.0)*moderationfactor+1.0)
	end
end

function char_rule(subset::Array{String,1},symbol_column::Symbol,b::Bool)
   if size(subset,1)==1
        #b ? (return (:($(parse(string(symbol_column))).==$(subset[1])))) : (return (:($(parse(string(symbol_column))).!=$(subset[1]))))
        b ? (return (:($(parse("dataframename[:$(symbol_column)]")).==$(subset[1])))) : (return (:($(parse("dataframename[:$(symbol_column)]")).!=$(subset[1]))))
    else
        r=:($(parse("dataframename[:$(symbol_column)]")).==$(subset[1]))
        for i=2:size(subset,1)
            rtmp=:($(parse("dataframename[:$(symbol_column)]")).==$(subset[i]))
            r=:($(rtmp) | $(r))
        end
        b ? (return r) : (:(!($(r))))
    end
end

function inverse_permutation(a::Array{T,1}) where {T}
  #returns the "inverse permutation" for a permutation vector such as srt=sortperm(rand(10))
  res=similar(a)
  for i=1:length(a)
    res[a[i]]=i
  end
  return res
end

function sort_splitlist(sd::Vector{Splitdef{T}}) where T<:Unsigned #::Array{Splitdef,1})
    ar=Array{Float64}(0)
	sizehint!(ar,size(sd,1))
    #todo, sort by splitvalue, feature_id, threshold and some interpretation of subset
    #otherwise there can be multiple splits with the same purity
    #I noted that the different settings with parallelization do not yield the exact same results
    for i=1:size(sd,1)
		#get rid of splits with splitvalue=-Inf
		if isfinite(sd[i].splitvalue)
			push!(ar,sd[i].splitvalue)
		end
    end
    p=sortperm(ar,rev=true,alg=QuickSort)
    return sd[p]
end

function subset_splitlist(sd::Vector{Splitdef{T}},minweight::Float64) where T<:Unsigned
	#outputs only these members that have enough weight
	res=Vector{Splitdef{T}}(0)
	for i=1:size(sd,1)
		if min(sd[i].weightr,sd[i].weightl)>=minweight
			push!(res,sd[i])
		end
	end
	return res
end

function index_vector_from_rulepath(rpArray::Array{Rulepath,1},numfeatures::Array{Float64,2},charfeatures::Array{String,2})
  nobs=max(size(numfeatures,1),size(charfeatures,1))
  res=Array{Bool}(nobs)
  fill!(res,true)
    for i=1:size(rpArray,1)
      rule=rpArray[i]
      if rule.featid>0
        error("this is not working anymore, to be modified")
		restmp=(numfeatures[:,rule.featid].<rule.featthresholdvalue)
        if !(rule.isLeftChild); restmp=.!restmp;end;
      else
        restmp=Array{Bool}(size(charfeatures,1))
        fill!(restmp,true)
        for j=1:size(rule.subset,1)
           restmp=restmp | (charfeatures[:,-rule.featid] .== rule.subset[j])
        end
      end
      #Apply "subsetting"
      res=(res & (restmp))
  end
  return res
end

function df_column_has_na(df::DataFrame,colindex::Int64)
  size(findin(isna(df[colindex]), true),1)
end

function defineleftrightPDA(left::BitVector,right::BitVector,mat_charfeatures::Array{String,2})
  number_of_char_features=size(mat_charfeatures,2)
  res_charfeatures_l=Array{PooledArray{String,UInt8}}(number_of_char_features)
  res_charfeatures_r=Array{PooledArray{String,UInt8}}(number_of_char_features)
  for i=1:number_of_char_features
    res_charfeatures_l[i]=PooledArray(mat_charfeatures[left,i],UInt8)
    res_charfeatures_r[i]=PooledArray(mat_charfeatures[right,i],UInt8)
  end
  return res_charfeatures_l,res_charfeatures_r
end

  function check_for_missing_data(dfIndata::DataFrame,const_shift_cols::Int64)
	#Ensure that data has no missing values
	namevec=names(dfIndata)
	for ii=1:size(dfIndata,2)
	  if ii!=const_shift_cols
      if any(ismissing,dfIndata[ii])
		idx=ismissing(dfIndata[ii])
		narows=find(idx)
		@assert length(narows)>0 "this should not have happend" #narows should be >0 by definition
		warn("Missing values detected. Rows= $(narows[1]), Column=$(ii):$(namevec[ii])")
        error("ABORT.")
        return false
      end
	  end
	end
	return true
end

function check_for_missing_data(mat::Array{Any,2},const_shift_cols::Int64,header::Array{AbstractString,1})
	#Ensure that data has no missing values
	for ii=1:size(mat,2)
	  if ii!=const_shift_cols
      col=view(mat,:,ii)
	    miss=col.==""
	    if sum(miss)>0
		  missingentries=findin(miss,[true])
		  missingentries=copy(missingentries[1:min(length(missingentries),5)])
        error("Column $(ii):$(header[ii]) contains missing values. Consider for instance these row numbers: $(join(missingentries,","))! Missing data is currently not supported. ABORT.")
        return false
      end
	  end
	end
	return true
end

function return_file_and_folders(ARGS)
  fallback_datafolder=convert(String,Main.dictGlobalSettings["Temp Folder"]);
  fallback_inputSettingsFile=convert(String,string(fallback_datafolder,"\\julia.settings.csv"))
  fallback_inputDataFile=convert(String,string(fallback_datafolder,"\\julia.csv"))
  fallback_outfileStringOnly=convert(String,"out_julia")
  fallback_outputFilename=convert(String,string(fallback_datafolder,fallback_outfileStringOnly))
	  #const this_hostname=gethostname();
	  fallback_datafolder=convert(String,Main.dictGlobalSettings["Temp Folder"])
	  #if (this_hostname=="ZURI20A7-CTO1WW");    fallback_datafolder=convert(String,"r:\\temp");  end
	  #if ((this_hostname=="MUNICH1")||(this_hostname=="MUNICH2")||(this_hostname=="MUNICH3"));    fallback_datafolder=convert(String,"C:\\temp");  end;
	function path_readable(s);    res=false;    try (res=isreadable(s));    catch (res=false);    end;    res;  end;
  #convert ARGS to String
  ASCIIStringARGS=Array{String}(0)
  try
	[ascii(x) for x in ARGS]
	ASCIIStringARGS=[strip(convert(String,x)) for x in ARGS]
	#it might be that parameters such as -v or -e ... where submitted; if so, remove these from the arglist
	while (length(ASCIIStringARGS)>0)&&(length(ASCIIStringARGS[1])>0)&&(ASCIIStringARGS[1][1]=='-')
		shift!(ASCIIStringARGS)
	end
	if ASCIIStringARGS[1]=="0"
		shift!(ASCIIStringARGS) #it may happen that the first entry is "0" (why? because of -p 0?"s
	end
  catch er
	strr="Failed to convert ARGS to ASCIIString"
	println(strr)
	@show ARGS
	@show er
	error(strr)
  end

  if size(ASCIIStringARGS)[1]>0
	   try parse(Int,ASCIIStringARGS[1])
		  #these are the settings if Julia is run in Juno or via the console
		   @show ASCIIStringARGS
		   warn("Using default input parameters!\r\n\r\n")
		   println("$(fallback_inputSettingsFile)")
		   println("$(fallback_inputDataFile)")
		   settingsFilename=fallback_inputSettingsFile
		   dataFilename=fallback_inputDataFile
		   datafolder=convert(String,fallback_datafolder)
		   outfileStringOnly=fallback_outfileStringOnly
		   outfilename=fallback_outputFilename
	catch
      #these are the settings if Julia is run through the SAS script
	  #this_file=ASCIIStringARGS[1]
	  settingsFilename=ASCIIStringARGS[1]
	  dataFilename=ASCIIStringARGS[2]
	  datafolder=ASCIIStringARGS[1][1:end-search(reverse(ASCIIStringARGS[1]),'\\')]
      outfileStringOnly=ASCIIStringARGS[3]
	  outfilename=convert(String,string(datafolder,"\\",outfileStringOnly))
    end
  else
    #these are the settings if Julia is run in Juno or via the console
     warn("Using default input parameters!\r\n\r\n")
	 println("$(fallback_inputSettingsFile)")
	 println("$(fallback_inputDataFile)")
	 settingsFilename=fallback_inputSettingsFile
	 dataFilename=fallback_inputDataFile
	 datafolder=convert(String,fallback_datafolder)
	 outfileStringOnly=fallback_outfileStringOnly
     outfilename=fallback_outputFilename
  end
  if datafolder[end]!='\\'
	datafolder=string(datafolder,'\\') #attach backslash at the end
  end

  return tuple([convert(String,x) for x in [settingsFilename,dataFilename,datafolder,outfilename,outfileStringOnly]]...)
  #return settingsFilename,dataFilename,datafolder,outfilename,outfileStringOnly
end

function reversals(vec;order=true)
  res=BitArray(length(vec))
  res[1]=0
  if order
    for i=2:length(vec)
      res[i]=vec[i]<vec[i-1]
    end
  else
    for i=2:length(vec)
      res[i]=vec[i]>vec[i-1]
    end
  end
  return res
end

function enforce_monotonicity!(vec::Array{T}) where T<:Number
  x=vec[1]
  i=j=0
  @inbounds while i<length(vec)
    while vec[j+1]<x
      vec[j+1]=x
      j+=1
      if j>=length(vec);break;end;
    end
     i+=1
     j=i
     x=vec[i]
  end
  return nothing
end

function lowessSmoothVector!(estimatedRatioPerScore::Array{Float64,1},span::Float64) #,linearregressionfn::Function)
  #Perform LOWESS/LOESS with degree 1 / weighted linear regression  estimatedRatioPerScore=copy(rawObservedRatioPerScore)
  @assert span>0.0
  #warn("it matters if we start at the bottom or top end! which ones is better?")
  nscores=size(estimatedRatioPerScore,1)
  @assert nscores>2 #for 2 or less scores smoothing is not meaningful
  intSpan=convert(Int,ceil((nscores)*span))
  intSpan=max(intSpan,2) #span should be at least two
  smoothed=intercept=slope=0.0

  #bottom end
  thispoint=0
  firstidx=max(1,thispoint-intSpan+1)
  lastidx=max(firstidx,min(nscores,thispoint+intSpan-1))
  x=float(collect(firstidx:lastidx))
  thisRange=copy(estimatedRatioPerScore[firstidx:lastidx]) #todo, determine if we need the copy command (do we modify current range later on ?)
  w=1.0.-((x.-thispoint)./intSpan).^2
  while thispoint<intSpan
    thispoint+=1
    push!(thisRange,estimatedRatioPerScore[thispoint+intSpan-1])
    push!(x,x[end]+1.0)
    newweight=1.0.-((x[1]-thispoint)/intSpan)^2
    unshift!(w,newweight)
    intercept,slope=mylinreg(x,thisRange,w)
	smoothed=float(thispoint)*slope+intercept
	estimatedRatioPerScore[thispoint]=smoothed
	#thisRange[thispoint]=smoothed	#this seems to make it worse
  end

  #middle part
  #thisRangeMid=length(thisRange)+1
  #wVec=FrequencyWeights(w)
  while thispoint<=nscores-intSpan
    thispoint+=1
    next=estimatedRatioPerScore[thispoint+intSpan-1]
    push!(thisRange,next)
    shift!(thisRange)
    push!(x,x[end]+1.0)
    shift!(x)
    intercept::Float64,slope::Float64=mylinreg(x,thisRange,w) #todo tbd, this line takes quite some time as it is used very often, this can certainly be improved!
    smoothed=float(thispoint)*slope+intercept
	estimatedRatioPerScore[thispoint]=smoothed
	#thisRange[thisRangeMid]=smoothed #this seems to make it worse
  end

  #Top end: consider smaller intSpan
  while thispoint<nscores
    thispoint+=1
    #remove the last element of each vector
    shift!(thisRange)
    shift!(x)
    pop!(w)
    intercept,slope=mylinreg(x,thisRange,w)
	smoothed=thispoint*slope+intercept
	estimatedRatioPerScore[thispoint]=smoothed
	#thisRange[end-(nscores-thispoint)]=smoothed #this seems to make it worse
  end

  return nothing
end

function mylinreg(x::Array{Float64,1},y::Array{Float64,1},w::Array{Float64,1})
	meanx=0.0
	meany=0.0
	meanxy=0.0
	meanx2=0.0
	wsum=0.0
	@inbounds for i=1:length(x)
		 xi=x[i]
		 wi=w[i]
		 yi=y[i]
		xwi=xi*wi
		meanx+=xwi
		meany+=yi*wi
		meanxy+=xwi*yi
		meanx2+=xwi*xi
		wsum+=wi
	end
	
	meanx/=wsum
	meany/=wsum
	meanxy/=wsum
	meanx2/=wsum

	beta::Float64 = (meanx2 == meanx * meany) ? 0.0 : (meanxy - meanx * meany) / (meanx2 - meanx * meanx)
	alpha::Float64 = meany - beta * meanx

	return alpha::Float64,beta::Float64
end
#=
@inline function mylinreg(x::Array{Float64,1},y::Array{Float64,1},w::Array{Float64,1})
	return mylinreg(x,y,FrequencyWeights(w))
end
=#
function mylinreg(x::Array{Float64,1},y::Array{Float64,1},w::W) where {W <: AbstractWeights}	
	meanx=0.0
	meany=0.0
	meanxy=0.0
	meanx2=0.0
	@inbounds for i=1:length(x)
		 xi=x[i]
		 wi=w[i]
		 yi=y[i]
		xwi=xi*wi
		meanx+=xwi
		meany+=yi*wi
		meanxy+=xwi*yi
		meanx2+=xwi*xi
	end
	meanx/=w.sum
	meany/=w.sum
	meanxy/=w.sum
	meanx2/=w.sum

	beta::Float64 = (meanx2 == meanx * meany) ? 0.0 : (meanxy - meanx * meany) / (meanx2 - meanx * meanx)
	alpha::Float64 = meany - beta * meanx

	return alpha::Float64,beta::Float64
end


function update_and_derive_scores!(scores::Vector{Int},maxRawRelativityPerScoreSorted::Array{Float64,1},relativities,idx::Vector{Int}) #::Array{Float64,1})	  
	  sz=size(maxRawRelativityPerScoreSorted,1) #it can happen that the rawrel is larger than what we observe in training	  	  
	  for i in idx					
        @inbounds scores[i]=min(searchsortedfirst(maxRawRelativityPerScoreSorted,relativities[i]),sz)
	  end	  
    return nothing
end


function update_and_derive_scores!(scores::Vector{Int},maxRawRelativityPerScoreSorted::Array{T,1},relativities::Vector{T}) where T <: Number	  
	sz=size(maxRawRelativityPerScoreSorted,1) #it can happen that the rawrel is larger than what we observe in training	  	  
	for i=1:length(scores)
	  @inbounds scores[i]=min(searchsortedfirst(maxRawRelativityPerScoreSorted,relativities[i]),sz)
	end	  
  return nothing
end

function sumByIncreasingIndex(x,list) #{T<:Number}(x::Array{T,1},list::Array{Int,1}) #think of a better name for this function
	res=zeros(eltype(x),length(list))
	current=1
	for k=1:length(list)
		endpoint=list[k]
		#endpoint can be 0 or m>0 here
		while current<=endpoint
			res[k]+=x[current]
			current+=1
		end
	end
	return res
end


function sumByIncreasingIndex(x,y,list) #{T<:Number}(x::Array{T,1},y::Array{T,1},list::Array{Int,1}) #think of a better name for this function
	res=zeros(eltype(x),length(list))
	res2=zeros(eltype(y),length(list))
	current=1
	for k=1:length(list)
		endpoint=list[k]
		#endpoint can be 0 or m>0 here
		while current<=endpoint
			res[k]+=x[current]
			res2[k]+=y[current]
			current+=1
		end
	end
	return res,res2
end

function deriveRowCountAndWeight(weight_srt,raw_rel_srt,uniqueRelativitiesSorted)
	rowCountPerUniqueRelativity=Vector{Int}(length(uniqueRelativitiesSorted))
	for i=1:length(uniqueRelativitiesSorted)
		@inbounds ui::Float64=uniqueRelativitiesSorted[i]
		@inbounds rowCountPerUniqueRelativity[i]=searchsortedlast(raw_rel_srt,ui)
	end
	weightPerUniqueRelativity=sumByIncreasingIndex(weight_srt,rowCountPerUniqueRelativity)
	cumulativeToIncremental!(rowCountPerUniqueRelativity)
return rowCountPerUniqueRelativity,weightPerUniqueRelativity
end

function sum_up_ratios!(rawObservedRatioPerScore,aggregatedModelledRatioPerScore,idx::Int,numerator_srt,numeratorEstimatedPerRow_srt,denominator_srt,previdx::Int,nextidx::Int)
	nest=0.0
	n=0.0
	d=0.0	
	@inbounds for i=previdx:nextidx
		nest += numeratorEstimatedPerRow_srt[i]
		n += numerator_srt[i]
		d += denominator_srt[i]
	end
	aggregatedModelledRatioPerScore[idx]=nest/d
	rawObservedRatioPerScore[idx]=n/d
	return nothing
end

function constructANDderiveScores!(trnidx::Vector{Int},validx::Vector{Int},sortvec_reused::Vector{Int},estimatedRatio::Array{Float64,1},relativities::Array{Float64,1},actualNumerator::Array{Float64,1},denominator::Array{Float64,1},weight::Array{Float64,1},meanobservedvalue::Float64,currentIteration::Int,sett::ModelSettings)
	deriveFitPerScoreFromObservedRatios=sett.deriveFitPerScoreFromObservedRatios
	nscores=sett.nscores
	smooth=sett.smoothEstimates
	print_details=sett.print_details
	boolSmoothingEnabled=!(smooth=="0")
	aggregatedModelledRatioPerScore,maxRawRelativityPerScoreSorted,MAPPINGSmoothedEstimatePerScore,vectorWeightPerScore,obsPerScore,rawObservedRatioPerScore,numPerScore,denomPerScore=constructScores!(deriveFitPerScoreFromObservedRatios,trnidx,validx,sortvec_reused,estimatedRatio,relativities,actualNumerator,denominator,weight,meanobservedvalue,nscores,print_details,boolSmoothingEnabled=boolSmoothingEnabled)	
return maxRawRelativityPerScoreSorted,MAPPINGSmoothedEstimatePerScore,vectorWeightPerScore,obsPerScore,rawObservedRatioPerScore,numPerScore,denomPerScore,nscores
end

function constructScores!(deriveFitPerScoreFromObservedRatios::Bool,trnidx::Vector{Int},validx::Vector{Int},srt::Vector{Int},estimatedRatioPerRow_full_vec::Array{Float64,1},raw_rel_full_vec::Array{Float64,1},numerator_full_vec::Array{Float64,1},denominator_full_vec::Array{Float64,1},weight_full_vec::Array{Float64,1},meanobservedvalue::Float64,nscores::Int,print_details::Bool;boolSmoothingEnabled::Bool=true)
	#@assert nscores>=nscoresPotentiallyReduced
	maxItersForSmoothing=0
    #todo,tbd this can probably be done more efficiently
	raw_rel=view(raw_rel_full_vec,trnidx)	
	weight=view(weight_full_vec,trnidx)
    numerator=view(numerator_full_vec,trnidx)
    denominator=view(denominator_full_vec,trnidx)
    estimatedRatioPerRow=view(estimatedRatioPerRow_full_vec,trnidx)
	
	#warn("remove this");writecsv("C:\\temp\\raw_rel.csv",raw_rel)
	#srt::Vector{Int}=sortperm(raw_rel,alg=QuickSort) #we could use RadixSort here, it uses about 40m of memory for 5m obs, AND I am unsure how easy it is to get the permutation (compared to sorting in place). But the algorithm is about 50% faster for random floats than QuickSort.
	#sortperm!(srt,raw_rel,alg=QuickSort)
	my_sortperm_you_should_regularily_check_if_sortperm_in_base_has_become_more_efficient!(srt,raw_rel) #we could use RadixSort here, it uses about 40m of memory for 5m obs, AND I am unsure how easy it is to get the permutation (compared to sorting in place). But the algorithm is about 50% faster for random floats than QuickSort.
	obs=length(raw_rel)
	raw_rel_srt=view(raw_rel,srt)
	#@show unique(raw_rel_srt) #this is moderated and has 'little' spread for small mf
	weight_srt=view(weight,srt)
	estimatedRatioPerRow_srt=view(estimatedRatioPerRow,srt)
	#@show unique(estimatedRatioPerRow_srt) #this is moderated and has 'little' spread for small mf
	denominator_srt=view(denominator,srt)
	numerator_srt=view(numerator,srt)

	#rowCountPerUniqueRelativity::Vector{Int64},weightPerUniqueRelativity::Vector{Float64}=deriveRowCountAndWeight(weight_srt,raw_rel_srt,uniqueRelativitiesSorted)
	weightPerUniqueRelativity=weight_srt
	
	numeratorEstimatedPerRow_srt=estimatedRatioPerRow_srt.*denominator_srt	
	wtot=sum(weight_srt)
	wperscore=wtot/nscores
	scoreEndPoints=zeros(Int,nscores)
	scoreEndPoints[end]=size(weight_srt,1)
	vectorWeightPerScore=Array{Float64}(nscores)

	# derive_scores_main_aggregation_step!(nscores,wperscore,relativitiesSorted,weight_srt,scoreEndPoints,vectorWeightPerScore)
	nscoresPotentiallyReducedTWOTimes=derive_scores_main_aggregation_step!(nscores,wperscore,raw_rel_srt,weight_srt,scoreEndPoints,vectorWeightPerScore)
	
	obsPerScore=copy(scoreEndPoints)
	cumulativeToIncremental!(obsPerScore)
	cumulativeToIncremental!(vectorWeightPerScore)

	#aggregate scores
	numPerScore,denomPerScore,maxRawRelativityPerScoreSorted,aggregatedModelledRatioPerScore,rawObservedRatioPerScore=aggregate_values_per_score(nscoresPotentiallyReducedTWOTimes,scoreEndPoints,raw_rel_srt,numerator_srt,denominator_srt,obs,numeratorEstimatedPerRow_srt)
	#smooth scores		
	#deriveFitPerScoreFromObservedRatios::Bool (if this is true then the fit per score will be based on the observed ratio per score (instead of the fitted ratio per score)) . Note that this does not influence the structure of the tree, but only the fitted ratio per score (and scoreband)	
	referenceVec = deriveFitPerScoreFromObservedRatios ? rawObservedRatioPerScore : aggregatedModelledRatioPerScore
	estimatedRatioPerScore=smooth_scores(referenceVec,print_details,boolSmoothingEnabled)		
	#insert gaps if necessary
	insert_gaps_to_scores!(wtot,nscores,nscoresPotentiallyReducedTWOTimes,aggregatedModelledRatioPerScore,maxRawRelativityPerScoreSorted,estimatedRatioPerScore,vectorWeightPerScore,obsPerScore,rawObservedRatioPerScore,numPerScore,denomPerScore)

	return aggregatedModelledRatioPerScore,maxRawRelativityPerScoreSorted,estimatedRatioPerScore,vectorWeightPerScore,obsPerScore,rawObservedRatioPerScore,numPerScore,denomPerScore
end

#(nscores,wperscore,relativitiesSorted,weight_srt,cumulativeNumberOfDistinctRawRelativitesPerScore,scoreEndPoints,vectorWeightPerScore)
function derive_scores_main_aggregation_step!(nscores,wperscore,relativitiesSorted,weight_srt,scoreEndPoints,vectorWeightPerScore)
	#todo/tbd this can be simplified if length(uniqueRelativitiesSorted)<sett.nscores then we do not need to do any aggregation at all.

#NOTE: (todo/tbd improve this) there are different ways to aggregate the relativities to scores. the approach below is a bottom up approach
#issues can arise when there is a mass (of weight) for a certain releativity towards the end of the scores
#because the approach below assigns at least one relativity to each score "one can run out of data" before all scores are filled up (in certain degenerate cases)
#an alternative would be a "binary" split starting in the middle

#note: if the estimates and the weight are degenerated then the vector scoreEndPoints can have the same value multiple times (i.e. it is constant for a few indices)
uniqueRelativitiesSortedCounter=0
current_cumulW=0.0
startPoint=1
# assert size(rowCountPerUniqueRelativity)==size(weightPerUniqueRelativity) #this is to avoid errors due to the inbounds macro below
stillRunning=true
#wi=weightPerUniqueRelativity[1]
wi=0.0
current_score=0
current_relativity=0.0
current_idx=0 #rowCountPerUniqueRelativity[1]

while stillRunning
current_score+=1
  	while current_idx<length(weight_srt) 
		current_idx+=1
		last_relativity=current_relativity
		@inbounds current_relativity=relativitiesSorted[current_idx]
		@inbounds wi::Float64=weight_srt[current_idx]::Float64 
		current_cumulW::Float64 += wi::Float64
		if current_cumulW>wperscore*current_score
			#-> weight per score is reached
			#we now need to ensure that the 'cut' is made at a point where the relativity changes (otherwise it will be impossible to define the scores in relation to the (moderated) relativities)
			while current_relativity==last_relativity
				current_idx+=1
				if current_idx>length(relativitiesSorted)
					break
				end
				last_relativity=current_relativity				
				@inbounds current_relativity=relativitiesSorted[current_idx]
				@inbounds wi=weight_srt[current_idx]::Float64
				current_cumulW += wi::Float64		
			end
			break
		end
 	end
	vectorWeightPerScore[current_score]=current_cumulW
	scoreEndPoints[current_score]=current_idx
	stillRunning=current_score<length(scoreEndPoints)
end

nscoresPotentiallyReducedTWOTimes::Int=0
if (length(scoreEndPoints)>1)&&(scoreEndPoints[end-1]>=scoreEndPoints[end])
#warn("Failed to uniformly distribute the scores! Scores will be degenerated. \r\nThis usually indicates that: (i) there are too many scores compared to the granularity of the model or (ii) the data is 'degenerated' (i.e. there is a mass of exposure for certain risk details)")
#Reduce the number of Scores again
	thismax=maximum(scoreEndPoints)
	endlocation::Int=searchsortedfirst(scoreEndPoints,convert(eltype(scoreEndPoints),thismax))
	resize!(scoreEndPoints,endlocation) #drop the tail end
	resize!(vectorWeightPerScore,endlocation) #drop the tail end
	nscoresPotentiallyReducedTWOTimes=length(scoreEndPoints)
else
	nscoresPotentiallyReducedTWOTimes=nscores::Int
end

scoreEndPoints[end]=min(length(weight_srt),scoreEndPoints[end])

return nscoresPotentiallyReducedTWOTimes::Int
end

function derive_scores_main_aggregation_step_old!(nscoresPotentiallyReduced,wperscore,uniqueRelativitiesSorted,weightPerUniqueRelativity,rowCountPerUniqueRelativity,cumulativeNumberOfDistinctRawRelativitesPerScore,scoreEndPoints,vectorWeightPerScore)
	#todo/tbd this can be simplified if length(uniqueRelativitiesSorted)<sett.nscores then we do not need to do any aggregation at all.

#NOTE: (todo/tbd improve this) there are different ways to aggregate the relativities to scores. the approach below is a bottom up approach
#issues can arise when there is a mass (of weight) for a certain releativity towards the end of the scores
#because the approach below assigns at least one relativity to each score "one can run out of data" before all scores are filled up (in certain degenerate cases)
#an alternative would be a "binary" split starting in the middle

#note: if the estimates and the weight are degenerated then the vector scoreEndPoints can have the same value multiple times (i.e. it is constant for a few indices)
uniqueRelativitiesSortedCounter=current_score=current_idx=0;current_cumulW=0.0;startPoint=1
# assert size(rowCountPerUniqueRelativity)==size(weightPerUniqueRelativity) #this is to avoid errors due to the inbounds macro below
stillRunning=true
weightOfThisRelativity=weightPerUniqueRelativity[1]
thisrowcount=rowCountPerUniqueRelativity[1]
while stillRunning
current_score+=1
  while uniqueRelativitiesSortedCounter<size(uniqueRelativitiesSorted,1)
	uniqueRelativitiesSortedCounter+=1
	@inbounds thisrowcount=rowCountPerUniqueRelativity[uniqueRelativitiesSortedCounter]
	current_idx+=thisrowcount
	@inbounds weightOfThisRelativity::Float64=weightPerUniqueRelativity[uniqueRelativitiesSortedCounter]::Float64 #um(view(weight_srt,startPoint:current_idx))
	current_cumulW::Float64 += weightOfThisRelativity::Float64
	if current_cumulW>wperscore*current_score
	  #decrease index by 1 if this puts us closer to the "expected" exposure for this score
	  if (current_idx-thisrowcount>startPoint) && (abs(current_cumulW-wperscore*current_score)>abs(current_cumulW-weightOfThisRelativity-wperscore*current_score))
		uniqueRelativitiesSortedCounter-=1
		current_cumulW-=weightOfThisRelativity
		current_idx-=thisrowcount
	  end
	  break
	end
  end
	cumulativeNumberOfDistinctRawRelativitesPerScore[current_score]=uniqueRelativitiesSortedCounter
	vectorWeightPerScore[current_score]=current_cumulW
	scoreEndPoints[current_score]=current_idx
	startPoint=current_idx+1
	stillRunning=current_score<length(scoreEndPoints)
end

nscoresPotentiallyReducedTWOTimes::Int=0
if (length(scoreEndPoints)>1)&&(scoreEndPoints[end-1]>=scoreEndPoints[end])
#warn("Failed to uniformly distribute the scores! Scores will be degenerated. \r\nThis usually indicates that: (i) there are too many scores compared to the granularity of the model or (ii) the data is 'degenerated' (i.e. there is a mass of exposure for certain risk details)")
#Reduce the number of Scores again
	thismax=maximum(scoreEndPoints)
	endlocation::Int=searchsortedfirst(scoreEndPoints,convert(eltype(scoreEndPoints),thismax))
	resize!(scoreEndPoints,endlocation) #drop the tail end
	resize!(vectorWeightPerScore,endlocation) #drop the tail end
	nscoresPotentiallyReducedTWOTimes=length(scoreEndPoints)
else
	nscoresPotentiallyReducedTWOTimes=nscoresPotentiallyReduced::Int
end

return nscoresPotentiallyReducedTWOTimes::Int
end

function aggregate_values_per_score(nscoresPotentiallyReducedTWOTimes,scoreEndPoints,raw_rel_srt,numerator_srt,denominator_srt,obs,numeratorEstimatedPerRow_srt)
	#NOTE: we consider the mean of two values here, to avoid floating point issues (there are no issues within Julia, but the C# Module may lead to slightly different results)
	maxRawRelativityPerScoreSorted=Float64[(raw_rel_srt[scoreEndPoints[i]]+raw_rel_srt[min(obs,scoreEndPoints[i]+1)])/2 for i=1:size(scoreEndPoints,1)]
	numPerScore::Vector{Float64},denomPerScore::Vector{Float64}=sumByIncreasingIndex(numerator_srt,denominator_srt,scoreEndPoints)
	#I think the weight should not influence the averaging process here (as it was already relevant to derive the raw_estimated_relativities during the construction of the trees)
	#raw_num_est_srt=raw_rel_srt.*denominator_srt*meanobservedvalue
	rawObservedRatioPerScore=zeros(Float64,nscoresPotentiallyReducedTWOTimes)
	aggregatedModelledRatioPerScore=zeros(Float64,nscoresPotentiallyReducedTWOTimes)
	#aggregatedModelledRatioPerScore_new=zeros(Float64,nscoresPotentiallyReducedTWOTimes)
	previdx=1
	for i=1:nscoresPotentiallyReducedTWOTimes
		nextidx=scoreEndPoints[i]
		if nextidx>previdx
		  #todo/tbd check with somebody if this approach is meaningful (are there better alternatives to derive a smoothed estimate?)
		  #rawObservedRatioPerScore[i]=			sum(view(numerator_srt,previdx:nextidx))/sum(view(denominator_srt,previdx:nextidx))
		  #aggregatedModelledRatioPerScore[i]=	sum(view(numeratorEstimatedPerRow_srt,previdx:nextidx))/sum(view(denominator_srt,previdx:nextidx)) #todo/tbd someone should check whether we do not have any redundancy here
		  sum_up_ratios!(rawObservedRatioPerScore,aggregatedModelledRatioPerScore,i,numerator_srt,numeratorEstimatedPerRow_srt,denominator_srt,previdx,nextidx)		 
		  #dd= aggregatedModelledRatioPerScore_new.-aggregatedModelledRatioPerScore
		  # dd[1:min(length(dd),200)]
		  # @assert maximum(abs.(dd))<1e-11 #aggregatedModelledRatioPerScore_new,aggregatedModelledRatioPerScore)
		end
		previdx=scoreEndPoints[i]
	end
	#@show 1,extrema(numeratorEstimatedPerRow_srt.-numerator_srt) #these are different
	#@show 2,extrema(aggregatedModelledRatioPerScore.-rawObservedRatioPerScore) #these are different
	return numPerScore,denomPerScore,maxRawRelativityPerScoreSorted,aggregatedModelledRatioPerScore,rawObservedRatioPerScore
end

function smooth_scores(rawObservedRatioPerScore,print_details,boolSmoothingEnabled)	
	#notably the smoothing here is performed regardless of the value of boolSmoothingEnabled
	#however the Excel (and other stats) will show different figures depending on the value of sett.smoothestimates ("0" or "1")
		maxItersForSmoothing=0
		#smooth the scores
		estimatedRatioPerScore=copy(rawObservedRatioPerScore)
	
		if print_details&&(length(estimatedRatioPerScore)<=2)
			printover("length(estimatedRatioPerScore) is 2 or less! No smoothing performed.")
		else
			span=0.0
			if length(estimatedRatioPerScore)>1000
				span=50/length(estimatedRatioPerScore)
			else
				span=0.05
			end
			
			if boolSmoothingEnabled
				maxItersForSmoothing=30
			else
				maxItersForSmoothing=2
			end
					i=0
					sumrev=typemax(Int)
					maxiterNotReached=true
					revs=max(10,convert(Int,round(size(estimatedRatioPerScore,1)/100)))					
					#arbitrary number of 10 passes which is always done (if there are still reversions)
					while (i<min(maxiterNotReached,10)&&sumrev>0)||((sumrev>revs)&&maxiterNotReached) #todo/tbd simplify this: What exactly do we want to achieve here!? what is our approach?
						i+=1
						maxiterNotReached=i<maxItersForSmoothing
						lowessSmoothVector!(estimatedRatioPerScore,span) #,mylinreg)
						rev=reversals(estimatedRatioPerScore)
						sumrev=sum(rev)
					end
					#maxiterNotReached ? println("Smoothing of estimates finished after $(i) recursions.") :  warn("Smoothing process aborted after $(i) iterations. There were still $(sumrev) reversals.")
					if print_details&&boolSmoothingEnabled #if smoothing is disabled this will not be of interest to the modeller
						maxiterNotReached ? nothing :  @info "Smoothing stopped after $(i) passes: $(sumrev) reversals remained."
					end
					enforce_monotonicity!(estimatedRatioPerScore)
		end
			return estimatedRatioPerScore
		end

		
function insert_gaps_to_scores!(wtot,nscores::Int,nscoresPotentiallyReducedTWOTimes::Int,aggregatedModelledRatioPerScore,maxRawRelativityPerScoreSorted,estimatedRatioPerScore,vectorWeightPerScore,obsPerScore,rawObservedRatioPerScore,numPerScore,denomPerScore)
	#NOTE: here we blow up the scores such that we have "nscores" (e.g. 1000) scores. This is only relevant when nscores>length(uniqueRelativitiesSorted) (=nscoresPotentiallyReduced) which is the case for the first few iterations or somewhat degenerated data sets (todo tbd: check this! is it meaningful?)
	nscoresToBeIncreased=nscoresPotentiallyReducedTWOTimes
	cumulw=0.0
	loopVariable=0
	wperscoreForActualNscores=wtot/nscores
	while nscoresToBeIncreased<nscores
		loopVariable+=1
		cumulw+=vectorWeightPerScore[loopVariable]
		while (cumulw>wperscoreForActualNscores*loopVariable)&&(nscoresToBeIncreased<nscores)
			#insert "gaps" into score vectors
			insert!(maxRawRelativityPerScoreSorted,loopVariable,maxRawRelativityPerScoreSorted[loopVariable])
			insert!(estimatedRatioPerScore,loopVariable,estimatedRatioPerScore[loopVariable])
			insert!(vectorWeightPerScore,loopVariable,0.0)
			insert!(obsPerScore,loopVariable,0.0)
			insert!(rawObservedRatioPerScore,loopVariable,rawObservedRatioPerScore[loopVariable])
			insert!(numPerScore,loopVariable,0.0)
			insert!(denomPerScore,loopVariable,0.0)
			nscoresToBeIncreased+=1
			loopVariable+=1
		end
	end
	return nothing
end

function cumulativeToIncremental!(x::Vector{T}) where {T <: Number}
	for i=length(x):-1:2
		x[i]-=x[i-1]
	end
	return nothing
end

function incrementalToCumulative!(x::Vector{T}) where {T <: Number}
	for i=2:length(x)
		x[i]+=x[i-1]
	end
	return nothing
end

#univariate graphs output
function createPredictorData(idx::Vector{Int},nameOfpredictorsSheet,mappings,candMatWOMaxValues,sett,scoresfull,actualNumeratorfull,denominatorfull,weightfull,finalEstimateForChartsfull,featuresfull)
	#Create tables (=data) and  Charts
		scores=view(scoresfull,idx)
		actualNumerator=view(actualNumeratorfull,idx)
		denominator=view(denominatorfull,idx)
		weight=view(weightfull,idx)
		finalEstimateForCharts=view(finalEstimateForChartsfull,idx)
		features=view(featuresfull,idx)

		colnames=String["Variable Name" "Weight" "Numerator" "Denominator" "Observed Ratio" "Estimate" "Difference" "Average  Score"]
		predictorsData=repmat([""],1,length(colnames))
		numeratorEst=denominator.*finalEstimateForCharts
		@assert length(numeratorEst)==length(weight)==length(denominator)==length(scores)==length(actualNumerator) #this is essential: the inbounds macro will lead to a serious error (within aggregate_data) if the vectors have unequal length
		#charts
		nfeat=size(features,2)
		predictorCharts=Array{Chart}(2*(nfeat)) #we have two charts for each feature
		chartscol="J"
		chartscol2="AJ"
		addchartrow=21
		headerrow=2
		currentchartrow=copy(headerrow)
		local thischart,lastchartrow
		#loop over features
		for i=1:nfeat #length(numfeatures)
			lastchartrow=currentchartrow			
			vname=deepcopy(sett.df_name_vector[i])
			feat=features[i]
			thischart=defineUnivariateChart(nameOfpredictorsSheet,nameOfpredictorsSheet,convert(typeof(nameOfpredictorsSheet),string(chartscol,currentchartrow)),vname,length(levels(feat)),1,8,2,headerrow)
			predictorCharts[i]=deepcopy(thischart)
			thischart=defineUnivariateChartWith2Lines(nameOfpredictorsSheet,nameOfpredictorsSheet,convert(typeof(nameOfpredictorsSheet),string(chartscol2,currentchartrow)),vname,length(levels(feat)),1,5,6,2,headerrow)
			predictorCharts[nfeat+i]=deepcopy(thischart)
			#todo,tbd maybe this can be replaced with feat.parent.pool
		    mp = features[i].pool  #this is the FULL pool, in the selected data only a subset might exist.		
			tbl=addPredictorData(mp,colnames,sett,scores,numeratorEst,actualNumerator,denominator,weight,finalEstimateForCharts,feat,i)
			predictorsData=vcat(predictorsData,tbl)
			currentchartrow=max(lastchartrow+addchartrow,size(predictorsData,1)+1)
			headerrow+=size(tbl,1)
		end		
	return predictorsData,predictorCharts
end

function replace_in_vector!(arr::Array{T,1},orig::T,new::T) where {T}
	for i=1:length(arr)
		if arr[i]==orig
			arr[i]=new
		end
	end
	return nothing
end

function addPredictorData(listOfValues,colnames,sett::ModelSettings,scores,numeratorEst,numerator,denominator,weight,finalEstimateForCharts,f,varnum::Int64)
	thisdata=deepcopy(colnames)
	thisdata[1]=deepcopy(sett.df_name_vector[varnum])
	valuelist,cnt,sumscores,sumnumeratorEst,sumnumerator,sumdenominator,sumweight=aggregate_data(f,scores,numeratorEst,numerator,denominator,weight)
	#@show typeof(f)
	#@show f.pool
	#@show valuelist
	if sum(isnan.(sumscores))!=0
		@show thisdata[1]
		@show f.ids
		@show f.values
		@show cnt
		@show sumscores
		@show sumnumeratorEst
		@show sumnumerator
		@show sumdenominator
		@show sumweight
		@assert false "Critical Error: Sumscores contains NaN values" #probably the @inbounds macro led to an issue when aggregating the data
	end
	#replace 0 entries with 1s (otherwise we divde by zero)
	cntmod=deepcopy(cnt)
	replace_in_vector!(cntmod,zero(eltype(cntmod)),one(eltype(cntmod)))
	avgscores=Int[convert(Int,x) for x in round.(sumscores./cntmod)]
	#listOfValues may contain values which do not occurr in the training data
	#therefore we need to slightly adjust the output in some instances

	local tbl	
	try
		observed_ratio=sumnumerator./sumdenominator
		fitted_ratio=sumnumeratorEst./sumdenominator
		tbl=hcat(valuelist,sumweight,sumnumerator,sumdenominator,observed_ratio,fitted_ratio,observed_ratio-fitted_ratio,avgscores)
	catch
		@show thisdata[1]
		try
			@show f.refs
		catch
			f.parent.refs[f.indices[1]]
		end
		try
			@show f.pool
		catch
			@show f.parent.pool
		end
		@show valuelist
		#@show size(f.ids)
		#@show size(f.values)
		@show size(valuelist)
		@show size(sumweight)
		@show size(sumnumerator)
		@show size(sumnumeratorEst)
		@show size(sumdenominator)
		@show size(avgscores)
		@show size(sumnumeratorEst./sumdenominator)
		save("c:\\temp\\mi.jld2","f",f,"listOfValues"	,listOfValues,	"colnames"	,colnames,	"sett"	,sett,	"scores"	,scores,	"numeratorEst"	,numeratorEst,	"numerator"	,numerator,	"denominator"	,denominator,	"weight"	,weight,	"finalEstimateForCharts"	,finalEstimateForCharts,	"varnum"	,varnum)
		@assert false "Critical Error when concatenating predictor statistics"
	end
	thisdata=vcat(thisdata,tbl,repmat([""],1,length(colnames)))
	return thisdata
end

function csharp_write_pub_str(fiostream::IOStream,name,value)
	write(fiostream,"public string ",name,'{'," get ",'{'," return ",DoubleQuote,value,DoubleQuote,"; } }\r\n")
	return nothing
end

function numberAsStringForCSharp(x)
	#todo tbd: if we have very small values things will fail here!
	#e.g. string(0.0000000000000000000000000000000000000002329) will become "2.329e-34" or similar, this might fail in c#
	s=replace(replace(string(x),'-',"MINUS"),'.','_')
	return s
end

function vba_get_signature(mappings,df_name_vector,number_of_num_features)
	res=""
	for i=1:number_of_num_features
		thisvarname=df_name_vector[i]		
		res=string(res,thisvarname," As Double,")		
	end
	for i=1:length(mappings)
		thisvarname=df_name_vector[number_of_num_features+i]		
		res=string(res,thisvarname," As String,")
	end	

	str_fn_evaluate_profile="\r\nFunction readDataAndDeriveRawRelativity() as Double\r\n\r\nApplication.Volatile\r\n'Define Variables\r\n"	
	for i=1:number_of_num_features
		thisvarname=df_name_vector[i]		
		str_fn_evaluate_profile=string(str_fn_evaluate_profile,"Dim ",thisvarname," As Double\r\n")		
	end
	for i=1:length(mappings)
		thisvarname=df_name_vector[number_of_num_features+i]		
		str_fn_evaluate_profile=string(str_fn_evaluate_profile,"Dim ",thisvarname," As String\r\n")
	end	

	str_fn_evaluate_profile=string(str_fn_evaluate_profile,"\r\n\r\n'Read Data from named ranges\r\n'NOTE: all strings are converted to uppercase here!\r\n")

	for i=1:number_of_num_features
		thisvarname=df_name_vector[i]		
		str_fn_evaluate_profile=string(str_fn_evaluate_profile,thisvarname," = Range(",DoubleQuote,thisvarname,DoubleQuote,").Value\r\n") #As Double\r\n")		
	end
	for i=1:length(mappings)
		thisvarname=df_name_vector[number_of_num_features+i]		
		str_fn_evaluate_profile=string(str_fn_evaluate_profile,thisvarname," = Ucase(Range(",DoubleQuote,thisvarname,DoubleQuote,").Text)\r\n")
	end	

	callString=""
	for i=1:number_of_num_features
		thisvarname=df_name_vector[i]		
		callString=string(callString,thisvarname,":=",thisvarname,", ")		
	end
	for i=1:length(mappings)
		thisvarname=df_name_vector[number_of_num_features+i]		
		callString=string(callString,thisvarname,":=",thisvarname,", ")
	end	
	#remove last two chars
	callString=callString[1:end-2]
	callString=string("raw_rel = deriveRawRelativity(",callString,")\r\n readDataAndDeriveRawRelativity = raw_rel\r\n\r\nEnd Function\r\n")

# resvba = deriveRawRelativity(GRP_NEW_NCD_NUMERIC:=GRP_NEW_NCD_NUMERIC, NEW_NCD_LAST_NUMERIC:=NEW_NCD_LAST_NUMERIC, GRP_CAR_AGE_NUMERIC:=GRP_CAR_AGE_NUMERIC, INC_YEAR:=INC_YEAR, SEAT:=SEAT, CAR_PRICE:=CAR_PRICE, AGE:=AGE, TONNAGE:=TONNAGE, HIGH_RISK_IND:=HIGH_RISK_IND, AIRBAG_JY:=AIRBAG_JY, DISP:=DISP, PROD_CODE:=PROD_CODE, BRAND_NAME_JY:=BRAND_NAME_JY, ORG_CLASS:=ORG_CLASS, CUST_TYPE:=CUST_TYPE, CHANNEL_GROUP:=CHANNEL_GROUP, NCD_BFD:=NCD_BFD, IS_NEW:=IS_NEW, MULTI_COV:=MULTI_COV, TRUCK_TYPE:=TRUCK_TYPE, SERIES_LEVEL:=SERIES_LEVEL, GRP_RENEW:=GRP_RENEW, GRP_SEX:=GRP_SEX)
	str_fn_evaluate_profile=string(str_fn_evaluate_profile,"\r\n\r\n",callString,"\r\n\r\n")

	#remove last comma from res
	return res[1:end-1],str_fn_evaluate_profile 
end

function write_and_create_boollist_char_vba(fiostream::IOStream,mappings::Array{Array{String,1}},df_name_vector::Array{String,1},number_of_num_features::Int,boolCharVarsUsedByModel::Array{Array{Bool,1},1})
	#note todo tbd check this: here we only consider the values which we observe in the training data
	#of course we would also know the values from the validation data, but that should not matter here, since any new data which the module
	#will see, might have even other/new values in it.
	#Does that make sense? Or should we include the values from the val data in this list?
		boollist=Array{Array{String,1}}(0)
		for i=1:length(mappings)
			thisvarname=df_name_vector[number_of_num_features+i]
			thisboollist=Array{String}(0)
			for j=1:length(mappings[i])
					str=create_custom_string(mappings[i][j]) #this can be any String, let us see what C# can handle
					#if C# isunable to handle the UTF8Strings, we need to convert it to some ASCIIString, while ensuring that we still have UNIQUENESS, todo /tbd check this
					thisbool=convert(String,string(thisvarname,"_",str))
					while in(thisbool,thisboollist) #need to ensure that the list is unique!
						thisbool=convert(String,string(thisbool,"_2"))
					end
					push!(thisboollist,thisbool)
					if boolCharVarsUsedByModel[i][j]
						write(fiostream,"Dim ",thisbool," as Boolean: ",thisbool," = False\r\n")
					end
			end
			@assert unique(thisboollist)==thisboollist
			push!(boollist,copy(thisboollist))
		end
		write(fiostream,"\r\n")
		return boollist
end

function write_and_create_boollist_num_vba(fiostream::IOStream,df_name_vector::Array{String,1},number_of_num_features::Int,candMatWOMaxValues::Array{Array{Float64,1},1},boolNumVarsUsedByModel::Array{Array{Bool,1},1})
#NOTE, here we assume that the numerical variables are "first" in the df_name_vector
#todo: only define booleans which are actually used by the model! this will keep the code smaller and likely more efficient
boollist=Array{Array{String,1}}(0)
	for i=1:number_of_num_features
		thisvarname=df_name_vector[i]
		thisboollist=Array{String}(0)
		for j=1:length(candMatWOMaxValues[i])
				thisbool=convert(String,string(thisvarname,"_isLE_",numberAsStringForCSharp(candMatWOMaxValues[i][j])))
				while in(thisbool,thisboollist) #need to ensure that the list is unique!
					thisbool=string(thisbool,"_2")
				end
				push!(thisboollist,thisbool)
				if boolNumVarsUsedByModel[i][j]
					write(fiostream,"Dim ",thisbool," as Boolean: ",thisbool," = False\r\n")
				end
		end
		@assert unique(thisboollist)==thisboollist
		push!(boollist,copy(thisboollist))
	end
	write(fiostream,"\r\n")
	return boollist
end

function vba_write_writeIterations(indent::Int,fiostream::IOStream,iteration::Int,bt::BoostedTree,boolListNum::Array{Array{String,1},1},boolListChar::Array{Array{String,1},1},df_name_vector::Array{String,1},candMatWOMaxValues::Array{Array{Float64,1},1},mappings::Array{Array{String,1},1},number_of_num_features)
	#write(fiostream,"private double Iteration$(iteration)(double dRawscore)\r\n{\r\n")
	tree=bt.trees[iteration]
	#orig_id=tree.featid
	#orig_id<0 ? this_id=number_of_num_features-orig_id : this_id=orig_id
	#write(fiostream," " ^ indent)
	#if orig_id>0
	#	write(fiostream,"if (",boolListNum[this_id][tree.subset[end]],")\r\n{\r\n")
	#else
#		write(fiostream,"if (",join(boolListChar[-orig_id][[tree.subset]],"||"),")\r\n{\r\n")
#	end
			vba_write_writeIterations_recursive(bt.moderationvector[iteration],indent,fiostream,tree,boolListNum,boolListChar,df_name_vector,candMatWOMaxValues,mappings,number_of_num_features)
#		write(fiostream," " ^ (indent-1))
#		write(fiostream," }\r\nelse\r\n{\r\n")
#			vba_write_writeIterations_recursive(indent+1,fiostream,tree.right,boolListNum,boolListChar,df_name_vector,candMatWOMaxValues,mappings)
#		write(fiostream," " ^ (indent-1))
#		write(fiostream,"}\r\n")

	#write(fiostream,"return dRawscore;\r\n}\r\n\r\n")
	return nothing
end

function vba_write_writeIterations_recursive(mdf::Float64,indent::Int,fiostream::IOStream,tree::Node{T},boolListNum::Array{Array{String,1},1},boolListChar::Array{Array{String,1},1},df_name_vector::Array{String,1},candMatWOMaxValues::Array{Array{Float64,1},1},mappings::Array{Array{String,1},1},number_of_num_features::Int) where T<:Unsigned 
	orig_id=tree.featid
	orig_id<0 ? this_id=number_of_num_features-orig_id : this_id=orig_id
	write(fiostream," " ^ indent)
	if orig_id>0
		write(fiostream,repeat("\t",indent),"If (",boolListNum[this_id][tree.subset[end]],") Then\r\n")
	else
		write(fiostream,repeat("\t",indent),"If (",join(boolListChar[-orig_id][collect(tree.subset)]," Or "),") Then\r\n")
	end
		#typeof(tree.left)!=Leaf ? write(fiostream,"\r\n{") : write(fiostream,'{')
		#write(fiostream,"EndIf\r\n")
			vba_write_writeIterations_recursive(mdf,indent+1,fiostream,tree.left,boolListNum,boolListChar,df_name_vector,candMatWOMaxValues,mappings,number_of_num_features)
		#if typeof(tree.left)!=Leaf
		#	write(fiostream," " ^ (indent-1))
			write(fiostream,repeat("\t",indent),"Else\r\n")
		#end
		#typeof(tree.left)!=Leaf ? write(fiostream,"\r\n{") : write(fiostream,'{')
		#write(fiostream,"\r\n")
				vba_write_writeIterations_recursive(mdf,indent+1,fiostream,tree.right,boolListNum,boolListChar,df_name_vector,candMatWOMaxValues,mappings,number_of_num_features)
		#if typeof(tree.right)!=Leaf
#		write(fiostream,"\r\n{") : write(fiostream,'{')
		#	write(fiostream," " ^ (indent-1))
			write(fiostream,repeat("\t",indent),"EndIf\r\n")
		#end
	return nothing
end

function vba_write_writeIterations_recursive(mdf::Float64,indent::Int,fiostream::IOStream,tree::Leaf,boolListNum::Array{Array{String,1},1},boolListChar::Array{Array{String,1},1},df_name_vector::Array{String,1},candMatWOMaxValues::Array{Array{Float64,1},1},mappings::Array{Array{String,1},1},number_of_num_features::Int)
	val=_moderate(tree.fitted,mdf)	
	write(fiostream,repeat("\t",indent),"dRawscore = dRawscore * $(val)\r\n") # }\r\n") #//Leafnumber=$(tree.id)")
	return nothing
end

function write_and_create_boollist_num(fiostream::IOStream,df_name_vector::Array{String,1},number_of_num_features::Int,candMatWOMaxValues::Array{Array{Float64,1},1},boolNumVarsUsedByModel::Array{Array{Bool,1},1})
#NOTE, here we assume that the numerical variables are "first" in the df_name_vector
#todo: only define booleans which are actually used by the model! this will keep the code smaller and likely more efficient
boollist=Array{Array{String,1}}(0)
	for i=1:number_of_num_features
		thisvarname=df_name_vector[i]
		thisboollist=Array{String}(0)
		for j=1:length(candMatWOMaxValues[i])
				thisbool=convert(String,string(thisvarname,"_isLE_",numberAsStringForCSharp(candMatWOMaxValues[i][j])))
				while in(thisbool,thisboollist) #need to ensure that the list is unique!
					thisbool=string(thisbool,"_2")
				end
				push!(thisboollist,thisbool)
				if boolNumVarsUsedByModel[i][j]
					write(fiostream,"private bool ",thisbool,"=false;\r\n")
				end
		end
		@assert unique(thisboollist)==thisboollist
		push!(boollist,copy(thisboollist))
	end
	write(fiostream,"\r\n")
	return boollist
end

function create_custom_string(instr::String)
	x=strip(instr) #leading/trailing blanks
	#replace(x,' ','_') #todo, find out why this is not done by the loop below
	for char in x #i=1:length(x)
		if !in(char,CSHARP_VALID_CHARS)
		#if !in(x[i],CSHARP_VALID_CHARS)
			x=replace(x,char,"_")
			#x=replace(x,x[i],"_")
		end
	end
	return x
end

function write_and_create_boollist_char(fiostream::IOStream,mappings::Array{Array{String,1}},df_name_vector::Array{String,1},number_of_num_features::Int,boolCharVarsUsedByModel::Array{Array{Bool,1},1})
#note todo tbd check this: here we only consider the values which we observe in the training data
#of course we would also know the values from the validation data, but that should not matter here, since any new data which the module
#will see, might have even other/new values in it.
#Does that make sense? Or should we include the values from the val data in this list?
	boollist=Array{Array{String,1}}(0)
	for i=1:length(mappings)
		thisvarname=df_name_vector[number_of_num_features+i]
		thisboollist=Array{String}(0)
		for j=1:length(mappings[i])
				str=create_custom_string(mappings[i][j]) #this can be any String, let us see what C# can handle
				#if C# isunable to handle the UTF8Strings, we need to convert it to some ASCIIString, while ensuring that we still have UNIQUENESS, todo /tbd check this
				thisbool=convert(String,string(thisvarname,"_",str))
				while in(thisbool,thisboollist) #need to ensure that the list is unique!
					thisbool=convert(String,string(thisbool,"_2"))
				end
				push!(thisboollist,thisbool)
				if boolCharVarsUsedByModel[i][j]
					write(fiostream,"private bool ",thisbool,"=false;\r\n")
				end
		end
		@assert unique(thisboollist)==thisboollist
		push!(boollist,copy(thisboollist))
	end
	write(fiostream,"\r\n")
	return boollist
end

function csharp_write_RequiredElementsProvided(fiostream::IOStream,df_name_vector::Array{String,1},boolVariablesUsed::Array{Bool,1})
	start=
	"""
	public bool RequiredElementsProvided()
		{
		// AllRequiredElementsProvided=true;
	"""
	write(fiostream,start)
	write(fiostream,"AllRequiredElementsProvided=")
	for i=1:length(df_name_vector)
		if boolVariablesUsed[i]
			str=df_name_vector[i]
			if i>1
				write(fiostream," && ")
			end
			write(fiostream,str,"_Provided")
		end
	end
	write(fiostream,";\r\n")
	for i=1:length(df_name_vector)
		if boolVariablesUsed[i]
			str=df_name_vector[i]
			write(fiostream,"if (",str,"_Provided==false) { ScoreErrors.Add(new ScoreError(ScoreError.ResponseCode.RequiredValueNotProvided, ",DoubleQuote,str,DoubleQuote,", string.Empty)); }\r\n")
		end
	end
	endstr=
	"""
	return AllRequiredElementsProvided;
	}

	"""
	write(fiostream,endstr)
	return nothing
end


function vba_write_evaluate_profile_fn()

	return nothing
end

function vba_write_booleans(fiostream::IOStream,boolListNum::Array{Array{String,1},1},boolListChar::Array{Array{String,1},1},df_name_vector::Array{String,1},candMatWOMaxValues::Array{Array{Float64,1},1},mappings::Array{Array{String,1},1},boolCharVarsUsedByModel::Array{Array{Bool,1},1},boolNumVarsUsedByModel::Array{Array{Bool,1},1},boolVariablesUsed::Array{Bool,1})
	start=
	"""
	'Definition of Booleans
	"""
	write(fiostream,start)
	#numeric vars
		for i=1:length(boolListNum)
			if boolVariablesUsed[i]
				#str=lowercase(df_name_vector[i])
				str=df_name_vector[i]
				#write(fiostream,"if (!Double.TryParse(strValue, NumberStyles.Number, numberFormatCulture, out dResult)) return ScoreError.ResponseCode.ValueIsNotAValidNumber;\r\n")
				write(fiostream,"'initially all LE booleans are set to false <-> all numeric values are infinite\r\n")
				write(fiostream,"'Note that VBA displays 'x <= 2.0' as 'x <= 2.#' if x is of type Double\r\n")
				for j=1:length(boolListNum[i])
					if boolNumVarsUsedByModel[i][j]
						write(fiostream,"if (",str," <= $(candMatWOMaxValues[i][j])) Then \r\n\t",boolListNum[i][j]," = True \r\n EndIf \r\n")						
					end
				end
				
			end
		end
	#character vars
		for i=1:length(boolListChar)
			i_index=length(boolListNum)+i
			if boolVariablesUsed[i_index]
				#str=lowercase(df_name_vector[i])
				str=df_name_vector[i_index]
				#write(fiostream,"case ",DoubleQuote,lowercase(str),DoubleQuote,':',"\r\n",str,"_Provided=true;\r\n")
				#write(fiostream,"if (strValue.Length==0) return ScoreError.ResponseCode.BlankValuesNotValid;\r\n")
				write(fiostream,"Select Case ",str,"\r\n")
				for j=1:length(boolListChar[i])
					if boolCharVarsUsedByModel[i][j]
						write(fiostream,"Case Ucase(",DoubleQuote,string(mappings[i][j]),DoubleQuote,")\r\n\t$(boolListChar[i][j]) = True \r\n")
					end
				end
				write(fiostream,"End Select\r\n\r\n")
			end
		end
	endstr=
	""" 
	'End of the definition of the Boolean variables

	'Start of the definition of the trees
	"""
	write(fiostream,endstr)
	return nothing
end

function csharp_write_ScoreErrorResponseCodeCheck(fiostream::IOStream,boolListNum::Array{Array{String,1},1},boolListChar::Array{Array{String,1},1},df_name_vector::Array{String,1},candMatWOMaxValues::Array{Array{Float64,1},1},mappings::Array{Array{String,1},1},boolCharVarsUsedByModel::Array{Array{Bool,1},1},boolNumVarsUsedByModel::Array{Array{Bool,1},1},boolVariablesUsed::Array{Bool,1})
	start=
	"""
	public ScoreError.ResponseCode Check(string strKey, string strValue)
		{
		switch(strKey.ToLower())
			{
	"""
	write(fiostream,start)
	#numeric vars
		for i=1:length(boolListNum)
			if boolVariablesUsed[i]
				#str=lowercase(df_name_vector[i])
				str=df_name_vector[i]
				write(fiostream,"case ",DoubleQuote,lowercase(str),DoubleQuote,':',"\r\n",str,"_Provided=true;\r\n",'{',"\r\n")
				write(fiostream,"if (strValue.Length==0) return ScoreError.ResponseCode.BlankValuesNotValid;\r\n double dResult = 0.0d;\r\n")
				write(fiostream,"if (!Double.TryParse(strValue, NumberStyles.Number, numberFormatCulture, out dResult)) return ScoreError.ResponseCode.ValueIsNotAValidNumber;\r\n")
				write(fiostream,"//initially all LE booleans are set to false <-> all numeric values are infinite\r\n")
				for j=1:length(boolListNum[i])
					if boolNumVarsUsedByModel[i][j]
						write(fiostream,"if (dResult<=$(candMatWOMaxValues[i][j])d) {",boolListNum[i][j],"=true;}\r\n")
					end
				end
				#write(fiostream,"// todo/tbd: implement a check: test whether the value is within min and maximum observed values\r\n")
				write(fiostream,"return ScoreError.ResponseCode.Passed;\r\n}\r\n")
				#=
					if (lResult >= minobserved && lResult <= maxobserved ) return ScoreError.ResponseCode.Passed;
								if (lResult < minobserved) return ScoreError.ResponseCode.ValueBelowModelBounds;
								if (lResult > 14 // maxobserved ) return ScoreError.ResponseCode.ValueAboveModelBounds;
								return ScoreError.ResponseCode.Passed;
							}
				=#
			end
		end
	#character vars
		for i=1:length(boolListChar)
			i_index=length(boolListNum)+i
			if boolVariablesUsed[i_index]
				#str=lowercase(df_name_vector[i])
				str=df_name_vector[i_index]
				write(fiostream,"case ",DoubleQuote,lowercase(str),DoubleQuote,':',"\r\n",str,"_Provided=true;\r\n")
				write(fiostream,"if (strValue.Length==0) return ScoreError.ResponseCode.BlankValuesNotValid;\r\n")
				write(fiostream,"switch(strValue)\r\n{\r\n")
				for j=1:length(boolListChar[i])
					if boolCharVarsUsedByModel[i][j]
						write(fiostream,"case ",DoubleQuote,string(mappings[i][j]),DoubleQuote,": $(boolListChar[i][j]) = true; return ScoreError.ResponseCode.Passed;\r\n")
					end
				end
				write(fiostream,"default: return ScoreError.ResponseCode.CategoricalValueNotAllowed;\r\n}\r\n")
			end
		end
	endstr=
	"""
			}
			return ScoreError.ResponseCode.Passed;
		}

	"""
	write(fiostream,endstr)
	return nothing
end

function csharp_write_elseif_ScoreMap(fiostream::IOStream,intArray_1_to_n::Array{Int,1},maxRawRelativityPerScoreSorted::Array{Float64,1},estimatesPerScore::Array{Float64,1})
	#NOTE: here we assume that the relativities are sorted in increasing order todo/tbd take this on a list of core assumptions
	minval=10
	if length(maxRawRelativityPerScoreSorted)>minval
		middle=fld(length(maxRawRelativityPerScoreSorted),2)
		middle=max(2,min(length(maxRawRelativityPerScoreSorted)-1,middle))
		write(fiostream,"if (a_dRawScore <= $(maxRawRelativityPerScoreSorted[middle])d)\r\n{")
			csharp_write_elseif_ScoreMap(fiostream,intArray_1_to_n[1:middle],maxRawRelativityPerScoreSorted[1:middle],estimatesPerScore[1:middle])
		write(fiostream,"}\r\nelse{\r\n")
			csharp_write_elseif_ScoreMap(fiostream,intArray_1_to_n[middle+1:end],maxRawRelativityPerScoreSorted[middle+1:end],estimatesPerScore[middle+1:end])
		write(fiostream,"}\r\n")
	else
		write(fiostream,"if (a_dRawScore <= $(maxRawRelativityPerScoreSorted[1])d) { iScore=$(intArray_1_to_n[1]); dEstimate= $(estimatesPerScore[1])d; }\r\n")
		for j=2:length(maxRawRelativityPerScoreSorted)
			write(fiostream,"else if (a_dRawScore <= $(maxRawRelativityPerScoreSorted[j])d) { iScore=$(intArray_1_to_n[j]); dEstimate= $(estimatesPerScore[j])d; }\r\n")
		end
		#Last else condition is "special"
	end
	return nothing
end


function sas_write_ScoreMap(fiostream::IOStream,maxRawRelativityPerScoreSorted::Array{Float64,1},estimatesPerScore::Array{Float64,1})
	start=
	"""

	/*map raw relativity to score and fitted value*/
		raw_relativity = 1.0;
		score = -1;
		fitted_value = 0.0;		
		
		array tmpnum{*} rel_mod_leaf_iter:; 
		do i_loop=1 to dim(tmpnum); 
			raw_relativity = raw_relativity * tmpnum[i_loop]; 
		end;
		
	"""
	write(fiostream,start)
	sas_write_elseif_ScoreMap(fiostream,"",0,collect(1:length(maxRawRelativityPerScoreSorted)),maxRawRelativityPerScoreSorted,estimatesPerScore)
	endstr=
	"""
	
	run;
	"""
	write(fiostream,endstr)
	return nothing
end

function sas_write_elseif_ScoreMap(fiostream::IOStream,prefix::String,indent::Int,intArray_1_to_n::Array{Int,1},maxRawRelativityPerScoreSorted::Array{Float64,1},estimatesPerScore::Array{Float64,1})
    #NOTE: here we assume that the relativities are sorted in increasing order todo/tbd take this on a list of core assumptions
    minval=10    
    if length(maxRawRelativityPerScoreSorted)>minval
        middle=fld(length(maxRawRelativityPerScoreSorted),2)
		middle=max(2,min(length(maxRawRelativityPerScoreSorted)-1,middle))
		valtmp=maxRawRelativityPerScoreSorted[middle]
        write(fiostream,repeat("\t",indent),prefix,"If (raw_relativity <= $(valtmp)) then do;\r\n")
            sas_write_elseif_ScoreMap(fiostream,"",indent+1,intArray_1_to_n[1:middle],maxRawRelativityPerScoreSorted[1:middle],estimatesPerScore[1:middle])        
		write(fiostream,repeat("\t",max(0,indent)),"eND;/* end of if clause - value <= $(valtmp)*/\r\n")		
		if length(intArray_1_to_n[middle+1:end])>0
			write(fiostream,repeat("\t",max(0,indent)),"ELSE DO;/* value > $(valtmp)*/\r\n")
			valtmp2=maxRawRelativityPerScoreSorted[middle+1:end][1]
        	sas_write_elseif_ScoreMap(fiostream,"",indent+1,intArray_1_to_n[middle+1:end],maxRawRelativityPerScoreSorted[middle+1:end],estimatesPerScore[middle+1:end])
			write(fiostream,repeat("\t",max(0,indent)),"End;/* end of else do clause - value > $(valtmp)*/\r\n")
		end
    else        
        write(fiostream,repeat("\t",max(0,indent)),prefix,"iF (raw_relativity <= $(maxRawRelativityPerScoreSorted[1])) then do;score = $(intArray_1_to_n[1]); fitted_value = $(estimatesPerScore[1]);END;\r\n")
        for j=2:length(maxRawRelativityPerScoreSorted)
            write(fiostream,repeat("\t",indent),"elSE if (raw_relativity <= $(maxRawRelativityPerScoreSorted[j])) then do;score = $(intArray_1_to_n[j]); fitted_value = $(estimatesPerScore[j]);EnD;\r\n")            
        end        
    end
    return nothing
end


function vba_write_ScoreMap(fiostream::IOStream,maxRawRelativityPerScoreSorted::Array{Float64,1},estimatesPerScore::Array{Float64,1})
	start=
	"""

	Function derive_score_from_raw_relativity(dRawScore as Double) as Integer 'as Collection
		'Dim cResult as New Collection
		Dim iScore as Integer
		Dim dEstimate as Double
		iScore = -1
		dEstimate = 0.0

		If dRawScore <= 0 Then
			'This should not happen
			derive_score_from_raw_relativity = -1
			Exit Function
		End If
	"""
	write(fiostream,start)
	vba_write_elseif_ScoreMap(fiostream,0,collect(1:length(maxRawRelativityPerScoreSorted)),maxRawRelativityPerScoreSorted,estimatesPerScore)
	endstr=
	"""

		'cResult.add CStr(iScore), "Score"
		'cResult.add CStr(dEstimate), "Estimate"

		derive_score_from_raw_relativity = iScore

	End Function

	"""
	write(fiostream,endstr)

	start2=
	"""

	Function derive_estimate_from_raw_relativity(dRawScore as Double) as Double		
		Dim iScore as Integer
		Dim dEstimate as Double
		iScore = -1
		dEstimate = 0.0

		If dRawScore <= 0 Then
        	'This should not happen
        	derive_estimate_from_raw_relativity = -999.999
        	Exit Function
    	End If    
	"""
	write(fiostream,start2)	
	vba_write_elseif_ScoreMap(fiostream,0,collect(1:length(maxRawRelativityPerScoreSorted)),maxRawRelativityPerScoreSorted,estimatesPerScore)
	endstr2=
	"""

		derive_estimate_from_raw_relativity = dEstimate

	End Function

	"""
	write(fiostream,endstr2)
	return nothing
end

function vba_write_elseif_ScoreMap(fiostream::IOStream,indent::Int,intArray_1_to_n::Array{Int,1},maxRawRelativityPerScoreSorted::Array{Float64,1},estimatesPerScore::Array{Float64,1})
	#NOTE: here we assume that the relativities are sorted in increasing order todo/tbd take this on a list of core assumptions
	minval=10	
	if length(maxRawRelativityPerScoreSorted)>minval
		middle=fld(length(maxRawRelativityPerScoreSorted),2)
		middle=max(2,min(length(maxRawRelativityPerScoreSorted)-1,middle))
		write(fiostream,repeat("\t",indent),"If (dRawScore <= $(maxRawRelativityPerScoreSorted[middle])) Then\r\n")
			vba_write_elseif_ScoreMap(fiostream,indent+1,intArray_1_to_n[1:middle],maxRawRelativityPerScoreSorted[1:middle],estimatesPerScore[1:middle])
		write(fiostream,repeat("\t",indent),"Else\r\n")
			vba_write_elseif_ScoreMap(fiostream,indent+1,intArray_1_to_n[middle+1:end],maxRawRelativityPerScoreSorted[middle+1:end],estimatesPerScore[middle+1:end])
		write(fiostream,repeat("\t",indent),"EndIf\r\n")
	else
		indent_score=repeat("\t",indent+1)
		write(fiostream,repeat("\t",indent),"If (dRawScore <= $(maxRawRelativityPerScoreSorted[1])) Then\r\n",indent_score,"iScore = $(intArray_1_to_n[1])\r\n",indent_score,"dEstimate = $(estimatesPerScore[1])\r\n")
		for j=2:length(maxRawRelativityPerScoreSorted)
			write(fiostream,repeat("\t",indent),"ElseIf (dRawScore <= $(maxRawRelativityPerScoreSorted[j])) Then\r\n",indent_score,"iScore=$(intArray_1_to_n[j])\r\n",indent_score,"dEstimate = $(estimatesPerScore[j])\r\n")
			if j==length(maxRawRelativityPerScoreSorted)
				write(fiostream,repeat("\t",indent),"Endif\r\n")
			end
		end
		#Last else condition is "special"
	end
	return nothing
end


function csharp_write_ScoreMap(fiostream::IOStream,maxRawRelativityPerScoreSorted::Array{Float64,1},estimatesPerScore::Array{Float64,1})
	start=
	"""
	private int ScoreMap(double a_dRawScore)
		{
			int iScore = -1;
			double dEstimate = 0.0d;
	"""
	write(fiostream,start)
	csharp_write_elseif_ScoreMap(fiostream,collect(1:length(maxRawRelativityPerScoreSorted)),maxRawRelativityPerScoreSorted,estimatesPerScore)
	endstr=
	"""
	FinalEstimate = dEstimate;
			return iScore;
		}

	"""
	write(fiostream,endstr)
	return nothing
end

function csharp_write_GenerateScore(fiostream::IOStream,bt::BoostedTree)
start=
"""
public int GenerateScore()
	{
		double dRawscore = 1.0d;
		TotalDefaultedIterations = 0;
"""
write(fiostream,start)
if global_debug_csharp
	write(fiostream,"double[] relarr = new double[$(length(bt.trees))];\r\n");
end
for i=1:length(bt.trees)
	write(fiostream,"dRawscore = Iteration$(i)(dRawscore);\r\n")
	if global_debug_csharp
	write(fiostream,"relarr[$(i-1)]=dRawscore;\r\n")
end
end

endstr=
	"""
	if (TotalDefaultedIterations > 0) { FinalScore = -1; return FinalScore; }
		RawScore = dRawscore;
		FinalScore = ScoreMap(dRawscore);
		return FinalScore;
	}

	"""
	write(fiostream,endstr)
return nothing
end


function csharp_write_writeIterations(indent::Int,fiostream::IOStream,iteration::Int,bt::BoostedTree,boolListNum::Array{Array{String,1},1},boolListChar::Array{Array{String,1},1},df_name_vector::Array{String,1},candMatWOMaxValues::Array{Array{Float64,1},1},mappings::Array{Array{String,1},1},number_of_num_features)
	write(fiostream,"private double Iteration$(iteration)(double dRawscore)\r\n{\r\n")
	tree=bt.trees[iteration]
	#orig_id=tree.featid
	#orig_id<0 ? this_id=number_of_num_features-orig_id : this_id=orig_id
	#write(fiostream," " ^ indent)
	#if orig_id>0
	#	write(fiostream,"if (",boolListNum[this_id][tree.subset[end]],")\r\n{\r\n")
	#else
#		write(fiostream,"if (",join(boolListChar[-orig_id][[tree.subset]],"||"),")\r\n{\r\n")
#	end
			csharp_write_writeIterations_recursive(bt.moderationvector[iteration],indent+1,fiostream,tree,boolListNum,boolListChar,df_name_vector,candMatWOMaxValues,mappings,number_of_num_features)
#		write(fiostream," " ^ (indent-1))
#		write(fiostream," }\r\nelse\r\n{\r\n")
#			csharp_write_writeIterations_recursive(indent+1,fiostream,tree.right,boolListNum,boolListChar,df_name_vector,candMatWOMaxValues,mappings)
#		write(fiostream," " ^ (indent-1))
#		write(fiostream,"}\r\n")

	write(fiostream,"return dRawscore;\r\n}\r\n\r\n")
	return nothing
end

function csharp_write_writeIterations_recursive(mdf::Float64,indent::Int,fiostream::IOStream,tree::Node{T},boolListNum::Array{Array{String,1},1},boolListChar::Array{Array{String,1},1},df_name_vector::Array{String,1},candMatWOMaxValues::Array{Array{Float64,1},1},mappings::Array{Array{String,1},1},number_of_num_features::Int) where T<:Unsigned 
	orig_id=tree.featid
	orig_id<0 ? this_id=number_of_num_features-orig_id : this_id=orig_id
	write(fiostream," " ^ indent)
	if orig_id>0
		write(fiostream,"\r\nif (",boolListNum[this_id][tree.subset[end]],")")
	else
		write(fiostream,"\r\nif (",join(boolListChar[-orig_id][collect(tree.subset)],"||"),")")
	end
		#typeof(tree.left)!=Leaf ? write(fiostream,"\r\n{") : write(fiostream,'{')
		write(fiostream,'{')
			csharp_write_writeIterations_recursive(mdf,indent+1,fiostream,tree.left,boolListNum,boolListChar,df_name_vector,candMatWOMaxValues,mappings,number_of_num_features)
		#if typeof(tree.left)!=Leaf
		#	write(fiostream," " ^ (indent-1))
			write(fiostream," }\r\nelse")
		#end
		#typeof(tree.left)!=Leaf ? write(fiostream,"\r\n{") : write(fiostream,'{')
		write(fiostream,'{')
				csharp_write_writeIterations_recursive(mdf,indent+1,fiostream,tree.right,boolListNum,boolListChar,df_name_vector,candMatWOMaxValues,mappings,number_of_num_features)
		#if typeof(tree.right)!=Leaf
#		write(fiostream,"\r\n{") : write(fiostream,'{')
		#	write(fiostream," " ^ (indent-1))
			write(fiostream,"}\r\n")
		#end
	return nothing
end

function csharp_write_writeIterations_recursive(mdf::Float64,indent::Int,fiostream::IOStream,tree::Leaf,boolListNum::Array{Array{String,1},1},boolListChar::Array{Array{String,1},1},df_name_vector::Array{String,1},candMatWOMaxValues::Array{Array{Float64,1},1},mappings::Array{Array{String,1},1},number_of_num_features::Int)
	val=_moderate(tree.fitted,mdf)
	write(fiostream," dRawscore*=$(val)d;") # }\r\n") #//Leafnumber=$(tree.id)")
	return nothing
end

function determine_used_variables(bt::BoostedTree)
	boolNumVarsUsedByModel=Array{Array{Bool,1}}(0)
	boolCharVarsUsedByModel=Array{Array{Bool,1}}(0)
	#settings	#number_of_char_features
	number_of_num_features=bt.settings.number_of_num_features
	number_of_char_features=bt.settings.number_of_char_features	
	intVarsUsed=bt.intVarsUsed
	for n=1:number_of_num_features 
		push!(boolNumVarsUsedByModel,zeros(Bool,length(intVarsUsed[1][1:number_of_num_features][n])))
	end
	for n=1:number_of_char_features
		push!(boolCharVarsUsedByModel,zeros(Bool,length(intVarsUsed[1][1+number_of_num_features:end][n])))
	end

	for iter=1:length(intVarsUsed)
		for n=1:number_of_num_features
			x=intVarsUsed[iter][1:number_of_num_features]
			boolNumVarsUsedByModel[n][:]=(x[n].>0).|boolNumVarsUsedByModel[n]
		end
		for n=1:number_of_char_features
			x=intVarsUsed[iter][1+number_of_num_features:end]
			boolCharVarsUsedByModel[n][:]=(x[n].>0).|boolCharVarsUsedByModel[n]
		end
	end

	boolVariablesUsed=zeros(Bool,length(intVarsUsed[1])) 
	number_of_nums=number_of_num_features
	for i=1:length(boolVariablesUsed)
		if i>number_of_nums
			j=i-number_of_nums
			if sum(boolCharVarsUsedByModel[j])>0
				boolVariablesUsed[i]=true
			end
		else
			if sum(boolNumVarsUsedByModel[i])>0
				boolVariablesUsed[i]=true
			end
		end
	end
return boolVariablesUsed,boolNumVarsUsedByModel,boolCharVarsUsedByModel
end

function write_vba_code(vectorOfLeafArrays::Array{Array{Leaf,1},1},estimatesPerScore::Array{Float64,1},candMatWOMaxValues::Array{Array{Float64,1},1},bt::BoostedTree,fileloc::String,settings::String,mappings::Array{Array{String,1},1},indent::Int,sett::ModelSettings)
	@assert size(estimatesPerScore)==size(bt.maxRawRelativityPerScoreSorted)
	number_of_num_features=sett.number_of_num_features
	df_name_vector=sett.df_name_vector
	iterations=size(bt.trees,1)

	fiostream=open(fileloc,"a")
	start_str="""
Function Score() As Integer
    Application.Volatile
    Dim rel As Double
    Dim sc As Integer
    rel = readDataAndDeriveRawRelativity()
    sc = -1
    sc = derive_score_from_raw_relativity(rel)
    Score = sc
End Function
	
	"""
	write(fiostream,start_str)
	sigstr,strEvaluteProfile=vba_get_signature(mappings,df_name_vector,number_of_num_features)
	write(fiostream,strEvaluteProfile)
	write(fiostream,"Function deriveRawRelativity")	
	write(fiostream,'(')	
	write(fiostream,sigstr)
	write(fiostream,")\r\n")	
	write(fiostream,"\r\n'Boosted Tree produced by Julia. Time: $(now())\r\n")
	#generate_csharpheader(fiostream)	
	write(fiostream,"\r\n")
	boolVariablesUsed,boolNumVarsUsedByModel,boolCharVarsUsedByModel=determine_used_variables(bt)
	
	#NOTE: these two lists are not really booleans, but merely strings which define booleans in the C# code, the notation may be slighly misleading.
	boolListNum=write_and_create_boollist_num_vba(fiostream,df_name_vector,number_of_num_features,candMatWOMaxValues,boolNumVarsUsedByModel)
	@assert length(boolListNum)==length(candMatWOMaxValues) #todo/tbd we could even check here whether each bool list has the right size (then again it should be fulfilled
	boolListChar=write_and_create_boollist_char_vba(fiostream,mappings,df_name_vector,number_of_num_features,boolCharVarsUsedByModel)
	@assert length(boolListChar)==length(mappings)
	#error("bk continue here")
	#csharp_write_RequiredElementsProvided(fiostream,df_name_vector,boolVariablesUsed)
	warn("todo, what if some variables are not as expected (unseen values for instance). Improve this / error handling similar to the generated SAS Code!")
	vba_write_booleans(fiostream,boolListNum,boolListChar,df_name_vector,candMatWOMaxValues,mappings,boolCharVarsUsedByModel,boolNumVarsUsedByModel,boolVariablesUsed)
	
	#for i=1:length(df_name_vector)
	#	if boolVariablesUsed[i]
	#		write(fiostream,"private bool ",df_name_vector[i],"_Provided=false;\r\n")
	#	end
	#end
	#write(fiostream,"\r\n")
	
	#todo tbd, take note of this somewhere and do not forget it!
	#for the c# code the last entry of maxRawRelativityPerScoreSorted is modified to be 50% (50% is arbitrary here) higher than the value which was derived by Julia
	#this means that any policy will get the maximal score if it is "worse" than any training policy
	#but if it is much worse, it will get the default score of -1
	modified_maxRawRelativityPerScoreSorted=deepcopy(bt.maxRawRelativityPerScoreSorted)
	modified_maxRawRelativityPerScoreSorted[end]=modified_maxRawRelativityPerScoreSorted[end]*1.5
	#csharp_write_GenerateScore(fiostream,bt)
	write(fiostream,"Dim dRawscore as Double\r\n 'Initialize Raw Score \r\n dRawscore=1.0\r\n\r\n")
	for i=1:length(bt.trees)
		write(fiostream,"\r\n'Iteration $(i)\r\n")
		vba_write_writeIterations(0,fiostream,i,bt,boolListNum,boolListChar,df_name_vector,candMatWOMaxValues,mappings,number_of_num_features)
	end
	#write end of the code
endstr="""
deriveRawRelativity = dRawscore 
End Function
"""
	write(fiostream,endstr)	
        vba_write_ScoreMap(fiostream,modified_maxRawRelativityPerScoreSorted,estimatesPerScore)	
	close(fiostream)

	return nothing
end	


function write_csharp_code(vectorOfLeafArrays::Array{Array{Leaf,1},1},estimatesPerScore::Array{Float64,1},candMatWOMaxValues::Array{Array{Float64,1},1},bt::BoostedTree,fileloc::String,settings::String,mappings::Array{Array{String,1},1},indent::Int,sett::ModelSettings)
	@assert size(estimatesPerScore)==size(bt.maxRawRelativityPerScoreSorted)
	number_of_num_features=sett.number_of_num_features
	df_name_vector=sett.df_name_vector
	iterations=size(bt.trees,1)
	fiostream=open(fileloc,"a")
	#settings2=copy(settings)
	#settings2=replace(settings2,"\n","\r\n")
	#settings2=replace(settings2,"\r\r\n","\r\n")
	#write(fiostream,"/* \r\nSettings: \r\n",settings2,"\r\n*/ \r\n")
	generate_csharpheader(fiostream)
	csharp_write_pub_str(fiostream,"AnalysisObjective","TBD_TODO_CompetitorPremium")
	csharp_write_pub_str(fiostream,"AnalysisName","TBD_TODO_SomeModelName")
	csharp_write_pub_str(fiostream,"AnalysisId","TBD_TODO_a_unique_model_id as a string (e.g. a hash?)")
	csharp_write_pub_str(fiostream,"AnalysisDate","TBD_TODO_October 21, 1984, 21:18:10")
	write(fiostream,"\r\n")

	boolVariablesUsed,boolNumVarsUsedByModel,boolCharVarsUsedByModel=determine_used_variables(bt)
	for i=1:length(df_name_vector)
		if boolVariablesUsed[i]
			write(fiostream,"private bool ",df_name_vector[i],"_Provided=false;\r\n")
		end
	end
	write(fiostream,"\r\n")
	#NOTE: these two lists are not really booleans, but merely strings which define booleans in the C# code, the notation may be slightly misleading.
	boolListNum=write_and_create_boollist_num(fiostream,df_name_vector,number_of_num_features,candMatWOMaxValues,boolNumVarsUsedByModel)
	@assert length(boolListNum)==length(candMatWOMaxValues) #todo/tbd we could even check here whether each bool list has the right size (then again it should be fulfilled
	boolListChar=write_and_create_boollist_char(fiostream,mappings,df_name_vector,number_of_num_features,boolCharVarsUsedByModel)
	@assert length(boolListChar)==length(mappings)
	csharp_write_RequiredElementsProvided(fiostream,df_name_vector,boolVariablesUsed)
	csharp_write_ScoreErrorResponseCodeCheck(fiostream,boolListNum,boolListChar,df_name_vector,candMatWOMaxValues,mappings,boolCharVarsUsedByModel,boolNumVarsUsedByModel,boolVariablesUsed)
	#todo tbd, take note of this somewhere and do not forget it!
	#for the c# code the last entry of maxRawRelativityPerScoreSorted is modified to be 50% (50% is arbitrary here) higher than the value which was derived by Julia
	#this means that any policy will get the maximal score if it is "worse" than any training policy
	#but if it is much worse, it will get the default score of -1
	modified_maxRawRelativityPerScoreSorted=deepcopy(bt.maxRawRelativityPerScoreSorted)
	modified_maxRawRelativityPerScoreSorted[end]=modified_maxRawRelativityPerScoreSorted[end]*1.5
	csharp_write_ScoreMap(fiostream,modified_maxRawRelativityPerScoreSorted,estimatesPerScore)
	csharp_write_GenerateScore(fiostream,bt)
	for i=1:length(bt.trees)
		csharp_write_writeIterations(0,fiostream,i,bt,boolListNum,boolListChar,df_name_vector,candMatWOMaxValues,mappings,number_of_num_features)
	end
	#write end of the code
	write(fiostream,"}\r\n}\r\n")
	close(fiostream)
end

function generate_csharpheader(fiostream)
	x=
	"""
	using System;
	using System.Collections.Generic;
	using System.Globalization;
	using System.Linq;
	using System.Text;

	namespace MillimanScoringModule
	{

		public class ScoreError
		{
			private ResponseCode rcCode;
			private string varname;
			private string value;

			public class ResponseCode
			{
				public static ResponseCode Passed;
				public static ResponseCode RequiredValueNotProvided;
				public static ResponseCode BlankValuesNotValid;
				public static ResponseCode ValueIsNotAValidNumber;
				public static ResponseCode CategoricalValueNotAllowed;
			}

			public ScoreError(ResponseCode rcCode, string p1, string p2)
			{
				this.rcCode = rcCode;
				this.varname = p1;
				this.value = p2;
			}
		}

		public class ScoreModule
		{
			public List<int> DefaultedIterations { get; set; }

			public int FinalScore { get; private set; }
			public double FinalEstimate { get; private set; }
			public double RawScore { get; private set; }
			public int TotalDefaultedIterations { get; set; }
			public bool AllRequiredElementsProvided { get; set; }
			public List<ScoreError> ScoreErrors { get; set; }

		   public int Score(Dictionary<string, string> DataElements)
			{
				int iScore = -1;
				ScoreErrors = new List<ScoreError>();
				ScoreError.ResponseCode rcCode;
				foreach (KeyValuePair<string, string> kvp in DataElements)
				{
					if ((rcCode = Check(kvp.Key, kvp.Value)) != ScoreError.ResponseCode.Passed)
					{
						ScoreErrors.Add(new ScoreError(rcCode, kvp.Key, kvp.Value));
					}
				}
				if (ScoreErrors.Count == 0 && RequiredElementsProvided() == true)
				{
					iScore=GenerateScore();
				}
				else
				{
					FinalEstimate = -1.0d;
					FinalScore = iScore;
				}
				return iScore;
			}

			private IFormatProvider numberFormatCulture = CultureInfo.CreateSpecificCulture("de-CH");

	"""
	write(fiostream,x)
	return nothing
end

promote_rule(::Type{Node}, ::Type{Leaf}) = Node
promote_rule(::Type{Leaf}, ::Type{Node}) = Node
promote_rule(::Type{Node{UInt8}}, ::Type{Leaf}) = Node{UInt8}
promote_rule(::Type{Node{UInt16}}, ::Type{Leaf}) = Node{UInt16}
promote_rule(::Type{Leaf}, ::Type{Node{UInt8}}) = Node{UInt8}
promote_rule(::Type{Leaf}, ::Type{Node{UInt16}}) = Node{UInt16}
depth(leaf::Leaf) = 0
function depth(tree::Node{T}) where T<:Unsigned  
return  1 + max(depth(tree.left), depth(tree.right))
end
depth(tree::Tree) = depth(tree.rootnode)
create_leaves_array(t::Tree)=create_leaves_array(t.rootnode)
function create_leaves_array(t::Node{T}) where T<:Unsigned 
    return vcat(create_leaves_array(t.left),create_leaves_array(t.right)) 
#warning! it is essential that this concept/order of
# number of the leaves is consistent with the apply_tree_fn for validation data where we derive the leaf number for val
end

create_leaves_array(x::Leaf)=[x]
number_of_nodes(t::Tree)=number_of_nodes(t.rootnode)
function number_of_nodes(t::Node{T}) where T<:Unsigned
    return 1+number_of_nodes(t.left)+number_of_nodes(t.right)
end
number_of_nodes(l::Leaf)=0


# ##send data to workers functions

function sendto_module(m::Module,p::Int; args...)
	for (nm, val) in args;
		@spawnat(p, eval(m, Expr(:(=), nm, val)))
	end
end
function sendto_module(m::Module,ps::Vector{Int}; args...)
	for p in ps
		sendto_module(m,p; args...)
	end
end

function sendto(p::Int; args...)
      for (nm, val) in args
          @spawnat(p, eval(Main, Expr(:(=), nm, val)))
      end
end

function sendto(ps::Vector{Int}; args...)
  for p in ps
    sendto(p; args...)
  end
end

function get_feature_pools(f::DataFrame)
	fp=Vector{Union{Vector{String},Vector{Float64}}}(0)
	@inbounds for i in 1:size(f,2)
		push!(fp,deepcopy(f[i].pool))
	end
	return fp
end