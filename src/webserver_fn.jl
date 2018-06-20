function update_folderlist!(listening_folders,folder_list_file)
# fld(mod(Int(now()),60000),1000) #time of this minute in seconds
	local farr
	try
		complete_list=Array{String}(0)
		res=readdlm(folder_list_file)
		res2=convert(Array{String,1},res[:,1])
		resize!(listening_folders,length(res2))
		listening_folders[:]=res2
		if true
			#include subfolders
			for i=1:length(listening_folders)		
				loopdir=listening_folders[i]
				if isdir(loopdir)
					push!(complete_list,loopdir)
					#isfile(file_with_folderlist)&&rm(file_with_folderlist)
					command=`cmd /C dir $(loopdir) /s /b /o:n /ad`  #command=`cmd /C dir $(loopdir) /s /b /o:n /ad >$(file_with_folderlist)` 
					try
						data=readall(command)					
						farr=split(data,"\r\n")
						if farr[end]==""
							pop!(farr)
						end
					catch eri
						@show eri
						farr=Array{String}(0)
					end
					#command=`cmd /C dir $(loopdir) /s /b /o:n /ad >$(file_with_folderlist)` 
					#run(command)					
					#f=open(file_with_folderlist)
					#linecount=countlines(f)
					#close(f)
					#if length(farr)>0
						append!(complete_list,farr)
					#end
				end
			end
		resize!(listening_folders,length(complete_list))
		listening_folders[:]=complete_list
		end
	catch er
		@show er
		warn("Unable to read and parse folder list file \n $(folder_list_file)")
		return nothing
	end
end

function test_module_on_workers()
	try
		#list of different test cases
		clist=["tree_tiny","boosting_tiny","bagging_tiny","boost_tiny_subs","boost_tiny_subs_withrepl"] #,"boosting_testmodel_small"]
		n=nprocs()
		ntests=length(clist)
		booleans=BitArray(n);fill!(booleans,true)
		i = 1
		nextidx() = (idx=i; i+=1; idx)
		#nextidx() = (idx=[i,jj]; i<n ? i+=1 : jj+=1; idx)
		@sync begin
			for p=1:n
				if p != myid() || n == 1			
					@async begin						
						while true
							idx = nextidx()							
							if idx > n
								break
							end
							for l=1:ntests
								chosenj=clist[l]
								result= remotecall_fetch(p, DecisionTrees.run_model, [string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_$(p)",chosenj)])
								thisb=(result!=String["false"]) #expected result is a list of output files
								@assert thisb
								if !thisb
									warn("Test failed: $(idx) $(p)")									
								end
								booleans[idx]=booleans[idx]&&thisb
							end							
						end			
					end
					#booleans[p] = remotecall_fetch(p,DecisionTrees.run_model,ARGS)
					#rref=remotecall_wait(p,DecisionTrees.run_model,ARGS)
				end
			end			
		end		
		
		failed=!(sum(booleans)==n)
		if failed!=0
			@show n
			@show ntests
			@show booleans
			@show sum(booleans)			
			error("something did not work")
		else
			@info "Decision Tree Module loaded on all workers"
		end
		return true
	catch eri
		@show eri
		return false
	end	
	
	error("this line should never be executed")
	return false
end

function startlistening(folder_list_file,julia_code_folder)
	julia_code_folder=add_backslashes(julia_code_folder)
	listening_folders=Array{String}(0)
	secs=2
	while true #intentional infinite loop
		try
			update_folderlist!(listening_folders,folder_list_file)
			add_backslashes!(listening_folders)	
			printover(string("Current Time: ",now()))
			for loopdir in listening_folders
				printover("Scanning $(loopdir)")
				scan_step(loopdir,julia_code_folder)
			end #loopdir loop
		catch err #an error can occur when new folders are created or renamed, thus we have a try/catch loop here
			@show err
			@info "Error. See above. Listenting will continue in $(secs)s"
			sleep(secs)
		end
	end
	return nothing
end

function isThisFileAlreadyRunningOnAProcess(hash::Int128)
	try
		command=`wmic path win32_process where "Name='julia.exe'"`
		a,b=open(command,"r")
		sleeptime=0.05
		slept=0.0		
		while (slept<3) #the process needs a bit of time to finish
			sleep(sleeptime)
			slept+=sleeptime			
			if b.exitcode==0;break;end;
		end				
		text=readall(a)
		searchfor=string("-e ",hash)		
		result=contains(text,searchfor)		
		return result
	catch err	
		@show err
	end
	return nothing
end

function alogfile_exists_result_is_already_there(loopdir,ff)	
	if size(ff,1)==0
		return false
	else
		for xx in ff
			if length(xx)>3
				if lowercase(xx)[end-3:end]==".log"
					printover("Existing log: $(loopdir)$(xx)")
					#sleep(0.1)
					return true
				end
			end
		end
	end	
	return false
end

function getdata_and_settingsfile(flist)
	settingsf="";dataf="";
	fl2=Array{eltype(flist)}(0)
	if size(flist,1)>0
		fl=[lowercase(x) for x in deepcopy(flist)]
		#*.jld are pushed first such that they are prioritized
		for x in fl
			if length(x)>4&&(fileending(x)==".jld2")
				push!(fl2,deepcopy(x))
			end
		end		
		for x in fl
			if length(x)>4&&(fileending(x)==".csv")
				push!(fl2,deepcopy(x))
			end
		end		
	end
	if size(fl2,1)>0
		for candidate1 in fl2
			settcand=string(deepcopy(candidate1)[1:end-length(fileending(candidate1))],".settings.csv")
			if in(settcand,fl2)
				return settcand,candidate1
			end
		end
	end	
	return settingsf,dataf
end

function scan_step(loopdir,julia_code_folder)
	#this is done for each directory loopdir
	settingsf,dataf,logFilename=return_settf_dataf_logf(loopdir)
	return run_these_files(settingsf,dataf,logFilename,julia_code_folder)
end

function return_settf_dataf_logf(loopdir)
	#this is done for each directory loopdir
	settingsf="";dataf="";logFilename="";
	if isdir(loopdir)	 
		flist=readdir(loopdir)
		# warn("ensure settings and data file exist!!")
		# if alogfile_exists_result_is_already_there(loopdir,flist)
			#do nothing
			#warn("we should write a *.log file which indicates that no model was run....")			
		#	return settingsf,dataf,logFilename
		# else 
		settingsf,dataf=getdata_and_settingsfile(flist)
		if length(dataf)>4
			completePath_wo_end=deepcopy(dataf)[1:end-length(fileending(dataf))]
			logFilename=string(loopdir,completePath_wo_end,".log")
			dataf=string(loopdir,dataf)
			settingsf=string(loopdir,settingsf)			
			if !isfile(logFilename)&&my_iswriteable(dataf)&&my_iswriteable(settingsf)
				return settingsf,dataf,logFilename
			else
				printover("Existing log (or no write access to files): $(logFilename)")
			end
		end				
		#completePath_wo_end=deepcopy(dataf)[1:end-length(fileending(dataf))]
		#logFilename=string(completePath_wo_end,".log")
		#return settingsf,dataf,logFilename
	end
	#else 
	#	return settingsf,dataf,logFilename
	#end # there is no log file and flist has size >0
	return "","",""				
end

function run_these_files(settingsf,dataf,logFilename,julia_code_folder)
	if isfile(settingsf)&&isfile(dataf)
		isfile(logFilename)&&rm(logFilename)
		intHashOfFile=Int128(hash(dataf)) #this is used to ensure that the same file is not submitted again while a process running on this file is already running
		harg=string("-e ",intHashOfFile)	
		#check if this file is already being processedif	
		if isThisFileAlreadyRunningOnAProcess(intHashOfFile)				
			printover("Already running: $(dataf)")			
		else	
			#start a julia session
			outFilename=string("out.",filename(deepcopy(dataf)[1:end-length(fileending(dataf))]))
			#invoke a julia session
			juliaExecutable="C:\\Juno\\resources\\app\\julia\\bin\\julia.exe"
			juliaProgram=string(julia_code_folder,"tryrun.jl")
			arguments=[settingsf dataf outFilename]		
			#run Julia
			process=`$(juliaExecutable) $(juliaProgram) $(harg) $(arguments)`
			# @async 
				remote_reference=run_julia_algorithm(process,logFilename,juliaProgram,harg,arguments)
			printover("Model is running. Check the log file for the status...")			
			#@info "Model run finished."
		end
	else 
		return nothing	
	end#settingsfile and datafile exist	
end

function run_julia_algorithm(process,logFilename,juliaProgram,harg,arguments)
	#version which invokes a separate julia process	
		#@async open(logFilename, "w") do f1; run(process, DevNull, f1, f1); end
	#remote call version with worker process
		rr = RemoteRef()
		rr=invoke_worker_process(rr,process,logFilename,juliaProgram,harg,arguments)
	return rr
	#Old version	
		#pipres=run(pipeline(process, stdout=logFilename, stderr=logFilename))
		#return nothing
end

function invoke_worker_process(rr,process,logFilename,juliaProgram,harg,arguments)	
	println("")
	@info "$(now()). Started: $(logFilename)"
	if nprocs()==1
			rr = this_runs_on_worker(process,logFilename,juliaProgram,harg,arguments)
	else
		@async begin
			rr = @spawn this_runs_on_worker(process,logFilename,juliaProgram,harg,arguments)
		end
		#It is crucial that we wait here until the LOG file is created
		#otherwise we span multiple julia processes which work on the same *.log file
		xn=8
		printover("Modelling finished. Pausing listening for $(xn) seconds.")
		sleep(xn)	
	end	
return rr
end

function this_runs_on_worker(process,logFilename,juliaProgram,harg,arguments)			
		printover("Starting modelling...")
		std_orig=stdout
		stderr_orig=stderr
		f1=open(logFilename, "w") 
			redirect_stdout(f1)
			redirect_stderr(f1)
			this_result=DecisionTrees.run_model(arguments)
		close(f1)
		redirect_stdout(std_orig)
		redirect_stderr(stderr_orig)
		gc();
		#create zip file
			try
				zipit(this_result,arguments,logFilename)		
			catch er
				@show er
				error("failed to zip data")
			end			
		printover("Modelling finished.")
	return this_result
end

function zipit(this_result::Array{T,1},arguments,logFilename) where {T <: AbstractString}
	dataf=arguments[2]
	outFilename=arguments[3]	
	@info "Zipping model results..."
	zipFilename=convert(String,string(dirname(dataf),'\\',outFilename,".7z"))
	#@show this_result
	if !isfile(logFilename)
		warn("Log file not found: $(logFilename)")		
	else
		if this_result==String["false"]	
				tmplogfile=convert(String,string(dirname(dataf),'\\',"failed_log.log"))		
				isfile(tmplogfile)&&rm(tmplogfile);
				cp(logFilename,tmplogfile);
				DecisionTrees.createZipFile(zipFilename,String[logFilename])	
		else
			#just to be on the save side: check if files exist
			flist=String[]
			push!(flist,logFilename)
			for i in this_result 
				if isfile(i)
					if lowercase(i)!=lowercase(zipFilename)
						push!(flist,deepcopy(i))
					end
				end
			end		
			DecisionTrees.createZipFile(zipFilename,flist)
		end
	end
	return nothing
end