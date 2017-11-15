addprocs(2);
#@everywhere thisdir=isinteractive() ? "R:\\TariffWatch New Repo\\algorithms\\Julia\\Code\\dev_webserver" : dirname(@__FILE__);

debugging=false

#define webserver functions
	@everywhere global const julia_code_folder=string(dictGlobalSettings["TW SVN Repo"],"algorithms\\Julia\\Code\\dev\\src\\");	
	@everywhere include(string(julia_code_folder,"webserver_fn.jl"));
#define julia module
	@everywhere push!(LOAD_PATH,julia_code_folder);
#load DTM on all processes
	@everywhere using DTM
#set folder list file
	folder_list_file=string(dictGlobalSettings["TW SVN Repo"],"algorithms\\Julia\\Code\\testJLfiles\\webserver_listening_directories.csv");
#test if module was loaded sucessfully
	debugging ? testresult=true : testresult=test_module_on_workers()
	@assert testresult	
	
if testresult 
	# listening_folders=["c:\\temp\\myf","c:\\temp\\myother","c:\\temp\\mixi\\","c:\\temp\\1\\","c:\\temp\\2\\","c:\\temp\\3\\","c:\\temp\\4\\"];
	startlistening(folder_list_file,julia_code_folder);
end