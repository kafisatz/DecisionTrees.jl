srand(1234)
tic() #for total  runtime

tic() #for the time to load all packages

println("Julia started - Version is $(VERSION) - Time: $(now())");
println("************this is the DEVELOPMENT environment************")
#Loading Packages and workers
	pdir=dirname(@__FILE__)	 # pdir=string(dictGlobalSettings["Julia Code Folder"],"\\src") # cd(pdir)
	println("Initializing... Program folder is: \n $(pdir)")
	@everywhere push!(LOAD_PATH,pdir)
	#if nprocs()>0;include("sendto.jl");sendto(workers(),pdir=pdir);end;
	#@everywhere include(string(pdir,"\\","totalsizeof.jl"))
	#@everywhere include(string(pdir,"\\","init_python.jl"))
	#@everywhere include(string(pdir,"\\","DecisionTrees.jl"))
	#@everywhere using Totalsizeof
	@everywhere using DTM
#report time	
	telapsed=toq()
	println("Time to start Julia and load packages: $(telapsed)")	
#Run the model
  #resbool=DecisionTrees.run_model(ARGS)
  resbool=run_model(ARGS)

telapsed=toq()
println("Total runtime of Julia: $(telapsed/60) minutes")

resbool==true ? println("Julia finished - Time: $(now())") : println("Julia failed -  Time: $(now())")    
resbool=isinteractive()||quit() #exit cmd and Julia if it is not run in REPL mode
