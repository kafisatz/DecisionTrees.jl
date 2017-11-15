this_dir=string(dictGlobalSettings["Julia Code Folder"],"\\src");cd(this_dir);push!(LOAD_PATH,this_dir)
chosenj="tree_small" #this runs a small version to compile all functions
using DTM
#using Graphics,ProfileView
using ProfileView

#this runs a small version to compile all functions
loc="C:\\Users\\tbd\\Documents\\Arbeitsordner (sync)\\Public NL\\Personal\\Bernhard\\Resources\\Milliman Software\\Predictive Analytics Paris\\Training\\data\\"
loc=string(dictGlobalSettings["Data Folder"],"JuliaTestData\\")
ARGS=[string(loc,"jmortgage.settings_profiling.csv") string(loc,"jmortgage.jld") string("out_jmortgage")];
resbool=DTM.run_model(ARGS) #small version to compile functions
@time resbool=DTM.run_model(ARGS)

#chosenj="julia200k";ARGS=[string("a:\\JuliaTestData\\",chosenj,".settings.csv") string("a:\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];
#@time resbool=DTM.run_model(ARGS)
#chosenj="julia200k_200s";ARGS=[string("a:\\JuliaTestData\\",chosenj,".settings2.csv") string("a:\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];
#juliaresult=DTM.run_model(ARGS)
chosenj="julia400k";ARGS=[string("a:\\JuliaTestData\\",chosenj,".settings2.csv") string("a:\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];
@time juliaresult=DTM.run_model(ARGS)
chosenj="julia400k";ARGS=[string("a:\\JuliaTestData\\",chosenj,".settings2.csv") string("a:\\JuliaTestData\\",chosenj,".jld") string("out_",chosenj)];
@time juliaresult=DTM.run_model(ARGS)

println("**************************************************\n\n\nStarting Profiling....")


Profile.clear()          # in case we have any previous profiling data
Profile.init(10^7, 0.01) #n,delay #Profiling settings
@profile juliaresult=DTM.run_model(ARGS)   # Profile the program
# profile() #this shows profile graphically in Juno

val=ProfileView.view()

#=
	@profile df=DTM.mysql_execute_query_DF_wo_NA(con,command);
	val=ProfileView.view()

	@profile df=DTM.mysql_execute_query(con,command);
	val=ProfileView.view()

	Profile.clear()          # in case we have any previous profiling data
	h=x*.23+9
	@profile juliaresult=DTM.gini(x,h)   # Profile the program
	val=ProfileView.view()
# profile() #this shows profile graphically in Juno

val=ProfileView.view()

=#
#=
	val2=ProfileView.view(Profile.fetch(); format = :tree,lidict = nothing, C = false, colorgc = true, fontsize = 12, combine = true)

	fi="c:\\temp\\myprof.txt"
	isfile(fi)&&rm(fi)
	f=open(fi,"w");
		Profile.print(f,format=:flat,combine=true,cols=1000,combine=true);
	close(f);
	f=open(fi,"w")
		Profile.print(f,format=:flat,combine=true,cols=1000,combine=false)
	close(f)
	d,l=Profile.retrieve();
=#
