number_of_threads=8 #needs to be equal to the number of cores (or processes which will "work")
addprocs(number_of_threads);
@everywhere this_dir=string(dictGlobalSettings["Julia Code Folder"],"\\src");
@everywhere push!(LOAD_PATH,this_dir)
using DTM

chosenj="boosting_minimalexample_multirow_settings";ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];
res=DTM.run_model(ARGS)

info("Finished multicore run,")

chosenj="boostmed_multirow";ARGS=[string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".settings.csv") string(dictGlobalSettings["Data Folder"],"\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];
res=DTM.run_model(ARGS)

ARGS=["e:\\temp\\xk\\boost_multir_test.settings.csv" "e:\\temp\\xk\\boost_multir_test.jld" "boost_multirow_out"];
res=DTM.run_model(ARGS)
