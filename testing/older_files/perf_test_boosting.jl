
chosenj="julia200k_200s_jldexists";ARGS=[string("a:\\JuliaTestData\\",chosenj,".settings.csv") string("a:\\JuliaTestData\\",chosenj,".jld") string("out_",chosenj)];

@time for i=1:3
	resbool=DTM.run_model(ARGS)
end
