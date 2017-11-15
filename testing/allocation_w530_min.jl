function my_sin(x)
zk=2*sin(x)
for i=1:99999999
	zk=2*sin(x)
end
return zk
end

@show my_sin(2.991)
my_sin(9)

Profile.clear()          # in case we have any previous profiling data
Profile.init()
#Profile.init(10^7, 0.01) #n,delay #Profiling settings

Profile.clear_malloc_data()
juliaresult=my_sin(3.221)

# Profile.print()

#exit()

#=

Executing .juliarc.jl... done
my_sin(2.991) = 0.3000482098478442

signal (11): SIGSEGV
while loading no file, in expression starting on line 0
crt_sig_handler at /home/Administrator/buildbot/slave/package_win6_2-x64/build/src/home/Administrator/bui
ldbot/slave/package_win6_2-x64/build/src\signals-win.c:89
_gnu_exception_handler at /usr/src/debug/mingw64-x86_64-runtime-4.0.6-1/crt\crt_handler.c:223
_C_specific_handler at C:\Windows\SYSTEM32\ntdll.dll (unknown line)
RtlDecodePointer at C:\Windows\SYSTEM32\ntdll.dll (unknown line)
RtlUnwindEx at C:\Windows\SYSTEM32\ntdll.dll (unknown line)
KiUserExceptionDispatcher at C:\Windows\SYSTEM32\ntdll.dll (unknown line)
ZNSt6localeaSERKS_ at C:\JuliaPro-0.5.1.1\Julia-0.5.1\bin\libstdc++-6.dll (unknown line)
ZNSt8ios_base7_M_initEv at C:\JuliaPro-0.5.1.1\Julia-0.5.1\bin\libstdc++-6.dll (unknown line)
ZNSt9basic_iosIcSt11char_traitsIcEE4initEPSt15basic_streambufIcS1_E at C:\JuliaPro-0.5.1.1\Julia-0.5.1\bi
n\libstdc++-6.dll (unknown line)
basic_ifstream at /home/Administrator/buildbot/slave/package_win6_2-x64/build/src/usr/lib/gcc/x86_64-w64-
mingw32/4.9.2/include/c++\fstream:473 [inlined]
write_log_data at /home/Administrator/buildbot/slave/package_win6_2-x64/build/src/home/Administrator/buil
dbot/slave/package_win6_2-x64/build/src\codegen.cpp:1448
jl_atexit_hook at /home/Administrator/buildbot/slave/package_win6_2-x64/build/src/home/Administrator/buil
dbot/slave/package_win6_2-x64/build/src\init.c:239
wmain at /home/Administrator/buildbot/slave/package_win6_2-x64/build/ui\repl.c:244
__tmainCRTStartup at /usr/src/debug/mingw64-x86_64-runtime-4.0.6-1/crt\crtexe.c:329
mainCRTStartup at /usr/src/debug/mingw64-x86_64-runtime-4.0.6-1/crt\crtexe.c:212
BaseThreadInitThunk at C:\Windows\system32\kernel32.dll (unknown line)
RtlUserThreadStart at C:\Windows\SYSTEM32\ntdll.dll (unknown line)
Allocations: 961961 (Pool: 961089; Big: 872); GC: 0

=#

#=
this_dir=string(dictGlobalSettings["Julia Code Folder"],"\\src");cd(this_dir);push!(LOAD_PATH,this_dir)
using DTM
#using Graphics,ProfileView

#this runs a small version to compile all functions
chosenj="boosting_small";ARGS=[string("a:\\JuliaTestData\\",chosenj,".settings.csv") string("a:\\JuliaTestData\\",chosenj,".csv") string("out_",chosenj)];
resbool=DTM.run_model(ARGS) #small version to compile functions	

chosenj="julia200k";ARGS=[string("a:\\JuliaTestData\\",chosenj,".settings.csv") string("a:\\JuliaTestData\\",chosenj,".jld") string("out_",chosenj)];
juliaresult=DTM.run_model(ARGS)   

Profile.clear_malloc_data()
chosenj="julia200k";ARGS=[string("a:\\JuliaTestData\\",chosenj,".settings.csv") string("a:\\JuliaTestData\\",chosenj,".jld") string("out_",chosenj)];
juliaresult=DTM.run_model(ARGS)   
=#

