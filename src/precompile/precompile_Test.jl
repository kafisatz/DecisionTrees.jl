function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(Tuple{typeof(Test.get_alignment), Test.DefaultTestSet, Int64})
    precompile(Tuple{typeof(Test.testset_beginend), Tuple{String, Expr}, Expr, LineNumberNode})
    precompile(Tuple{typeof(Test.ip_has_file_and_func), Ptr{Nothing}, String, Tuple{Symbol, Symbol}})
    precompile(Tuple{typeof(Test.record), Test.DefaultTestSet, Test.DefaultTestSet})
    precompile(Tuple{typeof(Test.get_testset)})
    precompile(Tuple{typeof(Test.get_test_counts), Test.DefaultTestSet})
    precompile(Tuple{typeof(Test.print_counts), Test.DefaultTestSet, Int64, Int64, Int64, Int64, Int64, Int64, Int64})
    precompile(Tuple{typeof(Test.ip_has_file_and_func), Ptr{Nothing}, String, Tuple{Symbol}})
    precompile(Tuple{typeof(Test.print_test_results), Test.DefaultTestSet, Int64})
    precompile(Tuple{typeof(Test.scrub_backtrace), Array{Union{Ptr{Nothing}, Base.InterpreterIP}, 1}})
    precompile(Tuple{typeof(Test.pop_testset)})
    precompile(Tuple{typeof(Test.eval_test), Expr, Expr, LineNumberNode})
    precompile(Tuple{typeof(Test.get_test_result), Expr, LineNumberNode})
    precompile(Tuple{typeof(Test.push_testset), Test.DefaultTestSet})
    precompile(Tuple{getfield(Test, Symbol("#@testset")), LineNumberNode, Module, Int})
    precompile(Tuple{typeof(Test.get_testset_depth)})
    precompile(Tuple{getfield(Test, Symbol("##13#16")), Base.GenericIOBuffer{Array{UInt8, 1}}})
    precompile(Tuple{typeof(Test.parse_testset_args), Tuple{String}})
    precompile(Tuple{typeof(Test._check_testset), Type{Int}, Expr})
    precompile(Tuple{typeof(Test.record), Test.DefaultTestSet, Test.Fail})
    precompile(Tuple{typeof(Test.filter_errors), Test.DefaultTestSet})
    precompile(Tuple{getfield(Test, Symbol("#@test")), LineNumberNode, Module, Int, Int})
    precompile(Tuple{typeof(Test.ip_has_file_and_func), Base.InterpreterIP, String, Tuple{Symbol, Symbol}})
    precompile(Tuple{typeof(Test.finish), Test.DefaultTestSet})
    precompile(Tuple{typeof(Test.do_test), Test.Returned, Expr})
    precompile(Tuple{typeof(Test.ip_has_file_and_func), Base.InterpreterIP, String, Tuple{Symbol}})
    precompile(Tuple{typeof(Test.test_expr!), String, Expr})
end