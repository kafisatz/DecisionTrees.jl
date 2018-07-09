function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(Tuple{typeof(DecisionTrees.poissonError), Array{Float64, 1}, Array{Float64, 1}})
    precompile(Tuple{typeof(DecisionTrees.subset_splitlist), Array{DecisionTrees.Splitdef{UInt8}, 1}, Float64})
    
    precompile(Tuple{typeof(DecisionTrees.absrel), Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}})
    precompile(Tuple{typeof(DecisionTrees.aggregate_data), DecisionTrees.PooledArraysDTM.PooledArray{String, UInt8, 1, Array{UInt8, 1}}, Base.SubArray{Int, 1, Array{Int, 1}, Tuple{Array{Int, 1}}, false}, Array{Float64, 1}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}})
    precompile(Tuple{typeof(DecisionTrees.aggregate_data), DecisionTrees.PooledArraysDTM.PooledArray{Float64, UInt8, 1, Array{UInt8, 1}}, Base.SubArray{Int, 1, Array{Int, 1}, Tuple{Array{Int, 1}}, false}, Array{Float64, 1}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}})
    precompile(Tuple{typeof(DecisionTrees.calculateSplitValue), DecisionTrees.DifferenceSplit, Symbol, Int, Array{UInt8, 1}, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}, Float64, DecisionTrees.MyGrayCodeSubsetsHALF})
    precompile(Tuple{typeof(DecisionTrees.InsertionSort!), Array{Float64, 1}, Array{Int, 1}, Int, Int})
    precompile(Tuple{typeof(DecisionTrees.sas_write_elseif_ScoreMap), Base.IOStream, String, Int, Array{Int, 1}, Array{Float64, 1}, Array{Float64, 1}})
    precompile(Tuple{typeof(DecisionTrees.errorhistogram), Array{Float64, 1}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}})
    precompile(Tuple{typeof(DecisionTrees.write_tree), Array{Array{Float64, 1}, 1}, DecisionTrees.BoostedTree, Int, Int, String, Array{String, 1}, Array{Array{String, 1}, 1}})
    precompile(Tuple{typeof(DecisionTrees.boosted_tree), DecisionTrees.DTMTable, DecisionTrees.ModelSettings})
    precompile(Tuple{typeof(DecisionTrees.buildStatistics), Array{String, 2}, Array{Any, 1}, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}})
    precompile(Tuple{typeof(DecisionTrees.removeUnionTypes!), DataFrames.DataFrame, Array{String, 1}})
    precompile(Tuple{typeof(DecisionTrees.run_model_actual), DecisionTrees.DTMTable, DecisionTrees.ModelSettings, String})
    precompile(Tuple{typeof(DecisionTrees.derive_scores_main_aggregation_step!), Int, Float64, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Array{Int, 1}, Array{Float64, 1}})
    precompile(Tuple{typeof(DecisionTrees.find_max_type), DataFrames.DataFrame})
    precompile(Tuple{typeof(DecisionTrees.gini_single_argument), Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}})
    precompile(Tuple{typeof(DecisionTrees.vba_write_elseif_ScoreMap), Base.IOStream, Int, Array{Int, 1}, Array{Float64, 1}, Array{Float64, 1}})
    
    precompile(Tuple{typeof(DecisionTrees.randint), Int})
    precompile(Tuple{typeof(DecisionTrees.defineTwoWayCharts), String, String, String, String, Int, Int, Int, Int, Int, String, String, String})
    precompile(Tuple{typeof(DecisionTrees.createTwoWayValidationCharts), Array{Int, 1}, Array{Int, 1}, Int, Array{String, 1}, Array{Array{String, 1}, 1}, Array{Array{Float64, 1}, 1}, DecisionTrees.ModelSettings, Array{Int, 1}, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}, DataFrames.DataFrame})
    precompile(Tuple{typeof(DecisionTrees.determine_used_variables), DecisionTrees.BoostedTree})
    precompile(Tuple{typeof(DecisionTrees.buildStatisticsInternal), Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}})
    precompile(Tuple{getfield(DecisionTrees, Symbol("#kw##defineRelativityChart")), NamedTuple{(:headerrow2, :headercol2, :datarow2, :xtitle, :ytitle, :xscale, :title, :datarow, :valuescol, :categoriescol, :valuescol2, :yscale), Tuple{Int, Int, Int, String, String, Float64, String, Int, Int, Int, Int, Float64}}, typeof(DecisionTrees.defineRelativityChart), String, String, String, Int, Int, Int})
    precompile(Tuple{typeof(DecisionTrees.correctSmallNegativeValues!), Array{Float64, 1}})
    precompile(Tuple{typeof(DecisionTrees.vba_get_signature), Array{Array{String, 1}, 1}, Array{String, 1}, Int})
    precompile(Tuple{typeof(DecisionTrees.calc_sum_squares), Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}})
    precompile(Tuple{getfield(DecisionTrees, Symbol("#kw##defineChartWithNSeries")), NamedTuple{(:categoriescol, :yaxisformat, :xscale), Tuple{Int, String, Float64}}, typeof(DecisionTrees.defineChartWithNSeries), String, String, String, String, Int, Int, Array{Int, 1}, Int, String, String, String})
    precompile(Tuple{typeof(DecisionTrees.calcErrors), Float64, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}})
    precompile(Tuple{typeof(DecisionTrees.csharp_write_elseif_ScoreMap), Base.IOStream, Array{Int, 1}, Array{Float64, 1}, Array{Float64, 1}})
    precompile(Tuple{typeof(DecisionTrees.addPredictorData), Array{String, 1}, Array{String, 2}, DecisionTrees.ModelSettings, Base.SubArray{Int, 1, Array{Int, 1}, Tuple{Array{Int, 1}}, false}, Array{Float64, 1}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, DecisionTrees.PooledArraysDTM.PooledArray{String, UInt8, 1, Array{UInt8, 1}}, Int})
    precompile(Tuple{typeof(DecisionTrees.addPredictorData), Array{Float64, 1}, Array{String, 2}, DecisionTrees.ModelSettings, Base.SubArray{Int, 1, Array{Int, 1}, Tuple{Array{Int, 1}}, false}, Array{Float64, 1}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, DecisionTrees.PooledArraysDTM.PooledArray{Float64, UInt8, 1, Array{UInt8, 1}}, Int})
    precompile(Tuple{typeof(DecisionTrees.addTariffEstimationStatsAndGraphs!), DecisionTrees.ExcelData, Array{Int, 1}, Array{Int, 1}, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 2}, Array{Float64, 2}})
    precompile(Tuple{typeof(DecisionTrees.sampleData!), Array{Int, 1}, Int, Array{Int, 1}})
    precompile(Tuple{typeof(DecisionTrees.write_tree_at_each_node!), Array{Array{Float64, 1}, 1}, DecisionTrees.Leaf, Int, Int, Base.IOStream, Array{String, 1}, Array{Array{String, 1}, 1}, String, Float64})
    precompile(Tuple{typeof(DecisionTrees.csharp_write_pub_str), Base.IOStream, String, String})
    precompile(Tuple{typeof(DecisionTrees.gini), Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Array{Float64, 1}})
    precompile(Tuple{typeof(DecisionTrees.defineScoreChart), String, String, String, Int, Int, Int, Int, Int})
    precompile(Tuple{typeof(DecisionTrees.fileroot), Base.Missing})
    precompile(Tuple{typeof(DecisionTrees.poissonError), Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Array{Float64, 1}})
    precompile(Tuple{typeof(DecisionTrees.createTrnValStats!), Array{Int, 1}, Array{Int, 1}, DecisionTrees.ModelSettings, String, String, DecisionTrees.ExcelData, DecisionTrees.Tree, Array{DecisionTrees.Leaf, 1}, Array{Float64, 1}, Array{Int, 1}, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}})
    precompile(Tuple{typeof(DecisionTrees.aggregate_values_per_score), Int, Array{Int, 1}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Int, Array{Float64, 1}})
    precompile(Tuple{getfield(DecisionTrees, Symbol("#kw##defineChartWithNSeries")), NamedTuple{(:categoriescol, :xaxisformat, :yaxisformat, :yscale, :xscale), Tuple{Int, String, String, Float64, Float64}}, typeof(DecisionTrees.defineChartWithNSeries), String, String, String, String, Int, Int, Array{Int, 1}, Int, String, String, String})
    precompile(Tuple{typeof(DecisionTrees.bagged_tree), DecisionTrees.DTMTable, DecisionTrees.ModelSettings})
    precompile(Tuple{typeof(DecisionTrees.initSettingsWhichAreTheSameForBoostingAndBagging), Array{Int, 1}, Array{Int, 1}, Array{Float64, 1}, Array{Float64, 1}, DecisionTrees.ModelSettings})
    precompile(Tuple{typeof(DecisionTrees.printover), String})
    precompile(Tuple{typeof(DecisionTrees.nodesize), DecisionTrees.Tree})
    precompile(Tuple{typeof(DecisionTrees.update_and_derive_scores!), Array{Int, 1}, Array{Float64, 1}, Array{Float64, 1}})
    precompile(Tuple{typeof(DecisionTrees.write_and_create_boollist_char_vba), Base.IOStream, Array{Array{String, 1}, 1}, Array{String, 1}, Int, Array{Array{Bool, 1}, 1}})
    precompile(Tuple{getfield(DecisionTrees, Symbol("##writeSasCode#38")), String, Int, typeof(identity), Array{Array{Float64, 1}, 1}, DecisionTrees.Tree, Int, String, Array{String, 1}, String, Array{Array{String, 1}, 1}, Float64})
    precompile(Tuple{typeof(DecisionTrees.__init__)})
    precompile(Tuple{typeof(DecisionTrees.write_and_create_boollist_num_vba), Base.IOStream, Array{String, 1}, Int, Array{Array{Float64, 1}, 1}, Array{Array{Bool, 1}, 1}})
    precompile(Tuple{typeof(DecisionTrees.defineUnivariateChartWith2Lines), String, String, String, String, Int, Int, Int, Int, Int, Int})
    precompile(Tuple{typeof(DecisionTrees.customInverseNormalized!), Array{Float64, 1}})
    precompile(Tuple{typeof(DecisionTrees.update_matched_vectors!), Array{Int, 1}, Array{Int, 1}, Base.BitArray{1}, Base.BitArray{1}, Array{UInt8, 1}, UInt8})
    precompile(Tuple{typeof(DecisionTrees.write_tree), Array{Array{Float64, 1}, 1}, DecisionTrees.Leaf, Int, Int, String, Array{String, 1}, Array{Array{String, 1}, 1}})
    precompile(Tuple{typeof(DecisionTrees.cleanString), String})
    
    precompile(Tuple{typeof(DecisionTrees.add_iteration_charts!), DecisionTrees.ExcelData, DecisionTrees.ModelSettings, Int})
    precompile(Tuple{typeof(DecisionTrees.vba_write_booleans), Base.IOStream, Array{Array{String, 1}, 1}, Array{Array{String, 1}, 1}, Array{String, 1}, Array{Array{Float64, 1}, 1}, Array{Array{String, 1}, 1}, Array{Array{Bool, 1}, 1}, Array{Array{Bool, 1}, 1}, Array{Bool, 1}})
    precompile(Tuple{typeof(DecisionTrees.pushsum!), Array{Float64, 1}})
    precompile(Tuple{typeof(DecisionTrees.excelLetter), Int})
    precompile(Tuple{typeof(DecisionTrees.csharp_write_ScoreErrorResponseCodeCheck), Base.IOStream, Array{Array{String, 1}, 1}, Array{Array{String, 1}, 1}, Array{String, 1}, Array{Array{Float64, 1}, 1}, Array{Array{String, 1}, 1}, Array{Array{Bool, 1}, 1}, Array{Array{Bool, 1}, 1}, Array{Bool, 1}})
    precompile(Tuple{typeof(DecisionTrees.write_tree), Array{Array{Float64, 1}, 1}, DecisionTrees.Leaf, Int, Array{String, 2}, Array{String, 2}, Int, Base.IOStream, Array{String, 1}, Array{Array{String, 1}, 1}})
    
    
    precompile(Tuple{typeof(DecisionTrees.calcIterationsPerCore), Int, Int})
    precompile(Tuple{typeof(DecisionTrees.insert_gaps_to_scores!), Float64, Int, Int, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}, Array{Int, 1}, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}})
    precompile(Tuple{typeof(DecisionTrees.wAverage), Array{Float64, 1}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Base.UnitRange{Int}}, true}})
    precompile(Tuple{typeof(DecisionTrees.build_listOfMeanResponse), DecisionTrees.DifferenceSplit, Array{Int, 1}, Array{Int, 1}, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}, Base.SubArray{String, 1, DecisionTrees.PooledArraysDTM.PooledArray{String, UInt8, 1, Array{UInt8, 1}}, Tuple{Array{Int, 1}}, false}, Array{UInt8, 1}, Float64})
    
    precompile(Tuple{typeof(DecisionTrees.variable_importance), Array{DecisionTrees.Leaf, 1}, Array{String, 1}, Int})
    precompile(Tuple{typeof(DecisionTrees.build_tree_iteration!), Array{Int, 1}, Array{Int, 1}, DecisionTrees.ModelSettings, DecisionTrees.Tree, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}, DataFrames.DataFrame, Int, Float64, Array{DecisionTrees.Rulepath{T} where T<:Unsigned, 1}, Bool, Int, Array{Float64, 1}})
    precompile(Tuple{typeof(DecisionTrees.checkIfSettingsAreValid), DecisionTrees.ModelSettings})
    precompile(Tuple{typeof(DecisionTrees._split), UInt8, Int, Array{Int, 1}, Array{Int, 1}, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}, Array{Symbol, 1}, DataFrames.DataFrame, Float64, Int, Float64, DecisionTrees.DifferenceSplit, Int, Int, Array{Int, 1}, Int, DecisionTrees.SortByMean})
    precompile(Tuple{typeof(DecisionTrees.defineChartWithNSeries), String, String, String, String, Int, Int, Array{Int, 1}, Int, String, String, String})
    precompile(Tuple{typeof(DecisionTrees.lrIndices), Array{Int, 1}, DecisionTrees.PooledArraysDTM.PooledArray{String, UInt8, 1, Array{UInt8, 1}}, Array{UInt8, 1}})
    precompile(Tuple{typeof(DecisionTrees.tariffEstStats), Float64, Array{Float64, 1}, Array{Float64, 2}, Float64, Array{Int, 1}, Array{Int, 1}, String})
    precompile(Tuple{typeof(DecisionTrees.createTrnValStatsForThisIteration), Array{String, 1}, Int, Array{Int, 1}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Int, 1, Array{Int, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Int, 1, Array{Int, 1}, Tuple{Array{Int, 1}}, false}, DecisionTrees.ModelSettings})
    
    precompile(Tuple{typeof(DecisionTrees.enforce_monotonicity!), Array{Float64, 1}})
    precompile(Tuple{typeof(DecisionTrees.variablesUsed), DecisionTrees.BoostedTree})
    precompile(Tuple{typeof(DecisionTrees.defineUnivariateChart), String, String, String, String, Int, Int, Int, Int, Int})
    precompile(Tuple{typeof(DecisionTrees.countsort!), Array{UInt8, 1}})
    precompile(Tuple{typeof(DecisionTrees.lrIndices), Array{Int, 1}, DecisionTrees.PooledArraysDTM.PooledArray{Float64, UInt8, 1, Array{UInt8, 1}}, Array{UInt8, 1}})
    precompile(Tuple{typeof(DecisionTrees.vba_write_writeIterations_recursive), Float64, Int, Base.IOStream, DecisionTrees.Leaf, Array{Array{String, 1}, 1}, Array{Array{String, 1}, 1}, Array{String, 1}, Array{Array{Float64, 1}, 1}, Array{Array{String, 1}, 1}, Int})
    precompile(Tuple{typeof(DecisionTrees.tariffEstStats), Float64, Array{Float64, 1}, Array{Float64, 1}, Float64, Array{Int, 1}, Array{Int, 1}, String, Int})
    precompile(Tuple{typeof(DecisionTrees.aggregateScores), Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Int, 1, Array{Int, 1}, Tuple{Array{Int, 1}}, false}, Array{Int, 1}})
    precompile(Tuple{typeof(DecisionTrees.sum_up_ratios!), Array{Float64, 1}, Array{Float64, 1}, Int, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Int, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Int, Int})
    precompile(Tuple{typeof(DecisionTrees.createPredictorData), Array{Int, 1}, String, Array{Array{String, 1}, 1}, Array{Array{Float64, 1}, 1}, DecisionTrees.ModelSettings, Array{Int, 1}, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}, DataFrames.DataFrame})
    
    precompile(Tuple{typeof(DecisionTrees.csharp_write_RequiredElementsProvided), Base.IOStream, Array{String, 1}, Array{Bool, 1}})
    precompile(Tuple{typeof(DecisionTrees.write_tree), Array{Array{Float64, 1}, 1}, DecisionTrees.Node{UInt8}, Int, Array{String, 2}, Array{String, 2}, Int, Base.IOStream, Array{String, 1}, Array{Array{String, 1}, 1}})
    precompile(Tuple{typeof(DecisionTrees.writeStatisticsFile!), String, DecisionTrees.ExcelData, Array{String, 1}})
    precompile(Tuple{typeof(DecisionTrees.print_some_settings), DecisionTrees.ModelSettings, Array{String, 1}})
    precompile(Tuple{typeof(DecisionTrees.calculateSplitValue), DecisionTrees.DifferenceSplit, Symbol, Int, Array{UInt8, 1}, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}, Float64, DecisionTrees.MyIncreasingSubsets})
    precompile(Tuple{typeof(DecisionTrees.lowessSmoothVector!), Array{Float64, 1}, Float64})
    precompile(Tuple{typeof(DecisionTrees.csharp_write_writeIterations), Int, Base.IOStream, Int, DecisionTrees.BoostedTree, Array{Array{String, 1}, 1}, Array{Array{String, 1}, 1}, Array{String, 1}, Array{Array{Float64, 1}, 1}, Array{Array{String, 1}, 1}, Int})
    precompile(Tuple{typeof(DecisionTrees.map_numdata_to_candidates), Array{Float64, 1}, Array{Float64, 1}})
    precompile(Tuple{typeof(DecisionTrees.addChartToWorkbook!), PyCall.PyCall.PyObject, PyCall.PyCall.PyObject, Base.Dict{AbstractString, Base.Dict{AbstractString, Any}}, String})
    precompile(Tuple{typeof(DecisionTrees.write_tree), Array{Array{Float64, 1}, 1}, DecisionTrees.Node{UInt8}, Int, Int, Base.IOStream, Array{String, 1}, Array{Array{String, 1}, 1}})
    precompile(Tuple{typeof(DecisionTrees.write_and_create_boollist_char), Base.IOStream, Array{Array{String, 1}, 1}, Array{String, 1}, Int, Array{Array{Bool, 1}, 1}})
    precompile(Tuple{typeof(DecisionTrees.mylinreg), Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}})
    precompile(Tuple{typeof(DecisionTrees.write_tree_at_each_node!), Array{Array{Float64, 1}, 1}, DecisionTrees.Node{UInt8}, Int, Int, Base.IOStream, Array{String, 1}, Array{Array{String, 1}, 1}, String, Float64})
    precompile(Tuple{typeof(DecisionTrees.graph), DecisionTrees.Node{UInt8}, Int, String, DecisionTrees.ModelSettings, String, Array{Array{String, 1}, 1}, Array{Array{Float64, 1}, 1}, Float64})
    precompile(Tuple{typeof(DecisionTrees.create_custom_string), String})
    precompile(Tuple{typeof(DecisionTrees.apply_tree_by_leaf_iteration!), Array{Int, 1}, DecisionTrees.Leaf, DataFrames.DataFrame, Array{Float64, 1}, Array{Int, 1}})
    precompile(Tuple{typeof(DecisionTrees.csharp_write_GenerateScore), Base.IOStream, DecisionTrees.BoostedTree})
    precompile(Tuple{typeof(DecisionTrees.variable_importance_internal), Array{DecisionTrees.Leaf, 1}, Array{String, 1}, Int})
    precompile(Tuple{typeof(DecisionTrees.gini_single_argument), Array{Float64, 1}})
    precompile(Tuple{typeof(DecisionTrees.defineCandidates), Array{Float64, 1}, Int})
    precompile(Tuple{typeof(DecisionTrees.graph), DecisionTrees.Leaf, Int, String, DecisionTrees.ModelSettings, String, Int, Array{Array{Float64, 1}, 1}, Float64})
    precompile(Tuple{typeof(DecisionTrees.cumulativeToIncremental!), Array{Int, 1}})
    precompile(Tuple{getfield(DecisionTrees, Symbol("#kw##defineRelativityChart")), NamedTuple{(:xscale,), Tuple{Float64}}, typeof(DecisionTrees.defineRelativityChart), String, String, String, Int, Int, Int})
    precompile(Tuple{typeof(DecisionTrees.map_these_values), Array{Float64, 1}, Array{Int, 1}})
    precompile(Tuple{typeof(DecisionTrees.graph), DecisionTrees.Tree})
    precompile(Tuple{typeof(DecisionTrees.calcErrors), Float64, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Float64, 1, Array{Float64, 2}, Tuple{Array{Int, 1}, Int}, false}})
    precompile(Tuple{typeof(DecisionTrees.total_n_subnodes), DecisionTrees.Node{UInt8}})
    precompile(Tuple{typeof(DecisionTrees.apply_tree_by_leaf_iteration!), Array{Int, 1}, DecisionTrees.Node{UInt16}, DataFrames.DataFrame, Array{Float64, 1}, Array{Int, 1}})
    precompile(Tuple{typeof(DecisionTrees.set_leaf_numbers!), DecisionTrees.Tree})
    precompile(Tuple{typeof(DecisionTrees.some_tree_settings), Array{Int, 1}, Array{Int, 1}, Array{Int, 1}, Array{Array{Float64, 1}, 1}, Array{Array{String, 1}, 1}, Float64, Array{Float64, 1}, Float64, Int})
    precompile(Tuple{typeof(DecisionTrees.numeratorsum), DecisionTrees.Node{UInt8}})
    precompile(Tuple{typeof(DecisionTrees.writeSasCode), Array{Float64, 1}, Array{Array{Float64, 1}, 1}, DecisionTrees.BoostedTree, Int, String, Array{String, 1}, String, Array{Array{String, 1}, 1}})
    precompile(Tuple{typeof(DecisionTrees.cumulativeToIncremental!), Array{Float64, 1}})
    precompile(Tuple{typeof(DecisionTrees.apply_tree_by_leaf_iteration!), Array{Int, 1}, DecisionTrees.Node{UInt8}, DataFrames.DataFrame, Array{Float64, 1}, Array{Int, 1}})
    precompile(Tuple{typeof(DecisionTrees.writeCsharpCode), Array{Array{DecisionTrees.Leaf, 1}, 1}, Array{Float64, 1}, Array{Array{Float64, 1}, 1}, DecisionTrees.BoostedTree, String, String, Array{Array{String, 1}, 1}, Int, DecisionTrees.ModelSettings})
    precompile(Tuple{typeof(DecisionTrees.quicksort!), Array{Float64, 1}, Array{Int, 1}, Int, Int})
    precompile(Tuple{typeof(DecisionTrees.createZipFile), String, Array{String, 1}, Bool})
    precompile(Tuple{typeof(DecisionTrees.write_and_create_boollist_num), Base.IOStream, Array{String, 1}, Int, Array{Array{Float64, 1}, 1}, Array{Array{Bool, 1}, 1}})
    precompile(Tuple{typeof(DecisionTrees.writeDFtoExcel), DecisionTrees.ExcelData, String, Int, Int, Bool, Bool})
    precompile(Tuple{typeof(DecisionTrees.defineWeightPerTree), Array{DecisionTrees.TreeWithErrorStats, 1}, String})
    precompile(Tuple{typeof(DecisionTrees.randomFeatureSelection), Int, Float64})
    precompile(Tuple{typeof(DecisionTrees.create_leaves_array), DecisionTrees.Tree})
    precompile(Tuple{typeof(DecisionTrees.bitflip_graycode_subsetsHALF), Array{UInt8, 1}})
    precompile(Tuple{typeof(DecisionTrees.apply_tree_by_leaf!), Array{Float64, 1}, Array{Int, 1}, Array{Int, 1}, DecisionTrees.Node{UInt8}, DataFrames.DataFrame})
    precompile(Tuple{typeof(DecisionTrees.createIntegerIndices), Int})
    precompile(Tuple{typeof(DecisionTrees.sortlists!), DecisionTrees.SortByMean, Array{UInt8, 1}, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}})
    precompile(Tuple{typeof(DecisionTrees.get_feature_pools), DataFrames.DataFrame})
    precompile(Tuple{typeof(DecisionTrees.build_listOfMeanResponse), DecisionTrees.DifferenceSplit, Array{Int, 1}, Array{Int, 1}, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}, Base.SubArray{Float64, 1, DecisionTrees.PooledArraysDTM.PooledArray{Float64, UInt8, 1, Array{UInt8, 1}}, Tuple{Array{Int, 1}}, false}, Array{UInt8, 1}, Float64})
    precompile(Tuple{typeof(DecisionTrees.nodesize), DecisionTrees.Node{UInt8}})
    precompile(Tuple{typeof(DecisionTrees.denominatorsum), DecisionTrees.Node{UInt8}})
    precompile(Tuple{typeof(DecisionTrees._split_feature), UInt8, Int, Array{Int, 1}, Array{Int, 1}, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}, Symbol, Int, Float64, DecisionTrees.DifferenceSplit, Int, Float64, Int, DecisionTrees.SortByMean})
    precompile(Tuple{typeof(DecisionTrees.smooth_scores), Array{Float64, 1}, Bool, Bool})
    precompile(Tuple{typeof(DecisionTrees.number_of_nodes), DecisionTrees.Node{UInt8}})
    precompile(Tuple{typeof(DecisionTrees.check_for_missing_data), DataFrames.DataFrame, Int})
    precompile(Tuple{typeof(DecisionTrees.build_tree!), Array{Int, 1}, Array{Int, 1}, Array{Array{Float64, 1}, 1}, Array{Array{String, 1}, 1}, DecisionTrees.ModelSettings, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}, DataFrames.DataFrame, Array{Float64, 1}})
    precompile(Tuple{typeof(DecisionTrees.maxdepth), DecisionTrees.Node{UInt8}})
    precompile(Tuple{typeof(DecisionTrees.myresort!), Array{UInt8, 1}, Array{Int, 1}})
    precompile(Tuple{typeof(DecisionTrees.myresort!), Array{Float64, 1}, Array{Int, 1}})
    precompile(Tuple{typeof(DecisionTrees.writeVbaCode), Array{Array{DecisionTrees.Leaf, 1}, 1}, Array{Float64, 1}, Array{Array{Float64, 1}, 1}, DecisionTrees.BoostedTree, String, String, Array{Array{String, 1}, 1}, Int, DecisionTrees.ModelSettings})
    precompile(Tuple{typeof(DecisionTrees.number_of_nodes), DecisionTrees.Tree})
    precompile(Tuple{typeof(DecisionTrees.increasing_subsets), Array{UInt8, 1}})
    precompile(Tuple{typeof(DecisionTrees.create_custom_dict), DataFrames.DataFrame})
    
    precompile(Tuple{typeof(DecisionTrees.sort_splitlist), Array{DecisionTrees.Splitdef{UInt8}, 1}})
    precompile(Tuple{typeof(DecisionTrees.is_categorical_column), Array{Int, 1}, Symbol})
    precompile(Tuple{typeof(DecisionTrees.writeStatistics), DecisionTrees.ExcelData, String, Bool, Bool})
    precompile(Tuple{typeof(DecisionTrees.is_categorical_column), Array{String, 1}, Symbol})
    precompile(Tuple{typeof(DecisionTrees.maxdepth), DecisionTrees.Tree})
    precompile(Tuple{typeof(DecisionTrees.writeAllFieldsToArray), DecisionTrees.ModelSettings})
    precompile(Tuple{typeof(DecisionTrees.fileroot), PyCall.PyCall.PyObject})
    precompile(Tuple{typeof(DecisionTrees.vba_write_writeIterations), Int, Base.IOStream, Int, DecisionTrees.BoostedTree, Array{Array{String, 1}, 1}, Array{Array{String, 1}, 1}, Array{String, 1}, Array{Array{Float64, 1}, 1}, Array{Array{String, 1}, 1}, Int})
    precompile(Tuple{typeof(DecisionTrees.write_tree), Array{Array{Float64, 1}, 1}, DecisionTrees.Tree, Int, Array{String, 2}, Array{String, 2}, Int, String, Array{String, 1}, Array{Array{String, 1}, 1}})
    precompile(Tuple{typeof(DecisionTrees.vba_write_ScoreMap), Base.IOStream, Array{Float64, 1}, Array{Float64, 1}})
    precompile(Tuple{typeof(DecisionTrees.normalized_gini), Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}})
    precompile(Tuple{typeof(DecisionTrees.sample_data_and_build_tree!), Array{Int, 1}, Array{Int, 1}, Array{Float64, 1}, Array{Array{Float64, 1}, 1}, Array{Array{String, 1}, 1}, DecisionTrees.ModelSettings, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}, DataFrames.DataFrame, Int, Int, Array{Int, 1}})
    precompile(Tuple{typeof(DecisionTrees.defineChartWith2Series), String, String, String, String, Int, Int, Int, Int})
    precompile(Tuple{typeof(DecisionTrees.get_sha1)})
    precompile(Tuple{typeof(DecisionTrees.my_sortperm_you_should_regularily_check_if_sortperm_in_base_has_become_more_efficient!), Array{Int, 1}, Base.SubArray{Float64, 1, Array{Float64, 1}, Tuple{Array{Int, 1}}, false}})
    precompile(Tuple{typeof(DecisionTrees.createSingleTreeExcel), Array{Int, 1}, Array{Int, 1}, DecisionTrees.ModelSettings, DecisionTrees.Tree, Array{DecisionTrees.Leaf, 1}, Array{Float64, 1}, Array{Int, 1}, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}})
    precompile(Tuple{typeof(DecisionTrees.dtm), DecisionTrees.DTMTable, DecisionTrees.ModelSettings})
    precompile(Tuple{typeof(DecisionTrees.create_leaves_array), DecisionTrees.Node{UInt8}})
    precompile(Tuple{typeof(DecisionTrees.fittedratio), DecisionTrees.Node{UInt8}})
end
