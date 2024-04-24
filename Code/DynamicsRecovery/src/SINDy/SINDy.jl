module SINDy 

export constructX
export TvDiffParams

include("data_ingestion.jl")

export Θ 

include("library_construction.jl")

export sparse_regression

include("sparse_regression.jl")

export active_SINDy, measure_accuracy, measure_sparsity, get_testing_X

include("active_algorithm.jl")

end # module SINDy

