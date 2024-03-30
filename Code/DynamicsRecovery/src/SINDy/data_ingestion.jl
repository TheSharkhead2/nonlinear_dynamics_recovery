using NoiseRobustDifferentiation

const TVDIFFPARAMS_ITERATIONS::Integer = 100
const TVDIFFPARAMS_ALPHA::Real = 1e-1
const TVDIFFPARAMS_SCALE::String = "small"
const TVDIFFPARAMS_EPSILON::Real = 1e-6
const TVDIFFPARAMS_PRECOND = "none"

"""
Struct holding parameters for total variation regularization

"""
struct TvDiffParams 
    iterations::Integer 
    α::Real 
    scale::String 
    ε::Real
    precond::String
end # struct TvDiffParams

"""
    TvDiffParams() -> TvDiffParams

Default initialization for `TvDiffParams`.

"""
TvDiffParams() = TvDiffParams(TVDIFFPARAMS_ITERATIONS, TVDIFFPARAMS_ALPHA, TVDIFFPARAMS_SCALE, TVDIFFPARAMS_EPSILON, TVDIFFPARAMS_PRECOND)

"""
    TvDiffParams(; iterations = 100, α = 1e-1, scale = "small", ε = 1e-6, precond = "none") -> TvDiffParams

Keyword argument initialization for `TvDiffParams`.

"""
TvDiffParams(; iterations::Integer = TVDIFFPARAMS_ITERATIONS, α::Real = TVDIFFPARAMS_ALPHA, scale::String = TVDIFFPARAMS_SCALE, ε::Real = TVDIFFPARAMS_EPSILON, precond::String = TVDIFFPARAMS_PRECOND) = TvDiffParams(iterations, α, scale, ε, precond)

"""
    constructX(data, derivatives) -> Tuple{Matrix{Float64}, Matrix{Float64}} 

Constructs data matrices ``X`` and ``\\dot X`` from given data. Note that 
```math
X = \\begin{bmatrix} x_1(t_1) & x_2(t_1) & \\cdots & x_n(t_1) \\\\
                     x_1(t_2) & x_2(t_2) & \\cdots & x_n(t_2) \\\\
                     \\vdots  & \\vdots  & \\ddots & \\vdots  \\\\
                     x_1(t_m) & x_2(t_m) & \\cdots & x_n(t_m) 
    \\end{bmatrix}
```

"""
function constructX(data::Vector{Vector{Float64}}, derivatives::Vector{Vector{Float64}})
    X::Matrix{Float64} = vcat(transpose.(data)...)
    Xdot::Matrix{Float64} = vcat(transpose.(derivatives)...)

    (X, Xdot)
end # function constructX

"""
    constructX(data, tvdiff_params) -> Tuple{Matrix{Float64}, Matrix{Float64}}  

Constructs data matries ``X`` and ``\\dot X`` from given data by numerically approximating ``\\dot X``.

"""
function constructX(data::Vector{Vector{Float64}}, tvdiff_params::TvDiffParams)
    # calculate derivatives using total variation regularization 
    derivatives = tvdiff(data, tvdiff_params.iterations, tvdiff_params.α, ε = tvdiff_params.ε, scale = tvdiff_params.scale, precond = tvdiff_params.precond)

    # return data and derivatives
    constructX(data, derivatives)
end # function constructX

"""
    constructX(data) -> Tuple{Matrix{Float64}, Matrix{Float64}}  

Construct data matrices ``X`` and ``\\dot X`` from given data by numerically approximating ``\\dot X``. Uses fully default parameters for total variation regularization.

"""
function constructX(data::Vector{Vector{Float64}}) 
    tvdiff_params = TvDiffParams()

    constructX(data, tvdiff_params)
end # function constructX
