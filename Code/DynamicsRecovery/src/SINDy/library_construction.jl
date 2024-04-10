using LinearAlgebra 
using Combinatorics 

using DynamicsRecovery.Utils

"""
    Θ(X; polynomials, functions) -> Matrix{Float64}

Constructs library of functions ``\\boldsymbol \\Theta (\\boldsymbol X)`` with desired function basis. `polynomials` is a vector of integers where each adds the corresponding nonlinearities of that degree. `functions` is a vector of functions that map from ``\\mathbb R^n`` to ``\\mathbb R``, representing a potential nonlinearity in the system. `functions1` is a vector of functions the map from ``\\mathbb R`` to ``\\mathbb R`` to simplify cases where you want nonlinearities in one variable applied to all variables.

"""
function Θ(X::Matrix{Float64}; 
    polynomials::Vector{T} = [], 
    functions = [],
    functions1 = []
) where T <: Integer
    if length(polynomials) == 0 && length(functions) == 0 && length(functions1) == 0 
        throw(ArgumentError("At least one candidate function needs to be supplied"))
    end # if

    # vector to hold all candidate functions 
    function_library::Vector{Matrix{Float64}} = []

    polynomial_functions = construct_polynomials(polynomials, size(X)[2])

    append!(functions, polynomial_functions)

    for func in functions 
        push!(function_library,
            hcat([func(X[i, :]...) for i in 1:size(X)[1]]) # hcat is to convert to matrix 
        )
    end # for func

    # for functions1, just apply the function to all items in X 
    for func1 in functions1
        push!(function_library, func1.(X))
    end # for func1

    hcat(function_library...)
end # function Θ

