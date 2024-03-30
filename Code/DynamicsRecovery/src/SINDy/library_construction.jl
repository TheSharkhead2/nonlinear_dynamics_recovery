using LinearAlgebra 

"""
    Θ(X; polynomials, functions) -> Matrix{Float64}

Constructs library of functions ``\\boldsymbol \\Theta (\\boldsymbol X)`` with desired function basis. `polynomials` is a vector of integers where each adds the corresponding nonlinearities of that degree. `functions` is a vector of functions that map from ``\\mathbb R^n`` to ``\\mathbb R``, representing a potential nonlinearity in the system. `functions1` is a vector of functions the map from ``\\mathbb R`` to ``\\mathbb R`` to simplify cases where you want nonlinearities in one variable applied to all variables.

"""
function Θ(X::Matrix{Float64}; 
    polynomials::Vector{Integer} = [], 
    functions::Vector{Function} = [],
    functions1::Vector{Function} = []
)
    if length(polynomials) == 0 && length(functions) == 0 && length(functions1) == 0 
        throw(ArgumentError("At least one candidate function needs to be supplied"))
    end # if

    # system dimension is important 
    dimension = size(X)[2]

    # vector to hold all candidate functions 
    function_library::Vector{Matrix{Float64}} = []

    # construct polynomial parts of the library 
    for n in polynomials 
        # n = 0 and n = 1 are special cases 
        if n == 0
            push!(function_library, ones(size(X))) # constants
        end # if 
        if n == 1 
            push!(function_library, X) # linear occurances 
        end # if
        
        

    end # for n
    
end # function Θ

