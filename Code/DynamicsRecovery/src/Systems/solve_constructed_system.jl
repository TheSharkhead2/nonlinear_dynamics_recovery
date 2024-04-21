using DifferentialEquations

using DynamicsRecovery.Utils

function solve_constructed_system(Ξ, u0, tspan; polynomials = [], functions = [], functions1 = [])

    found_system! = get_constructed_system(Ξ; polynomials = polynomials, functions = functions, functions1 = functions1)

    system = ODEProblem(found_system!, u0, tspan)
    # solve(system, Rosenbrock23())
end # function solve_constructed_system

function get_constructed_system(Ξ; polynomials = [], functions = [], functions1 = [])

    # add polynomials terms
    append!(functions, construct_polynomials(polynomials, size(Ξ)[2]))

    # add functions in functions1
    for func in functions1
        for n in 1:size(Ξ)[2]
            # call function on the nth variable
            f(x...) = func(x[n])

            push!(functions, f)
        end # for n
    end # for func

    # strip out only relevant functions
    relevant_functions::Vector{Vector{Any}} = []
    relevant_functions_coefficients::Vector{Vector{Any}} = []

    for n in 1:size(Ξ)[2]
        n_rel_funcs = []
        n_rel_funcs_coef = []
        
        for i in 1:length(functions)
            # skip if the coefficient is 0
            if Ξ[i, n] == 0
                continue
            end # if 

            push!(n_rel_funcs, functions[i])
            push!(n_rel_funcs_coef, Ξ[i, n])
        end # for i

        # add zero function if no functions where added
        if isempty(n_rel_funcs) 
            z0(x...) = 0.0 # zero map
            push!(n_rel_funcs, z0)
            push!(n_rel_funcs_coef, 0.)
        end # if


        push!(relevant_functions, n_rel_funcs)
        push!(relevant_functions_coefficients, n_rel_funcs_coef)
    end # for n

    dimension = size(Ξ)[2]

    function found_system!(du, u, p, t)
        evaluated_functions::Vector{Vector{Any}} = []

        for n in 1:dimension 
            eval_funcs_n = []

            # evaulate each function
            for i in eachindex(relevant_functions[n])
                push!(eval_funcs_n, relevant_functions[n][i](u...))
            end # for i 

            push!(evaluated_functions, eval_funcs_n)
        end # for n

        for n in 1:dimension
            du[n] = sum(evaluated_functions[n] .* relevant_functions_coefficients[n]) # multiply all functions by the coefficients and sum up all the functions
        end # for n 
    end # function found_system!

    found_system!
end # function get_constructed_system
