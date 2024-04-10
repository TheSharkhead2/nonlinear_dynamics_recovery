using Combinatorics

function print_constructed_system(Ξ; polynomials = [], functions = [], functions1 = [])
    n = size(Ξ)[2] # system dimension

    polynomial_terms = []

    for deg in polynomials
        # n = 0 is a special case 
        if deg == 0
            C(x...) = 1  # constant
            push!(functions,
                C
            )

            continue
        end # if 
        # all nth degree polynomials will be all combinations with replacement of the variables 
        nth_degree_polynomial_iter = with_replacement_combinations(collect(1:n), deg)  

        # add to collection of polynomials
        append!(polynomial_terms, collect(nth_degree_polynomial_iter))
    end # for deg

    # vector of equation strings 
    equations = []

    for k in 1:n
        equation_str = "\\dot u_" * string(k) * " = " 

        # go through functions first
        for i in eachindex(functions)
            # if you have a nonzero coefficient, add it 
            if Ξ[i, k] != 0
                equation_str *= string(Ξ[i, k]) * "f_$i + "
            end # if 
        end # for i 

        # go through polynomials
        for i in eachindex(polynomial_terms)
            poly_i = i + length(functions)

            if Ξ[poly_i, k] != 0 
                poly_term = prod(["u_$j" for j in polynomial_terms[i]])

                equation_str *= string(Ξ[poly_i, k]) * poly_term * " + " 
            end # if 
        end # for i

        # go through remaining 1-parameter functions 
        for i in eachindex(functions1)
            func_i = length(functions) + length(polynomial_terms) + (i-1)*n

            for j in 1:n
                if Ξ[func_i + j, k] != 0
                    equation_str *= string(Ξ[func_i + j, k]) * "g_$i(u_$j) + "
                end # if
            end # for j
        end # for i

        push!(equations, equation_str)
    end # for k
    
    equations
end # print_constructed_system
