using Combinatorics

function construct_polynomials(degrees::Vector{T}, n::S) where T <: Integer where S <: Integer
    functions = []
    for deg in degrees
        # n = 0 is a special case 
        if deg == 0
            C(x...) = 1  # constant
            push!(functions,
                C
            )

            continue
        end # if 
        # if deg == 1 
        #     I(x...) = [x...] # identity
        #     push!(functions, I)

        #     continue
        # end # if

        
        # all nth degree polynomials will be all combinations with replacement of the variables 
        nth_degree_polynomial_iter = with_replacement_combinations(collect(1:n), deg)  

        for poly in nth_degree_polynomial_iter 
            P(x...) = prod(x[poly])

            push!(functions, P)
        end # for poly

    end # for n

    functions
end # function construct_polynomials
