using LinearAlgebra

function sparse_regression(Θ::Matrix{Float64}, Xdot::Matrix{Float64}, λ::Float64, kmax::Integer) 
    # create initial guess for Xi as just regression
    Ξ = Θ\Xdot

    for _ in 1:kmax 
        # find coefficients that are less than λ 
        small_coefficients = abs.(Ξ) .< λ

        # set small coefficients to 0 
        Ξ[small_coefficients] .= 0

        # now individually run regression on non-zero coefficients
        for n in 1:(size(Xdot)[2]) # loop through all state dimensions 
            nonsmall = .!small_coefficients[:, n]
            
            Ξ[nonsmall, n] = Θ[:, nonsmall]\Xdot[:, n]
        end # for n
    end # for k 

    Ξ
end # function sparse_regression