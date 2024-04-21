using DifferentialEquations
using LinearAlgebra
using ProgressMeter
using StatsBase

using DynamicsRecovery.Systems

# function λ_sparsity_fitting(λ0, X::Matrix{Float64}, Xdot::Matrix{Float64}, X_test::Matrix{Float64}, Xdot_test::Matrix{{Float64}})
# end 

"""
    active_SINDy(system, dt, tspan0, u0, p, sol_method, η, control_param_range, tvdiff_params, λ0, kmax, desired_sparsity, λ_stepping; polynomials = [], functions = [], functions1 = [])    

Assumes that last(u) is the control parameter and that the intital condition for that control parameter is in u0.
"""
function active_SINDy(system!, dt, tspan0, u0, p, sol_method, η, control_param_range, tvdiff_params, λ0, kmax, desired_sparsity, λ_stepping; polynomials = [], functions = [], functions1 = [])
    # first generate initial data from initial position
    problem_initial = ODEProblem(system!, u0, tspan0, p)
    saveat_dt_initial = collect((tspan0[1]):dt:(tspan0[2])) # times to save data at 
    sol_initial = solve(problem_initial, sol_method, dt = dt, saveat = saveat_dt_initial) # solve system

    sol_with_noise = add_noise_to_sample(sol_initial.u, η)

    # solve initial SINDy system
    X_initial, Xdot_initial = constructX(sol_with_noise, tvdiff_params, sol_initial.t)

    Xtest, Xdottest, X, Xdot = get_testing_X(X_initial, Xdot_initial, 0.3)

    # TODO: pick smaller basis for initial solve ?

    Θ_initial = Θ(X; polynomials = polynomials, functions = functions, functions1 = functions1)
    Ξ_initial = sparse_regression(Θ_initial, Xdot, λ0, kmax)

    λ = λ0
    Ξ = Ξ_initial

    progress_meter = ProgressThresh(desired_sparsity; desc = "Minimizing:")
    # want to search for a reasonable sparisty parameter
    for _ in 1:1000 # large loop, change to parameter later? 
        Ξ_new = sparse_regression(Θ_initial, Xdot, λ + λ_stepping, kmax)

        # don't want to go under threshold as that could caused problems (theoretically?)
        if measure_sparsity(Ξ_new) > desired_sparsity
            break 
        end # if

        update!(progress_meter, measure_sparsity(Ξ_new)) # progress bar

        λ += λ_stepping
        Ξ = Ξ_new
        println(
            λ
        )
    end # for 

    println(measure_accuracy(Ξ, Xtest, Xdottest; polynomials = polynomials))

    Ξ
end # function active_SINDy

function measure_sparsity(Ξ)
    # simple measurement of "sparsity"
    return sum(Ξ .== 0)/length(Ξ)
end # function measure_sparsity

function measure_accuracy(Ξ, testX::Matrix{Float64}, testXdot::Matrix{Float64}; polynomials =[], functions = [], functions1 = [])
    system! = get_constructed_system(Ξ; polynomials = polynomials, functions = functions, functions1 = functions1)

    total_difference = 0.

    for (x, xdot) ∈ zip(eachrow(testX), eachrow(testXdot))
        du = zeros(size(xdot))

        system!(du, x, (0.), 0.) # zeros are nothing

        total_difference += norm(xdot - du)^2 
    end # for

    sqrt(total_difference)/size(testXdot)[1]
end # function measure_accuracy

function get_testing_X(X::Matrix{Float64}, Xdot::Matrix{Float64}, a::Float64)
    test_indicies = sample(1:size(X)[1], Int(floor(size(X)[1] * a)); replace = false, ordered = true)    

    X_rows = collect(eachrow(X))
    Xdot_rows = collect(eachrow(Xdot))

    Xtest = vcat(transpose.(
        X_rows[test_indicies]
    )...)

    Xdottest = vcat(transpose.(
        Xdot_rows[test_indicies]
    )...)

    Xtrain = vcat(transpose.(
        deleteat!(X_rows, test_indicies)
    )...)
    Xdottrain = vcat(transpose.(
        deleteat!(Xdot_rows, test_indicies)
    )...)

    (Xtrain, Xdottrain, Xtest, Xdottest)
end # function get_testing_X
