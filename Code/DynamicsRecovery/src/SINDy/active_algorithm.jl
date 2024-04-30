using DifferentialEquations
using LinearAlgebra
using ProgressMeter
using StatsBase

using DynamicsRecovery.Systems

function automatic_lambda_fitting(X::Matrix{Float64}, Xdot::Matrix{Float64}, X_test::Matrix{Float64}, Xdot_test::Matrix{Float64}, Θx; kmax = 10, expected_multiple = 2.5, epsilon = 1e-2, proportional_λ_steps = 10, stopping_lambda = 10., polynomials = [], functions = [], functions1 = [])
    # initial lambda, Xi, and error
    λ = 0. # start at 0
    Ξ= sparse_regression(Θx, Xdot, λ, kmax)
    ε = measure_accuracy(Ξ, X_test, Xdot_test; polynomials = polynomials, functions = functions, functions1 = functions1)
    previous_ε = ε 

    previous_λ = 0.
    # loop until we see an expected jump
    while (previous_ε * expected_multiple > ε || previous_λ == 0 || λ == 0)
        @info "λ = " * string(λ) * ", ε = " * string(ε)
        previous_λ = λ # save previous value 

        λ = minimum(abs.(Ξ[Ξ .!= 0])) + epsilon # to guarantee remove a coefficient, take the minimum value in Xi plus some tiny epsilon
        Ξ = sparse_regression(Θx, Xdot, λ, kmax)

        previous_ε = ε
        ε = measure_accuracy(Ξ, X_test, Xdot_test; polynomials = polynomials, functions = functions, functions1 = functions1)

        # failsafe
        if λ > stopping_lambda 
            return (λ, Ξ)
        end # if
    end # while

    @info "Found jump in error at λ = " * string(λ)

    initial_guess_min = previous_λ # save in case
    initial_guess_max = λ # save in case


    # more fine adjustment of λ
    λ_step = (λ - previous_λ)/proportional_λ_steps
    λ = previous_λ
    Ξ = sparse_regression(Θx, Xdot, λ, kmax)
    previous_Ξ = Ξ
    initial_guess_min_Ξ = Ξ # save in case
    ε = measure_accuracy(Ξ, X_test, Xdot_test; polynomials = polynomials, functions = functions, functions1 = functions1)
    while (previous_ε * expected_multiple > ε)
        previous_λ = λ 
        λ += λ_step

        previous_Ξ = Ξ
        Ξ = sparse_regression(Θx, Xdot, λ, kmax)

        previous_ε = ε 
        ε = measure_accuracy(Ξ, X_test, Xdot_test; polynomials = polynomials, functions = functions, functions1 = functions1)

        # failsafe 
        if λ > initial_guess_max 
            return (initial_guess_min, initial_guess_min_Ξ) # this is our best guess 
        end # if 
        # @info "λ = " string(λ) * ", ε = " * string(ε)
    end # while

    return (previous_λ, previous_Ξ) # we want before the jump
end # function automatic_lambda_fitting

function pareto_sparsity_accuracy(λ0, X::Matrix{Float64}, Xdot::Matrix{Float64}, X_test::Matrix{Float64}, Xdot_test::Matrix{Float64}, Θx; kmax = 10, λ_step_divider = 1000, polynomials = [], functions = [], functions1 = [])
    Ξ_initial = sparse_regression(Θx, Xdot, λ0, kmax)
    while maximum(abs.(Ξ_initial)) > 0.
        λ0 = maximum(abs.(Ξ_initial)) + λ0 # this rapidly grows, but just searching for a zero solution rapidly
        Ξ_initial = sparse_regression(Θx, Xdot, λ0, kmax)
    end # while 

    @info "Found zero solution with λ = " * string(λ0)
    
    λ_step = λ0 / λ_step_divider
    λ = λ0
    Ξ = Ξ_initial

    λ_history = []
    sparsity_history = []
    accuracy_history = []
    
    while λ > 0. 
        Ξ = sparse_regression(Θx, Xdot, λ, kmax)

        push!(λ_history, λ)
        push!(sparsity_history, measure_sparsity(Ξ))
        push!(accuracy_history, measure_accuracy(Ξ, X_test, Xdot_test; polynomials = polynomials, functions = functions, functions1 = functions1))

        λ -= λ_step
    end

    (λ_history, sparsity_history, accuracy_history)
end # function pareto_sparsity_accuracy

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

    total_difference = BigFloat.(zeros(size(testXdot)[2]))

    total_derivative = BigFloat.(zeros(size(testXdot)[2]))

    for (x, xdot) ∈ zip(eachrow(testX), eachrow(testXdot))
        du = zeros(size(xdot))

        system!(du, x, (0.), 0.) # zeros are nothing

        total_difference += (xdot - du).^2
        total_derivative += xdot .^2
    end # for

    total_derivative = map(x -> x == 0. ? 1 : x, total_derivative) # edge case when derivatives are all 0

    average_diff = (total_difference) ./ total_derivative

    norm(average_diff)
    # (sqrt(total_difference)/size(testXdot)[1]) / norm(average)
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
