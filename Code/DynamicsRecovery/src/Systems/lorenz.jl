using DifferentialEquations

function lorenz!(du, u, p, t)
    du[1] = p[1] * (u[2] - u[1])
    du[2] = u[1] * (p[2] - u[3]) - u[2]
    du[3] = u[1] * u[2] - p[3] * u[3]
end # function lorenz!

function lorenz_data(u0, tspan, p, dt)
    system = ODEProblem(lorenz!, u0, tspan, p) # define ODE problem
    times = collect((tspan[1]):dt:(tspan[2])) # constant timestamps

    solve(system, DP5(), dt = dt, saveat = times) # solve with Explicit Runge-Kutta Method.
end # function lorenz_data

function lorenz_data(u0, tspan, p)
    # define ODE system
    system = ODEProblem(lorenz!, u0, tspan, p)

    solution = solve(system, Rosenbrock23()) # solve using Tsit5

    solution
end # function lorenz_data
