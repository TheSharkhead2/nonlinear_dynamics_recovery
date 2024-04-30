### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 71998b4b-1bfc-4c3d-8c7e-5741de700162
using Pkg

# ╔═╡ 0e6831be-f2e1-11ee-1ecc-1f97692f0c45
begin
# using Pkg
Pkg.activate("./DynamicsRecovery")
using Revise
using DynamicsRecovery
end

# ╔═╡ aa29b69a-5750-4a5f-8e8c-fdc49d9a45b6
begin
	using DifferentialEquations
	using LinearAlgebra
	using Statistics
	using LaTeXStrings
	using Combinatorics
end

# ╔═╡ 1511b58f-0ec9-4dc3-90ec-d3169cf9189e
begin
	# Pkg.add("PyPlot")
	# Pkg.add("PlotlyBase")
	
	using Plots
	gr()
	# pyplot()
end

# ╔═╡ ee6d98bc-97dc-48f4-92fe-79a69afd51ca
begin
using Colors

function expand_palette(bg, cs; kwargs...)
    colors = palette(cs).colors.colors
    c = convert.(RGBA, distinguishable_colors(20, vcat(bg, colors); kwargs...))[2:end]
    return palette(c)
end
end

# ╔═╡ e4a6e6f9-43c6-489e-9ae2-451ee1f9e0df
begin 
	using PlotlyKaleido
	PlotlyKaleido.start()
end

# ╔═╡ a0a26a0d-a7e9-4168-996b-51f44fd96dad
# Pkg.add("LaTeXStrings")
# Pkg.add("PlotlyKaleido")
# Pkg.instantiate()

# ╔═╡ f7fb3eaf-0473-4074-a99c-892655fe12b7
begin
#inspired by nucleus theme for Atom
const dark_palette = [
    colorant"#FE4365", # red
    colorant"#eca25c", # orange
    colorant"#3f9778", # green
    colorant"#005D7F" # blue
]
const dark_bg = colorant"#000000"

const custom_theme = PlotTheme(Dict([
    :bg => dark_bg,
    :bginside => colorant"#000000",
    :fg => colorant"#ADB2B7",
    :fgtext => colorant"#FFFFFF",
    :fgguide => colorant"#FFFFFF",
    :fglegend => colorant"#FFFFFF",
    :legendfontcolor => colorant"#FFFFFF",
    :legendtitlefontcolor => colorant"#FFFFFF",
    :titlefontcolor => colorant"#FFFFFF",
    :palette => expand_palette(dark_bg, dark_palette; lchoices = [57], cchoices = [100]),
    :colorgradient => :fire])
)

add_theme(:custom, custom_theme)
end

# ╔═╡ c2016647-3b10-4041-a584-920ab5c116d1
# theme(:custom)
theme(:default)

# ╔═╡ 0d836d4e-d9d2-4fa1-9983-e541fb58ab04
function sin_cos_function_library(k)
	funcs = []
	for i in 1:k
		sink(x) = sin(i*x)
		cosk(x) = cos(i*x)
		append!(funcs, [sink, cosk])
	end
	funcs
end

# ╔═╡ 61ccdc43-c12f-4b59-b746-316c82b342a5
md"""
# Lorenz

"""

# ╔═╡ f2f7d603-4836-4010-95c7-10def9254590
begin

	lorenz_system = lorenz_data([-8., 7., 27.], (0.,100.), (10., 28., 8/3), 0.001)
	# lorenz_dt = lorenz_system.t[2] - lorenz_system.t[1]
	lorenz_system.t
end

# ╔═╡ 4cbdb8af-68ca-421f-853c-5bd8c4f4b464
plot(lorenz_system, idxs=(1,2,3), legend = false, xlabel = L"x", ylabel = L"y", zlabel = L"z", camera=(40,30), color=colorant"#00acff") # , xlabel = "x", ylabel = "y", zlabel = "z", axis = ([], false)

# ╔═╡ 7970ebf4-794c-41d0-979b-7e2f83163fad
# Plots.savefig("./lorenz_attractor.pdf")

# ╔═╡ b2a4fb51-90a0-49af-9a3e-f0d21a54495d
begin
	noisy_lorenz001 = add_noise_to_sample(lorenz_system.u, 0.01)
	noisy_lorenz01 = add_noise_to_sample(lorenz_system.u, 0.1)
end

# ╔═╡ 963a2ee7-fb5b-4ce9-88e1-da5d3f2ddbbd
noisy_lorenz1 = add_noise_to_sample(lorenz_system.u, 1.0)

# ╔═╡ d3c346a0-8970-4e12-aaf8-5884a0ef435c
X_lorenz_001, Xdot_lorenz_001 = constructX(noisy_lorenz001, lorenz_system.t)


# ╔═╡ d74cdb2f-7bfe-452c-b96b-1b9bf33c29f2
X_lorenz_01, Xdot_lorenz_01 = constructX(noisy_lorenz01, lorenz_system.t)

# ╔═╡ 3413048a-a378-414f-bd4a-f21beaf31289
X_lorenz_1, Xdot_lorenz_1 = constructX(noisy_lorenz1, TvDiffParams(; iterations = 100, α = 1e1, scale = "small", ε = 1e2, perform_tvdiff = true), lorenz_system.t)

# ╔═╡ 69042648-7393-453d-bd59-ad2b345e1ada
X_lorenz_001_tv, Xdot_lorenz_001_tv = constructX(noisy_lorenz001, TvDiffParams(; iterations = 10, α = 2e-1, scale = "small", ε = 1e6, perform_tvdiff = true), lorenz_system.t)

# ╔═╡ c948a893-a13c-45c2-a66b-8eff8ab85d68
X_lorenz_01_tv, Xdot_lorenz_01_tv = constructX(noisy_lorenz01, TvDiffParams(; iterations = 10, α = 1, scale = "small", ε = 1e6, perform_tvdiff = true), lorenz_system.t)

# ╔═╡ 369c1bca-dbed-4d84-9860-10b08e5e825c
begin 
	X_lorenz_001_td, Xdot_lorenz_001_td, X_lorenz_001_td_test, Xdot_lorenz_001_td_test = get_testing_X(X_lorenz_001_tv, Xdot_lorenz_001_tv, 0.1)
end

# ╔═╡ 71864910-7252-4408-be07-9bdd95192926
begin 
	X_lorenz_01_td, Xdot_lorenz_01_td, X_lorenz_01_td_test, Xdot_lorenz_01_td_test = get_testing_X(X_lorenz_01_tv, Xdot_lorenz_01_tv, 0.1)
end

# ╔═╡ afc779f4-0b87-4de6-b6ad-bc688b20b318
begin 
	X_lorenz_1_td, Xdot_lorenz_1_td, X_lorenz_1_td_test, Xdot_lorenz_1_td_test = get_testing_X(X_lorenz_1, Xdot_lorenz_1, 0.1)
end

# ╔═╡ a347d212-183c-4343-8171-1bc3a2ea7ef6
begin 
	Θ_lorenz_001_td = Θ(X_lorenz_001_td; polynomials=collect(0:7), functions1 = [sin_cos_function_library(10)...])

	Ξ_lorenz_001_td = sparse_regression(Θ_lorenz_001_td, Xdot_lorenz_001_td, 0.25, 10)
end

# ╔═╡ 383b7116-0612-4626-a318-3dc72ff5e35c
begin 
	lorenz_01_td_polynomials = collect(0:5)
	lorenz_01_td_functions1 = [sin_cos_function_library(10)...]
	
	Θ_lorenz_01_td = Θ(X_lorenz_01_td; polynomials=lorenz_01_td_polynomials, functions1 = lorenz_01_td_functions1)

	Ξ_lorenz_01_td = sparse_regression(Θ_lorenz_01_td, Xdot_lorenz_01_td, 0.4, 10)
end

# ╔═╡ b5b84157-9a8e-4cf7-97f0-5e8f168d15b9
begin 
	lorenz_1_td_polynomials = collect(0:5)
	lorenz_1_td_functions1 = [] #[sin_cos_function_library(10)...]
	
	Θ_lorenz_1_td = Θ(X_lorenz_1_td; polynomials=lorenz_1_td_polynomials, functions1 = lorenz_1_td_functions1)

	Ξ_lorenz_1_td = sparse_regression(Θ_lorenz_1_td, Xdot_lorenz_1_td, 0.416, 10)
end

# ╔═╡ 940b3efb-bbe9-4130-879d-44b12917e209
print_constructed_system(Ξ_lorenz_001_td; polynomials = collect(0:7), functions1 = [sin_cos_function_library(10)...])

# ╔═╡ 5bbf1cdc-7568-4485-b81f-2d0a9754244b
print_constructed_system(Ξ_lorenz_01_td; polynomials=lorenz_01_td_polynomials, functions1 = lorenz_01_td_functions1)

# ╔═╡ 83863d85-1276-407e-aef5-6ef41e802fe2
print_constructed_system(Ξ_lorenz_1_td; polynomials=lorenz_1_td_polynomials, functions1 = lorenz_1_td_functions1)

# ╔═╡ 1c29a5e0-014d-4f3f-be6f-560dd1533863
begin
lorenz_001_system_td = solve_constructed_system(Ξ_lorenz_001_td, [-8., 7, 27], (0., 100.); polynomials=collect(0:7), functions1 = [sin_cos_function_library(10)...])
lorenz_001_td_sol = solve(lorenz_001_system_td, Tsit5())
end

# ╔═╡ cb9c98c7-e57b-424d-980b-24013fbbeae9
begin
lorenz_01_system_td = solve_constructed_system(Ξ_lorenz_01_td, [-8., 7, 27], (0., 100.); polynomials=lorenz_01_td_polynomials, functions1 = lorenz_01_td_functions1)
lorenz_01_td_sol = solve(lorenz_01_system_td, Tsit5())
end

# ╔═╡ ca03bce6-6477-4235-83f5-9bfdbafa1739
begin
lorenz_1_system_td = solve_constructed_system(Ξ_lorenz_1_td, [-8., 7, 27], (0., 100.); polynomials=lorenz_1_td_polynomials, functions1 = lorenz_1_td_functions1)
lorenz_1_td_sol = solve(lorenz_1_system_td, Tsit5())
end

# ╔═╡ 84671725-6117-4ea1-80e7-2ebf578f7b79
plot(lorenz_001_td_sol, idxs=(1,2,3), legend = false, xlabel = L"x", ylabel = L"y", zlabel = L"z", camera=(40,30), color=colorant"#00acff") # , axis = ([], false)

# ╔═╡ 47e7f92a-ced3-4632-81c8-1309ef873637
# Plots.savefig("./found_lorenz_001.pdf")

# ╔═╡ 8e5c1398-c6ef-4b68-991d-d3045d3be84e
plot(lorenz_01_td_sol, idxs=(1,2,3), legend = false, xlabel = L"x", ylabel = L"y", zlabel = L"z", camera=(40,30), color=colorant"#00acff") # , axis = ([], false)

# ╔═╡ 6206f493-513c-4be8-bfee-6c8ca16609dc
# Plots.savefig("./found_lorenz_01.pdf")

# ╔═╡ 1c18d80f-0e77-4676-9302-42d5cb396184
plot(lorenz_1_td_sol, idxs=(1,2,3), legend = false, xlabel = L"x", ylabel = L"y", zlabel = L"z", camera=(40,30), color=colorant"#00acff") # , axis = ([], false)

# ╔═╡ 39d59c5c-38d0-49d9-ae6c-079ac9f81901
# Plots.savefig("./found_lorenz_1.pdf")

# ╔═╡ 3fc03792-830e-42f9-b4b7-fbd155cc4142
measure_accuracy(Ξ_lorenz_001_td, X_lorenz_001_td_test, Xdot_lorenz_001_td_test; polynomials = collect(0:7), functions1 = [sin_cos_function_library(10)...])

# ╔═╡ 5eb64ceb-aeac-4ea9-a6c6-ddf49937a860
measure_accuracy(Ξ_lorenz_01_td, X_lorenz_01_td_test, Xdot_lorenz_01_td_test; polynomials=lorenz_01_td_polynomials, functions1 = lorenz_01_td_functions1)

# ╔═╡ 25166eec-a00d-4c5d-8ac6-f70208706613
measure_accuracy(Ξ_lorenz_1_td, X_lorenz_1_td_test, Xdot_lorenz_1_td_test; polynomials=lorenz_1_td_polynomials, functions1 = lorenz_1_td_functions1)

# ╔═╡ 529c0dc4-b9e5-40d9-bcc0-872792298fc0
# lorenz_λ_history, lorenz_sparsity_history, lorenz_accuracy_history = pareto_sparsity_accuracy(5., X_lorenz_001_td, Xdot_lorenz_001_td, X_lorenz_001_td_test, Xdot_lorenz_001_td_test, Θ_lorenz_001_td; polynomials=collect(0:7), functions1 = [sin_cos_function_library(10)...])

# ╔═╡ 1a5465b5-bd3c-4c9f-9921-983f99492ed5
# offset_lorenz_sparsity = 1. .- lorenz_sparsity_history

# ╔═╡ 7d44b26e-b916-4145-8147-8c24df66b5cd
# negative_lorenz_accuracy = -1 .* lorenz_accuracy_history

# ╔═╡ b6350c06-cfb1-4ea4-bea6-1d95b1d3c5fd
# plot(offset_lorenz_sparsity, negative_lorenz_accuracy, xlabel = "Measured Sparsity", ylabel = L"-\varepsilon", legend = false)

# ╔═╡ d9050f03-98d6-477e-b18f-07034bf177b3
# Plots.savefig("./large_scale_error_sparsity.pdf")

# ╔═╡ b5a165c5-6371-4354-be1c-eedf8d322de4
# plot(lorenz_λ_history, lorenz_sparsity_history)

# ╔═╡ 3e02c732-d1b6-4768-86a8-fb244372c1bc
# plot(lorenz_λ_history, lorenz_accuracy_history)

# ╔═╡ 590df629-ec43-4a0d-8081-21a407e16df3

# plot(lorenz_λ_history[938:length(lorenz_λ_history)], lorenz_accuracy_history[938:length(lorenz_λ_history)])


# ╔═╡ c027b7db-7d2d-4fc7-89b9-ffb5d4501a75
# plot(lorenz_λ_history[938:length(lorenz_λ_history)], lorenz_sparsity_history[938:length(lorenz_λ_history)])

# ╔═╡ 0d0cd6f9-a0c7-49f1-b4c4-c8d76ffaa2a3
# plot(offset_lorenz_sparsity[938:length(lorenz_λ_history)], negative_lorenz_accuracy[938:length(lorenz_λ_history)], legend = false, xlabel = "Measured Sparsity", ylabel = L"-\varepsilon")

# ╔═╡ 020e0c54-955b-4117-9b3e-251cf3d7b7d5
# Plots.savefig("./small_scale_error_sparsity.pdf")

# ╔═╡ 5fa8c07b-5bd3-4136-95b2-09ea2eedbab1
begin 
	Ξ_lorenz_001_lambda_1 = sparse_regression(Θ_lorenz_001_td, Xdot_lorenz_001_td, 0.943433, 10)
	Ξ_lorenz_001_lambda_2 = sparse_regression(Θ_lorenz_001_td, Xdot_lorenz_001_td, 0.67388, 10)
	Ξ_lorenz_001_lambda_3 = sparse_regression(Θ_lorenz_001_td, Xdot_lorenz_001_td, 0.0599005, 10)
	Ξ_lorenz_001_lambda_4 = sparse_regression(Θ_lorenz_001_td, Xdot_lorenz_001_td, 0.0299502, 10)
end

# ╔═╡ 6246c8a2-638e-4eb6-88a9-2a335f86cd88
begin
	lorenz_001_system_lambda_1 = solve_constructed_system(Ξ_lorenz_001_lambda_1, [-8., 7, 27], (0., 100.); polynomials=collect(0:7), functions1 = [sin_cos_function_library(10)...])
	lorenz_001_td_sol_lambda_1 = solve(lorenz_001_system_lambda_1, Tsit5())

	lorenz_001_system_lambda_2 = solve_constructed_system(Ξ_lorenz_001_lambda_2, [-8., 7, 27], (0., 100.); polynomials=collect(0:7), functions1 = [sin_cos_function_library(10)...])
	lorenz_001_td_sol_lambda_2 = solve(lorenz_001_system_lambda_2, Tsit5())

	lorenz_001_system_lambda_3 = solve_constructed_system(Ξ_lorenz_001_lambda_3, [-8., 7, 27], (0., 100.); polynomials=collect(0:7), functions1 = [sin_cos_function_library(10)...])
	lorenz_001_td_sol_lambda_3 = solve(lorenz_001_system_lambda_3, Tsit5())

	lorenz_001_system_lambda_4 = solve_constructed_system(Ξ_lorenz_001_lambda_4, [-8., 7, 27], (0., 100.); polynomials=collect(0:7), functions1 = [sin_cos_function_library(10)...])
	lorenz_001_td_sol_lambda_4 = solve(lorenz_001_system_lambda_4, Tsit5())
end

# ╔═╡ d1a0d67e-a6a5-400b-b698-9257b56c6b99
print_constructed_system(Ξ_lorenz_001_lambda_1; polynomials=collect(0:7), functions1 = [sin_cos_function_library(10)...])

# ╔═╡ e61f7ee9-2f25-45e2-9952-7abc8123537e
print_constructed_system(Ξ_lorenz_001_lambda_2; polynomials=collect(0:7), functions1 = [sin_cos_function_library(10)...])

# ╔═╡ e70e3f6e-437e-4ee7-a4a8-76a72791fe37
print_constructed_system(Ξ_lorenz_001_lambda_3; polynomials=collect(0:7), functions1 = [sin_cos_function_library(10)...])

# ╔═╡ b2c01934-98a3-4daa-b957-ede0a6c3e20f
print_constructed_system(Ξ_lorenz_001_lambda_4; polynomials=collect(0:7), functions1 = [sin_cos_function_library(10)...])

# ╔═╡ 9cd257db-474d-4be6-be81-cd847be71cef
measure_accuracy(Ξ_lorenz_001_lambda_1, X_lorenz_001_td_test, Xdot_lorenz_001_td_test; polynomials=collect(0:7), functions1 = [sin_cos_function_library(10)...])

# ╔═╡ 438a33b0-2d93-4bc9-aab7-b789d89eb010
measure_accuracy(Ξ_lorenz_001_lambda_2, X_lorenz_001_td_test, Xdot_lorenz_001_td_test; polynomials=collect(0:7), functions1 = [sin_cos_function_library(10)...])

# ╔═╡ 284b77f9-1f84-4f79-84c5-b2299f8ff674
measure_accuracy(Ξ_lorenz_001_lambda_3, X_lorenz_001_td_test, Xdot_lorenz_001_td_test; polynomials=collect(0:7), functions1 = [sin_cos_function_library(10)...])

# ╔═╡ c0dc5518-edc3-4564-9f00-0632691c5a9c
measure_accuracy(Ξ_lorenz_001_lambda_4, X_lorenz_001_td_test, Xdot_lorenz_001_td_test; polynomials=collect(0:7), functions1 = [sin_cos_function_library(10)...])

# ╔═╡ f10bb9b2-0c9e-4c9b-8057-003a91da4128
plot(lorenz_001_td_sol_lambda_1, idxs=(1,2,3), legend = false, xlabel = L"x", ylabel = L"y", zlabel = L"z", camera=(40,30), color=colorant"#00acff") # , axis = ([], false)

# ╔═╡ 8c207abe-db98-43e0-9df8-5b662fac803f
# Plots.savefig("./lorenz_lambda_1.pdf")

# ╔═╡ 7a70e0fe-728e-411d-9d7c-adb5d39efdb1
plot(lorenz_001_td_sol_lambda_2, idxs=(1,2,3), legend = false, xlabel = L"x", ylabel = L"y", zlabel = L"z", camera=(40,30), color=colorant"#00acff") # , axis = ([], false)

# ╔═╡ 23f965a8-29e6-42e3-9fa5-7d02b20a1cc3
# Plots.savefig("./lorenz_lambda_2.pdf")

# ╔═╡ eb9455c5-238e-44e9-ab8b-5872d9b60d9e
plot(lorenz_001_td_sol_lambda_3, idxs=(1,2,3), legend = false, xlabel = L"x", ylabel = L"y", zlabel = L"z", camera=(40,30), color=colorant"#00acff") # , axis = ([], false)

# ╔═╡ 2a6f9358-b196-418e-adf7-1d6348c92ee1
# Plots.savefig("./lorenz_lambda_3.pdf")

# ╔═╡ 27efd12e-5608-40c1-889d-a247562ece3a
plot(lorenz_001_td_sol_lambda_4, idxs=(1,2,3), legend = false, xlabel = L"x", ylabel = L"y", zlabel = L"z", camera=(40,30), color=colorant"#00acff") # , axis = ([], false)

# ╔═╡ d59d0489-bd49-4ac8-a8eb-2a95ed1579b4
# Plots.savefig("./lorenz_lambda_4.pdf")

# ╔═╡ 7b4cb49e-d263-4711-bb1d-b90aa2e1df91
# Θ_lorenz_001 = Θ(X_lorenz_001; polynomials=[0, 1,2,3,4,5,6,7], functions1 = [sin_cos_function_library(10)...])

# ╔═╡ 12ec1247-f93a-4f49-a847-61dc9bf5bb0f
# Θ_lorenz_01 = Θ(X_lorenz_01; polynomials = [0,1,2,3,4,5,6,7]) # , functions1 = [sin_cos_function_library(10)...]

# ╔═╡ e22bfbb2-bb20-441d-8b40-82d5d5760d89
# Θ_lorenz_001_tv = Θ(X_lorenz_001_tv; polynomials=[0, 1,2,3,4,5,6,7], functions1 = [sin_cos_function_library(10)...])

# ╔═╡ 1c029b31-5d6a-423f-90fc-06c902977962
# Θ_lorenz_01_tv = Θ(X_lorenz_01_tv; polynomials=[0, 1,2,3,4,5,6,7], functions1 = [sin_cos_function_library(10)...])

# ╔═╡ 03c887da-dd03-49f0-9bcd-5bc33a619a7a
# Ξ_lorenz_001 = sparse_regression(Θ_lorenz_001, Xdot_lorenz_001, 0.25, 10)

# ╔═╡ e90d0073-3530-4662-80be-236675c648c3
# Ξ_lorenz_01 = sparse_regression(Θ_lorenz_01, Xdot_lorenz_01, 0.025, 10)

# ╔═╡ d04e4a91-1feb-4012-bff0-51dcf9a04f70
# Ξ_lorenz_001_tv = sparse_regression(Θ_lorenz_001_tv, Xdot_lorenz_001_tv, 0.25, 10)

# ╔═╡ bd8867fe-0a1c-4e3d-910d-32acba5ef2af
# Ξ_lorenz_01_tv = sparse_regression(Θ_lorenz_01_tv, Xdot_lorenz_01_tv, 0.25, 10)

# ╔═╡ 2a2d55ac-2a06-4f02-8847-82ba264bc40c
# print_constructed_system(Ξ_lorenz_001; polynomials=[0, 1,2,3,4,5,6,7], functions1 = [sin_cos_function_library(10)...])

# ╔═╡ 6f6c0019-ed7d-4a31-aa06-7291cd6cd532
# print_constructed_system(Ξ_lorenz_001_tv; polynomials=[0, 1,2,3,4,5,6,7], functions1 = [sin_cos_function_library(10)...])

# ╔═╡ 1682c040-1735-48d9-93a5-1a1900a1d1cb
# print_constructed_system(Ξ_lorenz_01_tv; polynomials=[0, 1,2,3,4,5,6,7], functions1 = [sin_cos_function_library(10)...])

# ╔═╡ 8eb49379-7b2a-4edb-924f-b6417450febf
# lorenz_001_system = solve_constructed_system(Ξ_lorenz_001, [-8., 7, 27], (0., 100.); polynomials=[0, 1,2,3,4,5,6,7], functions1 = [sin_cos_function_library(10)...])

# ╔═╡ bfeee94b-b98d-4345-8b8a-efd9eda09c22
# lorenz_01_system = solve_constructed_system(Ξ_lorenz_01, [-8., 7, 27], (0., 100.); polynomials=[0, 1,2,3,4,5,6,7]) # , functions1 = [sin_cos_function_library(10)...]

# ╔═╡ e00de518-1c38-452c-8b79-4adbfeecea53
# lorenz_001_tv_system = solve_constructed_system(Ξ_lorenz_001_tv, [-8., 7, 27], (0., 100.); polynomials=[0, 1,2,3,4,5,6,7], functions1 = [sin_cos_function_library(10)...])

# ╔═╡ e4fe9390-6955-405c-86cc-cc9defb9ca39
# lorenz_01_tv_system = solve_constructed_system(Ξ_lorenz_01_tv, [-8., 7, 27], (0., 1000.); polynomials=[0, 1,2,3,4,5,6,7], functions1 = [sin_cos_function_library(10)...])

# ╔═╡ c1ea010b-9549-43d7-8aa2-9b6a92c00d82
# lorenz_001_solution = solve(lorenz_001_system, DP5())

# ╔═╡ 0814a304-c439-40e5-84a9-4e501aaf2326
# lorenz_01_solution = solve(lorenz_01_system, DP5())

# ╔═╡ 3d9a773b-39a9-41bd-b9a9-4c79a33992de
# lorenz_001_tv_solution = solve(lorenz_001_tv_system, DP5())

# ╔═╡ b4a92457-6529-4ab9-872f-0bc2ba099ff4
# lorenz_01_tv_solution = solve(lorenz_01_tv_system, DP5())

# ╔═╡ 7a8cc3aa-5b30-46b9-83db-3f47d0377c13
# plot(lorenz_001_solution, idxs=(1,2,3), xlabel="x", ylabel="y", zlabel="z", camera=(40,30), legend = false)

# ╔═╡ e3f29ead-9e94-48e9-9dfa-a07fba5a0fef
# plot(lorenz_001_solution, idxs=(1,2,3), legend = false, xlabel = "", ylabel = "", zlabel = "", camera=(40,30), axis = ([], false), color=colorant"#00acff")

# ╔═╡ 9ced946d-5f6c-4801-a204-ab2a319cb861
# Plots.savefig("./found_lorenz_001.pdf")

# ╔═╡ f35be9d2-efea-4f0a-869d-25bd24f946a0
# plot(lorenz_01_solution, idxs=(1,2,3))

# ╔═╡ e0003431-ffef-429e-80c8-ba1c7f1774d2
# plot(lorenz_001_tv_solution, idxs=(1,2,3), xlabel="x", ylabel="y", zlabel="z", camera=(40,30), legend = false)

# ╔═╡ d333670f-c9c6-461f-81d8-84f7bb974d85
# plot(lorenz_01_tv_solution, idxs=(1,2,3), xlabel="x", ylabel="y", zlabel="z", camera=(40,30), legend = false)

# ╔═╡ 768ff383-57d4-4203-8c25-4cba99c993bf
begin
	hopf_solutions = []
	# save data and derivatives from multiple starting points
	hopf_data = []
	hopf_derivatives = []

	hopf_tspan = (0., 75.)
	hopf_ω = 1
	hopf_A = 1

	hopf_η = 0.005 # noise

	hopf_dt = 0.0025

	hopf_dt_saveat = collect((hopf_tspan[1]):hopf_dt:(hopf_tspan[2]))

	# parameter values from SINDy paper
	for μ in [-0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55]
		u0 = [2.,0] + rand(2)
		hopf_system = ODEProblem(hopf_normal!, u0, hopf_tspan, (μ,  hopf_ω, hopf_A))
		hopf_solution = solve(hopf_system, Rosenbrock23(), dt = hopf_dt, saveat=hopf_dt_saveat)

		# add "sensor noise"
		hopf_noisy = add_noise_to_sample(hopf_solution.u, hopf_η)

		# compute derivatives using tvdiff
		X_hopf, Xdot_hopf = constructX(hopf_noisy, TvDiffParams(; iterations = 20, α = 2e1, ε = 1e-3, perform_tvdiff = true, diff_kernel = "abs"), hopf_solution.t)
		
		# save
		push!(hopf_data, hcat(X_hopf, fill(μ, size(X_hopf)[1])))
		push!(hopf_derivatives, hcat(Xdot_hopf, fill(0, size(X_hopf)[1])))
		push!(hopf_solutions, hopf_solution)

		# don't do inner initial condition for case where we don't have limit cycle
		if μ < 0 
			continue
		end # if 

		u02 = [0.01, 0] + 0.01 *rand(2)

		# do same thing again with new initial condition 
		hopf_system2 = ODEProblem(hopf_normal!, u02, hopf_tspan, (μ,  hopf_ω, hopf_A))
		hopf_solution2 = solve(hopf_system2, Rosenbrock23(), dt = hopf_dt, saveat=hopf_dt_saveat)

		# add "sensor noise"
		hopf_noisy2 = add_noise_to_sample(hopf_solution2.u, hopf_η)

		# compute derivatives using tvdiff
		X_hopf2, Xdot_hopf2 = constructX(hopf_noisy2, TvDiffParams(; iterations = 20, α = 2e1, ε = 1e-3, perform_tvdiff = true, diff_kernel = "abs"), hopf_solution2.t)

		# save
		push!(hopf_data, hcat(X_hopf2, fill(μ, size(X_hopf2)[1])))
		push!(hopf_derivatives, hcat(Xdot_hopf2, fill(0, size(X_hopf2)[1])))

		push!(hopf_solutions, hopf_solution2)
	end # for μ


	X_hopf_combined = vcat(hopf_data...)
	Xdot_hopf_combined = vcat(hopf_derivatives...)
end

# ╔═╡ 6168c7b1-6f44-4f80-ba98-257e9d447e98
function hopf_normal_mu!(du, u, p, t)
    # μ: p[1], ω: p[2], A: p[3]
    
    du[1] = u[3]*u[1] - p[2]*u[2] - p[3]*u[1]*(u[1]^2 + u[2]^2)
    du[2] = p[2]*u[1] + u[3]*u[2] - p[3]*u[2]*(u[1]^2 + u[2]^2)
	du[3] = 0.0
end # function hopf_normal_mu!


# ╔═╡ 37dc2ad0-d356-4f78-9841-80836d9f0692
X_hopf_train, Xdot_hopf_train, X_hopf_test, Xdot_hopf_test = get_testing_X(X_hopf_combined, Xdot_hopf_combined, 0.1)

# ╔═╡ 3781cba0-2f81-471d-accc-de4d8546c11c
begin 
	hopf_polynomials = collect(0:5)
	hopf_functions1 = [] # [sin_cos_function_library(10)...]
	Θ_hopf = Θ(X_hopf_train; polynomials=hopf_polynomials, functions1 = hopf_functions1)
end

# ╔═╡ 20d272ea-c701-4edf-a393-a48f514465fc
begin 
	Ξ_hopf = sparse_regression(Θ_hopf, Xdot_hopf_train, 0.6, 10)
end

# ╔═╡ 261776eb-e90e-4471-8276-cb1b934cde28
print_constructed_system(Ξ_hopf; polynomials = hopf_polynomials, functions1 = hopf_functions1)

# ╔═╡ f66ceb82-0780-4149-a50d-5f6c1ffe6e6b
measure_accuracy(Ξ_hopf, X_hopf_test, Xdot_hopf_test; polynomials = hopf_polynomials, functions1 = hopf_functions1)

# ╔═╡ f6cb782c-934b-4342-b008-10be1fda2021
begin
	found_hopf_mu_neg02_1_system = solve_constructed_system(Ξ_hopf, [2, 0., -0.2], (0.0, 100.); polynomials=hopf_polynomials, functions1 = hopf_functions1)
	found_hopf_mu_neg015_1_system = solve_constructed_system(Ξ_hopf, [2, 0., -0.15], (0.0, 100.); polynomials=hopf_polynomials, functions1 = hopf_functions1)
	found_hopf_mu_neg01_1_system = solve_constructed_system(Ξ_hopf, [2, 0., -0.1], (0.0, 100.); polynomials=hopf_polynomials, functions1 = hopf_functions1)
	found_hopf_mu_neg005_1_system = solve_constructed_system(Ξ_hopf, [2, 0., -0.05], (0.0, 100.); polynomials=hopf_polynomials, functions1 = hopf_functions1)
	found_hopf_mu_0_1_system = solve_constructed_system(Ξ_hopf, [2, 0., 0.], (0.0, 100.); polynomials=hopf_polynomials, functions1 = hopf_functions1)
	found_hopf_mu_005_1_system = solve_constructed_system(Ξ_hopf, [2, 0., 0.05], (0.0, 100.); polynomials=hopf_polynomials, functions1 = hopf_functions1)
	found_hopf_mu_005_2_system = solve_constructed_system(Ξ_hopf, [0.1, 0., 0.05], (0.0, 100.); polynomials=hopf_polynomials, functions1 = hopf_functions1)
	found_hopf_mu_01_1_system = solve_constructed_system(Ξ_hopf, [2, 0., 0.1], (0.0, 100.); polynomials=hopf_polynomials, functions1 = hopf_functions1)
	found_hopf_mu_01_2_system = solve_constructed_system(Ξ_hopf, [0.1, 0., 0.1], (0.0, 100.); polynomials=hopf_polynomials, functions1 = hopf_functions1)
	found_hopf_mu_015_1_system = solve_constructed_system(Ξ_hopf, [2, 0., 0.15], (0.0, 100.); polynomials=hopf_polynomials, functions1 = hopf_functions1)
	found_hopf_mu_015_2_system = solve_constructed_system(Ξ_hopf, [0.1, 0., 0.15], (0.0, 100.); polynomials=hopf_polynomials, functions1 = hopf_functions1)
	found_hopf_mu_02_1_system = solve_constructed_system(Ξ_hopf, [2, 0., 0.2], (0.0, 100.); polynomials=hopf_polynomials, functions1 = hopf_functions1)
	found_hopf_mu_02_2_system = solve_constructed_system(Ξ_hopf, [0.1, 0., 0.2], (0.0, 100.); polynomials=hopf_polynomials, functions1 = hopf_functions1)
	found_hopf_mu_025_1_system = solve_constructed_system(Ξ_hopf, [2, 0., 0.25], (0.0, 100.); polynomials=hopf_polynomials, functions1 = hopf_functions1)
	found_hopf_mu_025_2_system = solve_constructed_system(Ξ_hopf, [0.1, 0., 0.25], (0.0, 100.); polynomials=hopf_polynomials, functions1 = hopf_functions1)
	found_hopf_mu_03_1_system = solve_constructed_system(Ξ_hopf, [2, 0., 0.3], (0.0, 100.); polynomials=hopf_polynomials, functions1 = hopf_functions1)
	found_hopf_mu_03_2_system = solve_constructed_system(Ξ_hopf, [0.1, 0., 0.3], (0.0, 100.); polynomials=hopf_polynomials, functions1 = hopf_functions1)

	found_hopf_mu_neg02_1_solution = solve(found_hopf_mu_neg02_1_system, Tsit5())
	found_hopf_mu_neg015_1_solution = solve(found_hopf_mu_neg015_1_system, Tsit5())
	found_hopf_mu_neg01_1_solution = solve(found_hopf_mu_neg01_1_system, Tsit5())
	found_hopf_mu_neg005_1_solution = solve(found_hopf_mu_neg005_1_system, Tsit5())
	found_hopf_mu_0_1_solution = solve(found_hopf_mu_0_1_system, Tsit5())
	found_hopf_mu_005_1_solution = solve(found_hopf_mu_005_1_system, Tsit5())
	found_hopf_mu_005_2_solution = solve(found_hopf_mu_005_2_system, Tsit5())
	found_hopf_mu_01_1_solution = solve(found_hopf_mu_01_1_system, Tsit5())
	found_hopf_mu_01_2_solution = solve(found_hopf_mu_01_2_system, Tsit5())
	found_hopf_mu_015_1_solution = solve(found_hopf_mu_015_1_system, Tsit5())
	found_hopf_mu_015_2_solution = solve(found_hopf_mu_015_2_system, Tsit5())
	found_hopf_mu_02_1_solution = solve(found_hopf_mu_02_1_system, Tsit5())
	found_hopf_mu_02_2_solution = solve(found_hopf_mu_02_2_system, Tsit5())
	found_hopf_mu_025_1_solution = solve(found_hopf_mu_025_1_system, Tsit5())
	found_hopf_mu_025_2_solution = solve(found_hopf_mu_025_2_system, Tsit5())
	found_hopf_mu_03_1_solution = solve(found_hopf_mu_03_1_system, Tsit5())
	found_hopf_mu_03_2_solution = solve(found_hopf_mu_03_2_system, Tsit5())
end

# ╔═╡ 379bba48-ebc0-415d-bd54-25707b97bab1
begin
	plot(found_hopf_mu_neg02_1_solution, idxs=(3, 1,2), color = :blue)
	plot!(found_hopf_mu_neg015_1_solution, idxs=(3, 1,2), color = :blue)
	plot!(found_hopf_mu_neg01_1_solution, idxs=(3, 1,2), color = :blue)
	plot!(found_hopf_mu_neg005_1_solution, idxs=(3, 1,2), color = :blue)
	plot!(found_hopf_mu_0_1_solution, idxs=(3, 1,2), color = :blue)
	plot!(found_hopf_mu_005_1_solution, idxs=(3, 1,2), color = :blue)
	plot!(found_hopf_mu_005_2_solution, idxs=(3, 1,2), color = :red)
	plot!(found_hopf_mu_01_1_solution, idxs=(3, 1,2), color = :blue)
	plot!(found_hopf_mu_01_2_solution, idxs=(3, 1,2), color = :red)
	plot!(found_hopf_mu_015_1_solution, idxs=(3, 1,2), color = :blue)
	plot!(found_hopf_mu_015_2_solution, idxs=(3, 1,2), color = :red)
	plot!(found_hopf_mu_02_1_solution, idxs=(3, 1,2), color = :blue)
	plot!(found_hopf_mu_02_2_solution, idxs=(3, 1,2), color = :red, zlabel = L"x", ylabel = L"y", xlabel = L"\mu", camera = (30,10), legend = false)
	# plot!(found_hopf_mu_025_1_solution, idxs=(3, 1,2), color = :blue)
	# plot!(found_hopf_mu_025_2_solution, idxs=(3, 1,2), color = :red, zlabel = L"x", ylabel = L"y", xlabel = L"\mu", camera = (30,10), legend = false)
	# plot!(found_hopf_mu_03_1_solution, idxs=(3, 1,2), color = :blue)
	# plot!(found_hopf_mu_03_2_solution, idxs=(3, 1,2), color = :red, zlabel = L"x", ylabel = L"y", xlabel = L"\mu", camera = (30,10), legend = false)
end

# ╔═╡ 51281685-dea0-4c19-b6c6-893b265b051c
# Plots.savefig("./found_hopf_bifurcation_0005_2.pdf")

# ╔═╡ ff4ae999-17c3-4487-abb8-28204aa7b73f
begin
	hopf_mu_neg02_1_system = ODEProblem(hopf_normal_mu!, [2., 0., -0.2], hopf_tspan, (-0.2,  hopf_ω, hopf_A))
	hopf_mu_neg015_1_system = ODEProblem(hopf_normal_mu!, [2., 0., -0.15], hopf_tspan, (-0.15,  hopf_ω, hopf_A))
	hopf_mu_neg01_1_system = ODEProblem(hopf_normal_mu!, [2., 0., -0.1], hopf_tspan, (-0.1,  hopf_ω, hopf_A))
	hopf_mu_neg005_1_system = ODEProblem(hopf_normal_mu!, [2., 0., -0.05], hopf_tspan, (-0.05,  hopf_ω, hopf_A))
	hopf_mu_0_1_system = ODEProblem(hopf_normal_mu!, [2., 0., 0.], hopf_tspan, (0,  hopf_ω, hopf_A))
	hopf_mu_005_1_system = ODEProblem(hopf_normal_mu!, [2., 0., 0.05], hopf_tspan, (0.05,  hopf_ω, hopf_A))
	hopf_mu_005_2_system = ODEProblem(hopf_normal_mu!, [0.1, 0., 0.05], hopf_tspan, (0.05,  hopf_ω, hopf_A))
	hopf_mu_01_1_system = ODEProblem(hopf_normal_mu!, [2., 0., 0.1], hopf_tspan, (0.1,  hopf_ω, hopf_A))
	hopf_mu_01_2_system = ODEProblem(hopf_normal_mu!, [0.1, 0., 0.1], hopf_tspan, (0.1,  hopf_ω, hopf_A))
	hopf_mu_015_1_system = ODEProblem(hopf_normal_mu!, [2., 0., 0.15], hopf_tspan, (0.15,  hopf_ω, hopf_A))
	hopf_mu_015_2_system = ODEProblem(hopf_normal_mu!, [0.1, 0., 0.15], hopf_tspan, (0.15,  hopf_ω, hopf_A))
	hopf_mu_02_1_system = ODEProblem(hopf_normal_mu!, [2., 0., 0.2], hopf_tspan, (0.2,  hopf_ω, hopf_A))
	hopf_mu_02_2_system = ODEProblem(hopf_normal_mu!, [0.1, 0., 0.2], hopf_tspan, (0.2,  hopf_ω, hopf_A))


	hopf_mu_neg02_1_solution = solve(hopf_mu_neg02_1_system, Tsit5())
	hopf_mu_neg015_1_solution = solve(hopf_mu_neg015_1_system, Tsit5())
	hopf_mu_neg01_1_solution = solve(hopf_mu_neg01_1_system, Tsit5())
	hopf_mu_neg005_1_solution = solve(hopf_mu_neg005_1_system, Tsit5())
	hopf_mu_0_1_solution = solve(hopf_mu_0_1_system, Tsit5())
	hopf_mu_005_1_solution = solve(hopf_mu_005_1_system, Tsit5())
	hopf_mu_005_2_solution = solve(hopf_mu_005_2_system, Tsit5())
	hopf_mu_01_1_solution = solve(hopf_mu_01_1_system, Tsit5())
	hopf_mu_01_2_solution = solve(hopf_mu_01_2_system, Tsit5())
	hopf_mu_015_1_solution = solve(hopf_mu_015_1_system, Tsit5())
	hopf_mu_015_2_solution = solve(hopf_mu_015_2_system, Tsit5())
	hopf_mu_02_1_solution = solve(hopf_mu_02_1_system, Tsit5())
	hopf_mu_02_2_solution = solve(hopf_mu_02_2_system, Tsit5())
end

# ╔═╡ 1b6e5c7b-6ea7-452e-b0f4-d569db4994e1
begin
	plot(hopf_mu_neg02_1_solution, idxs=(3, 1,2), color = :blue)
	plot!(hopf_mu_neg015_1_solution, idxs=(3, 1,2), color = :blue)
	plot!(hopf_mu_neg01_1_solution, idxs=(3, 1,2), color = :blue)
	plot!(hopf_mu_neg005_1_solution, idxs=(3, 1,2), color = :blue)
	plot!(hopf_mu_0_1_solution, idxs=(3, 1,2), color = :blue)
	plot!(hopf_mu_005_1_solution, idxs=(3, 1,2), color = :blue)
	plot!(hopf_mu_005_2_solution, idxs=(3, 1,2), color = :red)
	plot!(hopf_mu_01_1_solution, idxs=(3, 1,2), color = :blue)
	plot!(hopf_mu_01_2_solution, idxs=(3, 1,2), color = :red)
	plot!(hopf_mu_015_1_solution, idxs=(3, 1,2), color = :blue)
	plot!(hopf_mu_015_2_solution, idxs=(3, 1,2), color = :red)
	plot!(hopf_mu_02_1_solution, idxs=(3, 1,2), color = :blue)
	plot!(hopf_mu_02_2_solution, idxs=(3, 1,2), color = :red, zlabel = L"x", ylabel = L"y", xlabel = L"\mu", camera = (30,10), legend = false)
end

# ╔═╡ 1db7d38a-e095-49e7-aaac-db2c80285303
# Plots.savefig("./hopf_normal_actual.pdf")

# ╔═╡ c65722a9-d01e-4471-9c2c-e38f94dec6fc
md"""
# Stability testing
"""

# ╔═╡ de88e9ef-3f97-4285-9949-12c2a6481bcd
begin
	lorenz_polys_1 = collect(0:7)
	lorenz_polys_2 = collect(0:5)
	lorenz_polys_3 = collect(0:3)
	lorenz_polys_4 = collect(0:2)
	lorenz_polys_5 = collect(0:7)
	lorenz_polys_6 = collect(0:7)
	lorenz_func1_1 = []
	lorenz_func1_2 = []
	lorenz_func1_3 = []
	lorenz_func1_4 = []
	lorenz_func1_5 = [sin_cos_function_library(10)..., exp]
	lorenz_func1_6 = [exp]
	
	Θ_lorenz_01_tv_1 = Θ(X_lorenz_01_td; polynomials=lorenz_polys_1)
	Θ_lorenz_01_tv_2 = Θ(X_lorenz_01_td; polynomials=lorenz_polys_2)
	Θ_lorenz_01_tv_3 = Θ(X_lorenz_01_td; polynomials=lorenz_polys_3)
	Θ_lorenz_01_tv_4 = Θ(X_lorenz_01_td; polynomials=lorenz_polys_4)
	Θ_lorenz_01_tv_5 = Θ(X_lorenz_01_td; polynomials=lorenz_polys_5, functions1 = lorenz_func1_5)
	Θ_lorenz_01_tv_6 = Θ(X_lorenz_01_td; polynomials=lorenz_polys_6, functions1 = lorenz_func1_6)
end

# ╔═╡ b1636e3d-14c9-423c-8425-1cdd4ec5ce3f
begin 
	lorenz_polys_7 = [0]
	lorenz_polys_8 = [0]
	lorenz_polys_9 = [0]
	
	lorenz_func1_7 = [sin_cos_function_library(10)...]
	lorenz_func1_8 = [sin_cos_function_library(50)...]
	lorenz_func1_9 = [sin_cos_function_library(100)...]

	Θ_lorenz_01_tv_7 = Θ(X_lorenz_01_td; polynomials=lorenz_polys_7, functions1 = lorenz_func1_7)
	Θ_lorenz_01_tv_8 = Θ(X_lorenz_01_td; polynomials=lorenz_polys_8, functions1 = lorenz_func1_8)
	Θ_lorenz_01_tv_9 = Θ(X_lorenz_01_td; polynomials=lorenz_polys_9, functions1 = lorenz_func1_9)
end

# ╔═╡ c20472a4-502f-465e-99bd-9b2b2529f923
Ξ_lorenz_01_tv_1 = sparse_regression(Θ_lorenz_01_tv_1, Xdot_lorenz_01_td, 0.25, 10)

# ╔═╡ 19efb23b-cd9b-4d17-808b-e41fb68c6f0a
measure_accuracy(Ξ_lorenz_01_tv_1, X_lorenz_01_td_test, Xdot_lorenz_01_td_test; polynomials=lorenz_polys_1)

# ╔═╡ 28959896-8c9a-4419-849f-df73b6b3efd0
Ξ_lorenz_01_tv_2 = sparse_regression(Θ_lorenz_01_tv_2, Xdot_lorenz_01_td, 0.25, 10)

# ╔═╡ 4a3b17ca-5489-4598-8ddf-e8934f229c8d
measure_accuracy(Ξ_lorenz_01_tv_2, X_lorenz_01_td_test, Xdot_lorenz_01_td_test; polynomials=lorenz_polys_2)

# ╔═╡ fbc9e493-3626-4a32-b6e6-76e3fd2347ba
Ξ_lorenz_01_tv_3 = sparse_regression(Θ_lorenz_01_tv_3, Xdot_lorenz_01_td, 0.25, 10)

# ╔═╡ 60ee4f2a-441a-42d7-9b12-c8407c153d96
measure_accuracy(Ξ_lorenz_01_tv_3, X_lorenz_01_td_test, Xdot_lorenz_01_td_test; polynomials=lorenz_polys_3)

# ╔═╡ ba00669c-f5d0-400c-b373-52c7c3f4bfcd
Ξ_lorenz_01_tv_4 = sparse_regression(Θ_lorenz_01_tv_4, Xdot_lorenz_01_td, 0.25, 10)

# ╔═╡ 3416bc4d-dfae-4d70-8597-a6d0b0e1a9d6
measure_accuracy(Ξ_lorenz_01_tv_4, X_lorenz_01_td_test, Xdot_lorenz_01_td_test; polynomials=lorenz_polys_4)

# ╔═╡ a951f5dc-ce5a-4394-8f39-6a8f76b4f536
# Ξ_lorenz_01_tv_5 = sparse_regression(Θ_lorenz_01_tv_5, Xdot_lorenz_01_tv, 0.25, 10)

# ╔═╡ f1ddab03-bd72-4bf6-b4ce-db57b6203df1
# Ξ_lorenz_01_tv_6 = sparse_regression(Θ_lorenz_01_tv_6, Xdot_lorenz_01_tv, 0.25, 10)

# ╔═╡ 3398e448-c159-4e75-8147-26c0aab6fda2
Ξ_lorenz_01_tv_7 = sparse_regression(Θ_lorenz_01_tv_7, Xdot_lorenz_01_td, 0.25, 10)

# ╔═╡ 5875a499-1c42-41ff-89f6-370562585d1f
measure_accuracy(Ξ_lorenz_01_tv_7, X_lorenz_01_td_test, Xdot_lorenz_01_td_test; polynomials=lorenz_polys_7, functions1 = lorenz_func1_7)

# ╔═╡ 34d575f7-1bfa-439b-92e7-005a1c7c9902
Ξ_lorenz_01_tv_8 = sparse_regression(Θ_lorenz_01_tv_8, Xdot_lorenz_01_td, 0.25, 10)

# ╔═╡ 2e9c1c02-8ce2-44b5-af91-fb7dae42a64b
measure_accuracy(Ξ_lorenz_01_tv_8, X_lorenz_01_td_test, Xdot_lorenz_01_td_test; polynomials=lorenz_polys_8, functions1 = lorenz_func1_8)

# ╔═╡ f694350a-c4ea-482a-8d4c-f59345e077b5
# Ξ_lorenz_01_tv_9 = sparse_regression(Θ_lorenz_01_tv_9, Xdot_lorenz_01_tv, 0.25, 10)

# ╔═╡ 8e8a48ec-54e4-411d-a5ce-ab6d07db890f
lorenz_01_tv_system_1 = solve_constructed_system(Ξ_lorenz_01_tv_1, [-8., 7, 27], (0., 1000.); polynomials = lorenz_polys_1, functions1 = lorenz_func1_1)

# ╔═╡ c72fe074-475d-4992-a07d-6feef81aa78b
lorenz_01_tv_system_2 = solve_constructed_system(Ξ_lorenz_01_tv_2, [-8., 7, 27], (0., 1000.); polynomials = lorenz_polys_2, functions1 = lorenz_func1_2)

# ╔═╡ 6fca07d1-960c-4cf7-89b6-25c86c6ef06e
lorenz_01_tv_system_3 = solve_constructed_system(Ξ_lorenz_01_tv_3, [-8., 7, 27], (0., 1000.); polynomials = lorenz_polys_3, functions1 = lorenz_func1_3)

# ╔═╡ 7405be92-6b5c-4a1e-a982-71f67ffd8388
lorenz_01_tv_system_4 = solve_constructed_system(Ξ_lorenz_01_tv_4, [-8., 7, 27], (0., 1000.); polynomials = lorenz_polys_4, functions1 = lorenz_func1_4)

# ╔═╡ d886ae5c-54ca-4c41-a945-0dee9a66f469
# lorenz_01_tv_system_5 = solve_constructed_system(Ξ_lorenz_01_tv_5, [-8., 7, 27], (0., 1000.); polynomials = lorenz_polys_5, functions1 = lorenz_func1_5)

# ╔═╡ 0bcac161-3a8e-448f-8a3e-46a27b00430d
# lorenz_01_tv_system_6 = solve_constructed_system(Ξ_lorenz_01_tv_6, [-8., 7, 27], (0., 1000.); polynomials = lorenz_polys_6, functions1 = lorenz_func1_6)

# ╔═╡ cc9b8403-6a80-497f-aa68-6e24df23ba55
lorenz_01_tv_system_7 = solve_constructed_system(Ξ_lorenz_01_tv_7, [-8., 7, 27], (0., 1000.); polynomials = lorenz_polys_7, functions1 = lorenz_func1_7)

# ╔═╡ d232946b-2e3a-4435-94cc-1bcc9b93030e
lorenz_01_tv_system_8 = solve_constructed_system(Ξ_lorenz_01_tv_8, [-8., 7, 27], (0., 1000.); polynomials = lorenz_polys_8, functions1 = lorenz_func1_8)

# ╔═╡ 7ef9229f-4d3f-48d4-b1ae-f2665edaa56d
print_constructed_system(Ξ_lorenz_01_tv_1; polynomials = lorenz_polys_1, functions1 = lorenz_func1_1)

# ╔═╡ ea9821b3-6209-4f33-9f9b-096963622622
begin
	lorenz_1_sol = solve(lorenz_01_tv_system_1, Tsit5())
	plot(lorenz_1_sol, idxs=(1,2,3), xlabel=L"x", ylabel=L"y", zlabel=L"z", camera=(40,30), legend = false)
end

# ╔═╡ 84cfffc2-8116-4258-b21e-ec3a781c50d6
# Plots.savefig("./lorenz_1_sol.pdf")

# ╔═╡ 04047955-f751-488b-b5f4-59c549227c8d
print_constructed_system(Ξ_lorenz_01_tv_2; polynomials = lorenz_polys_2, functions1 = lorenz_func1_2)

# ╔═╡ f050abbf-14ad-4946-b741-41263a842ef8
begin
	lorenz_2_sol = solve(lorenz_01_tv_system_2, Tsit5())
	plot(lorenz_2_sol, idxs=(1,2,3), xlabel=L"x", ylabel=L"y", zlabel=L"z", camera=(40,30), legend = false)
end

# ╔═╡ e94bfe49-355b-40bc-bd9c-af32ceeb9e72
# Plots.savefig("./lorenz_2_sol.pdf")

# ╔═╡ 9ddf90fc-2cf3-480c-8ea8-74c46176325f
print_constructed_system(Ξ_lorenz_01_tv_3; polynomials = lorenz_polys_3, functions1 = lorenz_func1_3)

# ╔═╡ 6cffde55-f336-413d-808c-2689cbc56ce5
begin
	lorenz_3_sol = solve(lorenz_01_tv_system_3, Tsit5())
	plot(lorenz_3_sol, idxs=(1,2,3), xlabel=L"x", ylabel=L"y", zlabel=L"z", camera=(40,30), legend = false)
end

# ╔═╡ a46d0b09-de77-4ccb-80e0-56d530a5c0b9
# Plots.savefig("./lorenz_3_sol.pdf")

# ╔═╡ 094e375d-a718-4b69-b996-0449cdaa9332
print_constructed_system(Ξ_lorenz_01_tv_4; polynomials = lorenz_polys_4, functions1 = lorenz_func1_4)

# ╔═╡ 2f911b2f-3f92-4505-be88-d8616833d812
begin
	lorenz_4_sol = solve(lorenz_01_tv_system_4, Tsit5())
	plot(lorenz_4_sol, idxs=(1,2,3), xlabel=L"x", ylabel=L"y", zlabel=L"z", camera=(40,30), legend = false)
end

# ╔═╡ 1f57d65d-3844-4f10-8972-bf20e76c067d
# Plots.savefig("./lorenz_4_sol.pdf")

# ╔═╡ 14fb43e1-9473-4fb9-9fdf-840f0e123d54
# print_constructed_system(Ξ_lorenz_01_tv_5; polynomials = lorenz_polys_5, functions1 = lorenz_func1_5)

# ╔═╡ 1c0d3997-0f39-48c0-8004-1a0f785a2174
# print_constructed_system(Ξ_lorenz_01_tv_6; polynomials = lorenz_polys_6, functions1 = lorenz_func1_6)

# ╔═╡ bfc51a88-1d1a-4fb4-8f5c-4e89a854eaeb
print_constructed_system(Ξ_lorenz_01_tv_7; polynomials = lorenz_polys_7, functions1 = lorenz_func1_7)

# ╔═╡ d9ff35ff-61bb-48ab-a5e6-74fead0d4593
begin
	lorenz_7_sol = solve(lorenz_01_tv_system_7, Tsit5())
	plot(lorenz_7_sol, idxs=(1,2,3), xlabel=L"x", ylabel=L"y", zlabel=L"z", camera=(40,30), legend = false)
end

# ╔═╡ 8a0a5622-4828-4749-ba07-fc2ed554872a
# Plots.savefig("./lorenz_7_sol.pdf")

# ╔═╡ 97dc0b67-a357-40d3-9ccf-edade7c5d409
print_constructed_system(Ξ_lorenz_01_tv_8; polynomials = lorenz_polys_8, functions1 = lorenz_func1_8)

# ╔═╡ b8e77d92-a72e-4129-bf52-b4cec66dc559
begin
	lorenz_8_sol = solve(lorenz_01_tv_system_8, Tsit5())
	plot(lorenz_8_sol, idxs=(1,2,3), xlabel=L"x", ylabel=L"y", zlabel=L"z", camera=(40,30), legend = false)
end

# ╔═╡ b74f69c5-75a8-4ba5-9349-af9d369729a5
# Plots.savefig("./lorenz_8_sol.pdf")

# ╔═╡ 2df97f42-ec7f-4c6f-9431-3b69a8cb8024
# begin
# 	test1 = ODEProblem(lorenz!, [-8., 7., 27.], (0.,100.), (10., 28., 8/3))
# 	saving = collect(0:0.01:100)
# 	res = solve(test1, DP5(), dt = 0.01, saveat=saving)

# 	res.u
# end

# ╔═╡ daa284d3-f9bc-45b7-9c53-df70ab60c204
# begin 
# 	test = active_SINDy(lorenz!, 0.001, (0.,100.), [-8., 7., 27.], (10., 28., 8/3), Tsit5(), 0.01, (0., 100.), TvDiffParams(), 0.0001, 10, 0.98, 0.0001; polynomials = collect(0:5))
# end

# ╔═╡ 9d481e6c-feef-4f85-a9ea-d7def332580d
# begin 
# 	test_system = solve_constructed_system(test, [-8., 7., 27.], (0., 100.); polynomials = collect(0:5))
# 	test_sol = solve(test_system, Tsit5())
# end

# ╔═╡ ca23d1c7-4efd-4097-a6e4-cf6e5cd84fa5
# plot(test_sol, idxs=(1,2,3))

# ╔═╡ 78670281-9464-43a9-8a27-d5e0d05066dc
# print_constructed_system(test; polynomials = collect(0:5))

# ╔═╡ 6083eb6b-1606-4f23-922c-f1516313b064
lorenz_001_error_test = automatic_lambda_fitting(X_lorenz_001_td, Xdot_lorenz_001_td, X_lorenz_001_td_test, Xdot_lorenz_001_td_test, Θ_lorenz_001_td;polynomials=collect(0:7), functions1 = [sin_cos_function_library(10)...], expected_multiple = 2.5)

# ╔═╡ 667d7262-8668-4045-b25a-91f483f23326
lorenz_01_error_test = automatic_lambda_fitting(X_lorenz_01_td, Xdot_lorenz_01_td, X_lorenz_01_td_test, Xdot_lorenz_01_td_test, Θ_lorenz_01_td; polynomials=lorenz_01_td_polynomials, functions1 = lorenz_01_td_functions1)

# ╔═╡ af1e66c9-c59e-466e-8bd1-1e2e963873f6
lorenz_1_error_test = automatic_lambda_fitting(X_lorenz_1_td, Xdot_lorenz_1_td, X_lorenz_1_td_test, Xdot_lorenz_1_td_test, Θ_lorenz_1_td; polynomials=lorenz_1_td_polynomials, functions1 = lorenz_1_td_functions1)

# ╔═╡ 07559337-2d5a-4fcb-8d77-7e8f15cb370a
print_constructed_system(lorenz_001_error_test[2]; polynomials = collect(0:7), functions1 = [sin_cos_function_library(10)...])

# ╔═╡ 8333a191-c1af-41a1-9ace-8f2cdbd457af
measure_accuracy(lorenz_001_error_test[2], X_lorenz_001_td_test, Xdot_lorenz_001_td_test; polynomials = collect(0:7), functions1 = [sin_cos_function_library(10)...])

# ╔═╡ c624dcf0-c0c3-41e8-a703-db61c948f3ae
print_constructed_system(lorenz_01_error_test[2]; polynomials=lorenz_01_td_polynomials, functions1 = lorenz_01_td_functions1)

# ╔═╡ 39e184ba-b5c3-4c8c-84f2-56566fd3c415
measure_accuracy(lorenz_01_error_test[2], X_lorenz_01_td_test, Xdot_lorenz_01_td_test; polynomials=lorenz_01_td_polynomials, functions1 = lorenz_01_td_functions1)

# ╔═╡ 299ce6da-d9a7-4882-9a4b-0d74272653f2
print_constructed_system(lorenz_1_error_test[2]; polynomials=lorenz_1_td_polynomials, functions1 = lorenz_1_td_functions1)

# ╔═╡ d4e722d5-7abe-44ea-8c41-074c7ea4303a
measure_accuracy(lorenz_1_error_test[2], X_lorenz_1_td_test, Xdot_lorenz_1_td_test; polynomials=lorenz_1_td_polynomials, functions1 = lorenz_1_td_functions1)

# ╔═╡ 5bb9e6d0-f322-498f-80a9-09cce2303b35
md"""
# Three Body Problem 
"""

# ╔═╡ 34056977-9406-458a-a7bc-191898555d97
function three_body!(du, u, p, t)
	# 18 dimensional system as 9 positions (x, y, z) for each body and 9 velocities  
	# first nine entries of du are positions, giving the x, y, z for the first body, and then the x, y, z for the second, etc. 
	# Last nine are the velocities in the same order
	#
	# Equations from: https://en.wikipedia.org/wiki/Three-body_problem
	
	(G, m_1, m_2, m_3) = p # unpack p

	# create position and velocity vectors
	r_1 = u[1:3]
	r_2 = u[4:6]
	r_3 = u[7:9]
	dot_r_1 = u[10:12]
	dot_r_2 = u[13:15]
	dot_r_3 = u[16:18]

	# compute accelerations
	ddot_r_1 = -1 * G * m_2 * (r_1 - r_2)/(norm(r_1 - r_2)^3) - G * m_3 * (r_1 - r_3)/(norm(r_1-r_3)^3)
	ddot_r_2 = -1 * G * m_3 * (r_2 - r_3)/(norm(r_2 - r_3)^3) - G * m_1 * (r_2 - r_1)/(norm(r_2 - r_1)^3)
	ddot_r_3 = -1 * G * m_1 * (r_3 - r_1)/(norm(r_3 - r_1)^3) - G * m_2 * (r_3 - r_2)/(norm(r_3 - r_2)^3)

	# derivatives of position will be velocities
	du[1] = u[10]
	du[2] = u[11]
	du[3] = u[12]
	du[4] = u[13]
	du[5] = u[14]
	du[6] = u[15]
	du[7] = u[16]
	du[8] = u[17]
	du[9] = u[18]

	# derivatives of velocities will be accelerations 
	du[10] = ddot_r_1[1]
	du[11] = ddot_r_1[2]
	du[12] = ddot_r_1[3]
	du[13] = ddot_r_2[1]
	du[14] = ddot_r_2[2]
	du[15] = ddot_r_2[3]
	du[16] = ddot_r_3[1]
	du[17] = ddot_r_3[2]
	du[18] = ddot_r_3[3]
end # three_body!

# ╔═╡ 2578b346-653f-4f9c-b89a-d66011db25af
begin 
	# wikipedia suggests interesting initial condition is identical bodies at vertices of scalene triangle (zero initial velocity)
	tbp_u0 = [
		0., 0., 0., # body 1 pos
		3., 4., sqrt(11), # body 2 pos
		5., 7., 1., # body 3 pos
		0., 0., 0., # body 1 vel 
		0., 0., 0., # body 2 vel
		0., 0., 0.  # body 3 vel
	]
	
	tbp_tspan = (0., 1000.)

	tbp_p = (1e-1, 1., 1., 1.) # G = 1 now
	
	tbp_sys = ODEProblem(three_body!, tbp_u0, tbp_tspan, tbp_p)
	
	tbp_dt = 0.001
	tbp_times = collect((tbp_tspan[1]):tbp_dt:(tbp_tspan[2]))

	tbp_sol = solve(tbp_sys, Tsit5(), dt = tbp_dt, saveat = tbp_times)
end

# ╔═╡ 52775ee9-1119-4328-829d-b1f5c95045fa
begin 
	plot(tbp_sol, idxs = (1,2,3), label = "body 1", color = colorant"#b616fa")
	plot!(tbp_sol, idxs = (4,5,6), label = "body 2", color = colorant"#1693fa")
	plot!(tbp_sol, idxs = (7,8,9), label = "body 3", color = colorant"#e66f07", xlabel = L"x", ylabel = L"y", zlabel = L"z")
end

# ╔═╡ 3334a8d6-ee5a-4d1e-ab3a-ecaa53b5de75
noisy_tbp = add_noise_to_sample(tbp_sol.u, 0.001)

# ╔═╡ 8adac501-90b6-4d6e-8766-14c8cd146e29
# X_tbp, Xdot_tbp = constructX(noisy_tbp, TvDiffParams(; iterations = 10, α = 2e-1, scale = "small", ε = 1e6, perform_tvdiff = true), tbp_sol.t)

# ╔═╡ 35184db0-dc43-48bc-a60b-28fbbbdbaa56
# X_tb, Xdot_tb, X_tb_test, Xdot_tb_test = get_testing_X(X_tbp, Xdot_tbp, 0.1)

# ╔═╡ 2047d73f-1434-4501-a117-116305eef590
function generate_inverse_32_functions()
	functions = []
	pairs_of_variables = combinations(collect(0:2), 2)
	for pair in pairs_of_variables 
		for index in 1:9
			# function f(u...)
			# 	uvec = collect(u)

			# 	println(uvec[(pair[1]+1):(pair[1]+3)])
			# 	println(uvec[(pair[2]+1):(pair[2]+3)])

			# 	println(norm(uvec[(pair[1]+1):(pair[1]+3)] - uvec[(pair[2]+1):(pair[2]+3)]))
			# end
			
			f(u...) = u[index] / ( norm(collect(u[(pair[1]+1):(pair[1]+3)]) - collect(u[(pair[2]+1):(pair[2]+3)]))^3 )
			push!(functions, f)
		end # for index 
	end # for pair 
	functions
end

# ╔═╡ c3b901d7-8cbc-4624-af35-69dd37c1f67d
# Θ_tbp = Θ(X_tb; polynomials = collect(0:2), functions = [generate_inverse_32_functions()...])

# ╔═╡ a7c0ef42-6c1b-42bc-bf98-a3cd1f26b71a
# tbp_fit = automatic_lambda_fitting(X_tb, Xdot_tb, X_tb_test, Xdot_tb_test, Θ_tbp; polynomials=collect(0:2), functions = [generate_inverse_32_functions()...])

# ╔═╡ Cell order:
# ╠═71998b4b-1bfc-4c3d-8c7e-5741de700162
# ╠═a0a26a0d-a7e9-4168-996b-51f44fd96dad
# ╠═0e6831be-f2e1-11ee-1ecc-1f97692f0c45
# ╠═aa29b69a-5750-4a5f-8e8c-fdc49d9a45b6
# ╠═1511b58f-0ec9-4dc3-90ec-d3169cf9189e
# ╠═ee6d98bc-97dc-48f4-92fe-79a69afd51ca
# ╠═f7fb3eaf-0473-4074-a99c-892655fe12b7
# ╠═c2016647-3b10-4041-a584-920ab5c116d1
# ╠═e4a6e6f9-43c6-489e-9ae2-451ee1f9e0df
# ╠═0d836d4e-d9d2-4fa1-9983-e541fb58ab04
# ╟─61ccdc43-c12f-4b59-b746-316c82b342a5
# ╠═f2f7d603-4836-4010-95c7-10def9254590
# ╠═4cbdb8af-68ca-421f-853c-5bd8c4f4b464
# ╠═7970ebf4-794c-41d0-979b-7e2f83163fad
# ╠═b2a4fb51-90a0-49af-9a3e-f0d21a54495d
# ╠═963a2ee7-fb5b-4ce9-88e1-da5d3f2ddbbd
# ╠═d3c346a0-8970-4e12-aaf8-5884a0ef435c
# ╠═d74cdb2f-7bfe-452c-b96b-1b9bf33c29f2
# ╠═3413048a-a378-414f-bd4a-f21beaf31289
# ╠═69042648-7393-453d-bd59-ad2b345e1ada
# ╠═c948a893-a13c-45c2-a66b-8eff8ab85d68
# ╠═369c1bca-dbed-4d84-9860-10b08e5e825c
# ╠═71864910-7252-4408-be07-9bdd95192926
# ╠═afc779f4-0b87-4de6-b6ad-bc688b20b318
# ╠═a347d212-183c-4343-8171-1bc3a2ea7ef6
# ╠═383b7116-0612-4626-a318-3dc72ff5e35c
# ╠═b5b84157-9a8e-4cf7-97f0-5e8f168d15b9
# ╠═940b3efb-bbe9-4130-879d-44b12917e209
# ╠═5bbf1cdc-7568-4485-b81f-2d0a9754244b
# ╠═83863d85-1276-407e-aef5-6ef41e802fe2
# ╠═1c29a5e0-014d-4f3f-be6f-560dd1533863
# ╠═cb9c98c7-e57b-424d-980b-24013fbbeae9
# ╠═ca03bce6-6477-4235-83f5-9bfdbafa1739
# ╠═84671725-6117-4ea1-80e7-2ebf578f7b79
# ╠═47e7f92a-ced3-4632-81c8-1309ef873637
# ╠═8e5c1398-c6ef-4b68-991d-d3045d3be84e
# ╠═6206f493-513c-4be8-bfee-6c8ca16609dc
# ╠═1c18d80f-0e77-4676-9302-42d5cb396184
# ╠═39d59c5c-38d0-49d9-ae6c-079ac9f81901
# ╠═3fc03792-830e-42f9-b4b7-fbd155cc4142
# ╠═5eb64ceb-aeac-4ea9-a6c6-ddf49937a860
# ╠═25166eec-a00d-4c5d-8ac6-f70208706613
# ╠═529c0dc4-b9e5-40d9-bcc0-872792298fc0
# ╠═1a5465b5-bd3c-4c9f-9921-983f99492ed5
# ╠═7d44b26e-b916-4145-8147-8c24df66b5cd
# ╠═b6350c06-cfb1-4ea4-bea6-1d95b1d3c5fd
# ╠═d9050f03-98d6-477e-b18f-07034bf177b3
# ╠═b5a165c5-6371-4354-be1c-eedf8d322de4
# ╠═3e02c732-d1b6-4768-86a8-fb244372c1bc
# ╠═590df629-ec43-4a0d-8081-21a407e16df3
# ╠═c027b7db-7d2d-4fc7-89b9-ffb5d4501a75
# ╠═0d0cd6f9-a0c7-49f1-b4c4-c8d76ffaa2a3
# ╠═020e0c54-955b-4117-9b3e-251cf3d7b7d5
# ╠═5fa8c07b-5bd3-4136-95b2-09ea2eedbab1
# ╠═6246c8a2-638e-4eb6-88a9-2a335f86cd88
# ╠═d1a0d67e-a6a5-400b-b698-9257b56c6b99
# ╠═e61f7ee9-2f25-45e2-9952-7abc8123537e
# ╠═e70e3f6e-437e-4ee7-a4a8-76a72791fe37
# ╠═b2c01934-98a3-4daa-b957-ede0a6c3e20f
# ╠═9cd257db-474d-4be6-be81-cd847be71cef
# ╠═438a33b0-2d93-4bc9-aab7-b789d89eb010
# ╠═284b77f9-1f84-4f79-84c5-b2299f8ff674
# ╠═c0dc5518-edc3-4564-9f00-0632691c5a9c
# ╠═f10bb9b2-0c9e-4c9b-8057-003a91da4128
# ╠═8c207abe-db98-43e0-9df8-5b662fac803f
# ╠═7a70e0fe-728e-411d-9d7c-adb5d39efdb1
# ╠═23f965a8-29e6-42e3-9fa5-7d02b20a1cc3
# ╠═eb9455c5-238e-44e9-ab8b-5872d9b60d9e
# ╠═2a6f9358-b196-418e-adf7-1d6348c92ee1
# ╠═27efd12e-5608-40c1-889d-a247562ece3a
# ╠═d59d0489-bd49-4ac8-a8eb-2a95ed1579b4
# ╠═7b4cb49e-d263-4711-bb1d-b90aa2e1df91
# ╠═12ec1247-f93a-4f49-a847-61dc9bf5bb0f
# ╠═e22bfbb2-bb20-441d-8b40-82d5d5760d89
# ╠═1c029b31-5d6a-423f-90fc-06c902977962
# ╠═03c887da-dd03-49f0-9bcd-5bc33a619a7a
# ╠═e90d0073-3530-4662-80be-236675c648c3
# ╠═d04e4a91-1feb-4012-bff0-51dcf9a04f70
# ╠═bd8867fe-0a1c-4e3d-910d-32acba5ef2af
# ╠═2a2d55ac-2a06-4f02-8847-82ba264bc40c
# ╠═6f6c0019-ed7d-4a31-aa06-7291cd6cd532
# ╠═1682c040-1735-48d9-93a5-1a1900a1d1cb
# ╠═8eb49379-7b2a-4edb-924f-b6417450febf
# ╠═bfeee94b-b98d-4345-8b8a-efd9eda09c22
# ╠═e00de518-1c38-452c-8b79-4adbfeecea53
# ╠═e4fe9390-6955-405c-86cc-cc9defb9ca39
# ╠═c1ea010b-9549-43d7-8aa2-9b6a92c00d82
# ╠═0814a304-c439-40e5-84a9-4e501aaf2326
# ╠═3d9a773b-39a9-41bd-b9a9-4c79a33992de
# ╠═b4a92457-6529-4ab9-872f-0bc2ba099ff4
# ╠═7a8cc3aa-5b30-46b9-83db-3f47d0377c13
# ╠═e3f29ead-9e94-48e9-9dfa-a07fba5a0fef
# ╠═9ced946d-5f6c-4801-a204-ab2a319cb861
# ╠═f35be9d2-efea-4f0a-869d-25bd24f946a0
# ╠═e0003431-ffef-429e-80c8-ba1c7f1774d2
# ╠═d333670f-c9c6-461f-81d8-84f7bb974d85
# ╠═768ff383-57d4-4203-8c25-4cba99c993bf
# ╠═6168c7b1-6f44-4f80-ba98-257e9d447e98
# ╠═37dc2ad0-d356-4f78-9841-80836d9f0692
# ╠═3781cba0-2f81-471d-accc-de4d8546c11c
# ╠═20d272ea-c701-4edf-a393-a48f514465fc
# ╠═261776eb-e90e-4471-8276-cb1b934cde28
# ╠═f66ceb82-0780-4149-a50d-5f6c1ffe6e6b
# ╠═f6cb782c-934b-4342-b008-10be1fda2021
# ╠═379bba48-ebc0-415d-bd54-25707b97bab1
# ╠═51281685-dea0-4c19-b6c6-893b265b051c
# ╠═ff4ae999-17c3-4487-abb8-28204aa7b73f
# ╠═1b6e5c7b-6ea7-452e-b0f4-d569db4994e1
# ╠═1db7d38a-e095-49e7-aaac-db2c80285303
# ╟─c65722a9-d01e-4471-9c2c-e38f94dec6fc
# ╠═de88e9ef-3f97-4285-9949-12c2a6481bcd
# ╠═b1636e3d-14c9-423c-8425-1cdd4ec5ce3f
# ╠═c20472a4-502f-465e-99bd-9b2b2529f923
# ╠═19efb23b-cd9b-4d17-808b-e41fb68c6f0a
# ╠═28959896-8c9a-4419-849f-df73b6b3efd0
# ╠═4a3b17ca-5489-4598-8ddf-e8934f229c8d
# ╠═fbc9e493-3626-4a32-b6e6-76e3fd2347ba
# ╠═60ee4f2a-441a-42d7-9b12-c8407c153d96
# ╠═ba00669c-f5d0-400c-b373-52c7c3f4bfcd
# ╠═3416bc4d-dfae-4d70-8597-a6d0b0e1a9d6
# ╠═a951f5dc-ce5a-4394-8f39-6a8f76b4f536
# ╠═f1ddab03-bd72-4bf6-b4ce-db57b6203df1
# ╠═3398e448-c159-4e75-8147-26c0aab6fda2
# ╠═5875a499-1c42-41ff-89f6-370562585d1f
# ╠═34d575f7-1bfa-439b-92e7-005a1c7c9902
# ╠═2e9c1c02-8ce2-44b5-af91-fb7dae42a64b
# ╠═f694350a-c4ea-482a-8d4c-f59345e077b5
# ╠═8e8a48ec-54e4-411d-a5ce-ab6d07db890f
# ╠═c72fe074-475d-4992-a07d-6feef81aa78b
# ╠═6fca07d1-960c-4cf7-89b6-25c86c6ef06e
# ╠═7405be92-6b5c-4a1e-a982-71f67ffd8388
# ╠═d886ae5c-54ca-4c41-a945-0dee9a66f469
# ╠═0bcac161-3a8e-448f-8a3e-46a27b00430d
# ╠═cc9b8403-6a80-497f-aa68-6e24df23ba55
# ╠═d232946b-2e3a-4435-94cc-1bcc9b93030e
# ╠═7ef9229f-4d3f-48d4-b1ae-f2665edaa56d
# ╠═ea9821b3-6209-4f33-9f9b-096963622622
# ╠═84cfffc2-8116-4258-b21e-ec3a781c50d6
# ╠═04047955-f751-488b-b5f4-59c549227c8d
# ╠═f050abbf-14ad-4946-b741-41263a842ef8
# ╠═e94bfe49-355b-40bc-bd9c-af32ceeb9e72
# ╠═9ddf90fc-2cf3-480c-8ea8-74c46176325f
# ╠═6cffde55-f336-413d-808c-2689cbc56ce5
# ╠═a46d0b09-de77-4ccb-80e0-56d530a5c0b9
# ╠═094e375d-a718-4b69-b996-0449cdaa9332
# ╠═2f911b2f-3f92-4505-be88-d8616833d812
# ╠═1f57d65d-3844-4f10-8972-bf20e76c067d
# ╠═14fb43e1-9473-4fb9-9fdf-840f0e123d54
# ╠═1c0d3997-0f39-48c0-8004-1a0f785a2174
# ╠═bfc51a88-1d1a-4fb4-8f5c-4e89a854eaeb
# ╠═d9ff35ff-61bb-48ab-a5e6-74fead0d4593
# ╠═8a0a5622-4828-4749-ba07-fc2ed554872a
# ╠═97dc0b67-a357-40d3-9ccf-edade7c5d409
# ╠═b8e77d92-a72e-4129-bf52-b4cec66dc559
# ╠═b74f69c5-75a8-4ba5-9349-af9d369729a5
# ╠═2df97f42-ec7f-4c6f-9431-3b69a8cb8024
# ╠═daa284d3-f9bc-45b7-9c53-df70ab60c204
# ╠═9d481e6c-feef-4f85-a9ea-d7def332580d
# ╠═ca23d1c7-4efd-4097-a6e4-cf6e5cd84fa5
# ╠═78670281-9464-43a9-8a27-d5e0d05066dc
# ╠═6083eb6b-1606-4f23-922c-f1516313b064
# ╠═667d7262-8668-4045-b25a-91f483f23326
# ╠═af1e66c9-c59e-466e-8bd1-1e2e963873f6
# ╠═07559337-2d5a-4fcb-8d77-7e8f15cb370a
# ╠═8333a191-c1af-41a1-9ace-8f2cdbd457af
# ╠═c624dcf0-c0c3-41e8-a703-db61c948f3ae
# ╠═39e184ba-b5c3-4c8c-84f2-56566fd3c415
# ╠═299ce6da-d9a7-4882-9a4b-0d74272653f2
# ╠═d4e722d5-7abe-44ea-8c41-074c7ea4303a
# ╟─5bb9e6d0-f322-498f-80a9-09cce2303b35
# ╠═34056977-9406-458a-a7bc-191898555d97
# ╠═2578b346-653f-4f9c-b89a-d66011db25af
# ╠═52775ee9-1119-4328-829d-b1f5c95045fa
# ╠═3334a8d6-ee5a-4d1e-ab3a-ecaa53b5de75
# ╠═8adac501-90b6-4d6e-8766-14c8cd146e29
# ╠═35184db0-dc43-48bc-a60b-28fbbbdbaa56
# ╠═2047d73f-1434-4501-a117-116305eef590
# ╠═c3b901d7-8cbc-4624-af35-69dd37c1f67d
# ╠═a7c0ef42-6c1b-42bc-bf98-a3cd1f26b71a
