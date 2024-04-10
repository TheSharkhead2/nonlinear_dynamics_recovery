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
	plotly()
	# pyplot()
end

# ╔═╡ 51281685-dea0-4c19-b6c6-893b265b051c
begin 
	using PlotlyKaleido
	PlotlyKaleido.start()
	Plots.savefig("./found_hopf_bifurcation_0005.pdf")
end

# ╔═╡ a0a26a0d-a7e9-4168-996b-51f44fd96dad
# Pkg.add("LaTeXStrings")
# Pkg.add("PlotlyKaleido")

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
plot(lorenz_system, idxs=(1,2,3), legend = false, xlabel = "x", ylabel = "y", zlabel = "z", camera=(40,30))

# ╔═╡ b2a4fb51-90a0-49af-9a3e-f0d21a54495d
begin
	noisy_lorenz001 = add_noise_to_sample(lorenz_system.u, 0.01)
	noisy_lorenz01 = add_noise_to_sample(lorenz_system.u, 0.1)
end

# ╔═╡ d3c346a0-8970-4e12-aaf8-5884a0ef435c
X_lorenz_001, Xdot_lorenz_001 = constructX(noisy_lorenz001, lorenz_system.t)


# ╔═╡ d74cdb2f-7bfe-452c-b96b-1b9bf33c29f2
X_lorenz_01, Xdot_lorenz_01 = constructX(noisy_lorenz01, lorenz_system.t)

# ╔═╡ 69042648-7393-453d-bd59-ad2b345e1ada
X_lorenz_001_tv, Xdot_lorenz_001_tv = constructX(noisy_lorenz001, TvDiffParams(; iterations = 10, α = 2e-1, scale = "small", ε = 1e6, perform_tvdiff = true), lorenz_system.t)

# ╔═╡ c948a893-a13c-45c2-a66b-8eff8ab85d68
X_lorenz_01_tv, Xdot_lorenz_01_tv = constructX(noisy_lorenz01, TvDiffParams(; iterations = 10, α = 1, scale = "small", ε = 1e6, perform_tvdiff = true), lorenz_system.t)

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

# ╔═╡ 7b4cb49e-d263-4711-bb1d-b90aa2e1df91
Θ_lorenz_001 = Θ(X_lorenz_001; polynomials=[0, 1,2,3,4,5,6,7], functions1 = [sin_cos_function_library(10)...])

# ╔═╡ 12ec1247-f93a-4f49-a847-61dc9bf5bb0f
Θ_lorenz_01 = Θ(X_lorenz_01; polynomials = [0,1,2,3,4,5,6,7]) # , functions1 = [sin_cos_function_library(10)...]

# ╔═╡ e22bfbb2-bb20-441d-8b40-82d5d5760d89
Θ_lorenz_001_tv = Θ(X_lorenz_001_tv; polynomials=[0, 1,2,3,4,5,6,7], functions1 = [sin_cos_function_library(10)...])

# ╔═╡ 1c029b31-5d6a-423f-90fc-06c902977962
Θ_lorenz_01_tv = Θ(X_lorenz_01_tv; polynomials=[0, 1,2,3,4,5,6,7], functions1 = [sin_cos_function_library(10)...])

# ╔═╡ 03c887da-dd03-49f0-9bcd-5bc33a619a7a
Ξ_lorenz_001 = sparse_regression(Θ_lorenz_001, Xdot_lorenz_001, 0.25, 10)

# ╔═╡ e90d0073-3530-4662-80be-236675c648c3
Ξ_lorenz_01 = sparse_regression(Θ_lorenz_01, Xdot_lorenz_01, 0.025, 10)

# ╔═╡ d04e4a91-1feb-4012-bff0-51dcf9a04f70
Ξ_lorenz_001_tv = sparse_regression(Θ_lorenz_001_tv, Xdot_lorenz_001_tv, 0.25, 10)

# ╔═╡ bd8867fe-0a1c-4e3d-910d-32acba5ef2af
Ξ_lorenz_01_tv = sparse_regression(Θ_lorenz_01_tv, Xdot_lorenz_01_tv, 0.25, 10)

# ╔═╡ 2a2d55ac-2a06-4f02-8847-82ba264bc40c
print_constructed_system(Ξ_lorenz_001; polynomials=[0, 1,2,3,4,5,6,7], functions1 = [sin_cos_function_library(10)...])

# ╔═╡ 6f6c0019-ed7d-4a31-aa06-7291cd6cd532
print_constructed_system(Ξ_lorenz_001_tv; polynomials=[0, 1,2,3,4,5,6,7], functions1 = [sin_cos_function_library(10)...])

# ╔═╡ 1682c040-1735-48d9-93a5-1a1900a1d1cb
print_constructed_system(Ξ_lorenz_01_tv; polynomials=[0, 1,2,3,4,5,6,7], functions1 = [sin_cos_function_library(10)...])

# ╔═╡ 8eb49379-7b2a-4edb-924f-b6417450febf
lorenz_001_system = solve_constructed_system(Ξ_lorenz_001, [-8., 7, 27], (0., 100.); polynomials=[0, 1,2,3,4,5,6,7], functions1 = [sin_cos_function_library(10)...])

# ╔═╡ bfeee94b-b98d-4345-8b8a-efd9eda09c22
lorenz_01_system = solve_constructed_system(Ξ_lorenz_01, [-8., 7, 27], (0., 100.); polynomials=[0, 1,2,3,4,5,6,7]) # , functions1 = [sin_cos_function_library(10)...]

# ╔═╡ e00de518-1c38-452c-8b79-4adbfeecea53
lorenz_001_tv_system = solve_constructed_system(Ξ_lorenz_001_tv, [-8., 7, 27], (0., 100.); polynomials=[0, 1,2,3,4,5,6,7], functions1 = [sin_cos_function_library(10)...])

# ╔═╡ e4fe9390-6955-405c-86cc-cc9defb9ca39
lorenz_01_tv_system = solve_constructed_system(Ξ_lorenz_01_tv, [-8., 7, 27], (0., 1000.); polynomials=[0, 1,2,3,4,5,6,7], functions1 = [sin_cos_function_library(10)...])

# ╔═╡ c1ea010b-9549-43d7-8aa2-9b6a92c00d82
lorenz_001_solution = solve(lorenz_001_system, DP5())

# ╔═╡ 0814a304-c439-40e5-84a9-4e501aaf2326
# lorenz_01_solution = solve(lorenz_01_system, DP5())

# ╔═╡ 3d9a773b-39a9-41bd-b9a9-4c79a33992de
lorenz_001_tv_solution = solve(lorenz_001_tv_system, DP5())

# ╔═╡ b4a92457-6529-4ab9-872f-0bc2ba099ff4
lorenz_01_tv_solution = solve(lorenz_01_tv_system, DP5())

# ╔═╡ 7a8cc3aa-5b30-46b9-83db-3f47d0377c13
plot(lorenz_001_solution, idxs=(1,2,3), xlabel="x", ylabel="y", zlabel="z", camera=(40,30), legend = false)

# ╔═╡ f35be9d2-efea-4f0a-869d-25bd24f946a0
# plot(lorenz_01_solution, idxs=(1,2,3))

# ╔═╡ e0003431-ffef-429e-80c8-ba1c7f1774d2
plot(lorenz_001_tv_solution, idxs=(1,2,3), xlabel="x", ylabel="y", zlabel="z", camera=(40,30), legend = false)

# ╔═╡ d333670f-c9c6-461f-81d8-84f7bb974d85
plot(lorenz_01_tv_solution, idxs=(1,2,3), xlabel="x", ylabel="y", zlabel="z", camera=(40,30), legend = false)

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
		u0 = [2.,0]
		hopf_system = ODEProblem(hopf_normal!, u0, hopf_tspan, (μ,  hopf_ω, hopf_A))
		hopf_solution = solve(hopf_system, Rosenbrock23(), dt = hopf_dt, saveat=hopf_dt_saveat)

		# add "sensor noise"
		hopf_noisy = add_noise_to_sample(hopf_solution.u, hopf_η)

		# compute derivatives using tvdiff
		X_hopf, Xdot_hopf = constructX(hopf_noisy, TvDiffParams(; iterations = 5, α = 2, ε = 1e-3, perform_tvdiff = true, diff_kernel = "abs"), hopf_solution.t)
		
		# save
		push!(hopf_data, hcat(X_hopf, fill(μ, size(X_hopf)[1])))
		push!(hopf_derivatives, hcat(Xdot_hopf, fill(0, size(X_hopf)[1])))
		push!(hopf_solutions, hopf_solution)

		# don't do inner initial condition for case where we don't have limit cycle
		if μ < 0 
			continue
		end # if 

		u02 = [0.01, 0]

		# do same thing again with new initial condition 
		hopf_system2 = ODEProblem(hopf_normal!, u02, hopf_tspan, (μ,  hopf_ω, hopf_A))
		hopf_solution2 = solve(hopf_system2, Rosenbrock23(), dt = hopf_dt, saveat=hopf_dt_saveat)

		# add "sensor noise"
		hopf_noisy2 = add_noise_to_sample(hopf_solution2.u, hopf_η)

		# compute derivatives using tvdiff
		X_hopf2, Xdot_hopf2 = constructX(hopf_noisy2, TvDiffParams(; iterations = 5, α = 2, ε = 1e-3, perform_tvdiff = true, diff_kernel = "abs"), hopf_solution2.t)

		# save
		push!(hopf_data, hcat(X_hopf2, fill(μ, size(X_hopf2)[1])))
		push!(hopf_derivatives, hcat(Xdot_hopf2, fill(0, size(X_hopf2)[1])))

		push!(hopf_solutions, hopf_solution2)
	end # for μ


	X_hopf_combined = vcat(hopf_data...)
	Xdot_hopf_combined = vcat(hopf_derivatives...)
end

# ╔═╡ 3781cba0-2f81-471d-accc-de4d8546c11c
begin 
	hopf_polynomials = collect(0:7)
	hopf_functions1 = [] # [sin_cos_function_library(10)...]
	Θ_hopf = Θ(X_hopf_combined; polynomials=hopf_polynomials, functions1 = hopf_functions1)

end

# ╔═╡ 20d272ea-c701-4edf-a393-a48f514465fc
begin 
	Ξ_hopf = sparse_regression(Θ_hopf, Xdot_hopf_combined, 0.25, 10)
end

# ╔═╡ 261776eb-e90e-4471-8276-cb1b934cde28
print_constructed_system(Ξ_hopf; polynomials = hopf_polynomials, functions1 = hopf_functions1)

# ╔═╡ f6cb782c-934b-4342-b008-10be1fda2021
begin
	found_hopf_mu_neg015_1_system = solve_constructed_system(Ξ_hopf, [2, 0., -0.2], (0.0, 100.); polynomials=hopf_polynomials, functions1 = hopf_functions1)
	found_hopf_mu_neg01_1_system = solve_constructed_system(Ξ_hopf, [2, 0., -0.1], (0.0, 100.); polynomials=hopf_polynomials, functions1 = hopf_functions1)
	found_hopf_mu_neg005_1_system = solve_constructed_system(Ξ_hopf, [2, 0., -0.05], (0.0, 100.); polynomials=hopf_polynomials, functions1 = hopf_functions1)
	found_hopf_mu_005_1_system = solve_constructed_system(Ξ_hopf, [2, 0., 0.05], (0.0, 100.); polynomials=hopf_polynomials, functions1 = hopf_functions1)
	found_hopf_mu_005_2_system = solve_constructed_system(Ξ_hopf, [0.1, 0., 0.05], (0.0, 100.); polynomials=hopf_polynomials, functions1 = hopf_functions1)
	found_hopf_mu_01_1_system = solve_constructed_system(Ξ_hopf, [2, 0., 0.1], (0.0, 100.); polynomials=hopf_polynomials, functions1 = hopf_functions1)
	found_hopf_mu_01_2_system = solve_constructed_system(Ξ_hopf, [0.1, 0., 0.1], (0.0, 100.); polynomials=hopf_polynomials, functions1 = hopf_functions1)
	found_hopf_mu_022_1_system = solve_constructed_system(Ξ_hopf, [2, 0., 0.22], (0.0, 100.); polynomials=hopf_polynomials, functions1 = hopf_functions1)
	found_hopf_mu_022_2_system = solve_constructed_system(Ξ_hopf, [0.1, 0., 0.22], (0.0, 100.); polynomials=hopf_polynomials, functions1 = hopf_functions1)
	
	found_hopf_mu_neg015_1_solution = solve(found_hopf_mu_neg015_1_system, DP5())
	found_hopf_mu_neg01_1_solution = solve(found_hopf_mu_neg01_1_system, DP5())
	found_hopf_mu_neg005_1_solution = solve(found_hopf_mu_neg005_1_system, DP5())
	found_hopf_mu_005_1_solution = solve(found_hopf_mu_005_1_system, DP5())
	found_hopf_mu_005_2_solution = solve(found_hopf_mu_005_2_system, DP5())
	found_hopf_mu_01_1_solution = solve(found_hopf_mu_01_1_system, DP5())
	found_hopf_mu_01_2_solution = solve(found_hopf_mu_01_2_system, DP5())
	found_hopf_mu_022_1_solution = solve(found_hopf_mu_022_1_system, DP5())
	found_hopf_mu_022_2_solution = solve(found_hopf_mu_022_2_system, DP5())
end

# ╔═╡ 379bba48-ebc0-415d-bd54-25707b97bab1
begin
	plot(found_hopf_mu_neg015_1_solution, idxs=(3, 1,2), color = :blue)
	plot!(found_hopf_mu_neg01_1_solution, idxs=(3, 1,2), color = :blue)
	plot!(found_hopf_mu_neg005_1_solution, idxs=(3, 1,2), color = :blue)
	plot!(found_hopf_mu_005_1_solution, idxs=(3, 1,2), color = :blue)
	plot!(found_hopf_mu_005_2_solution, idxs=(3, 1,2), color = :red)
	plot!(found_hopf_mu_01_1_solution, idxs=(3, 1,2), color = :blue)
	plot!(found_hopf_mu_01_2_solution, idxs=(3, 1,2), color = :red)
	plot!(found_hopf_mu_022_1_solution, idxs=(3, 1,2), color = :blue)
	plot!(found_hopf_mu_022_2_solution, idxs=(3, 1,2), color = :red, zlabel = "x", ylabel = "y", xlabel = "mu", camera = (30,10), legend = false)

	
end

# ╔═╡ 2df97f42-ec7f-4c6f-9431-3b69a8cb8024
# begin
# 	test1 = ODEProblem(lorenz!, [-8., 7., 27.], (0.,100.), (10., 28., 8/3))
# 	saving = collect(0:0.01:100)
# 	res = solve(test1, DP5(), dt = 0.01, saveat=saving)

# 	res.u
# end

# ╔═╡ daa284d3-f9bc-45b7-9c53-df70ab60c204


# ╔═╡ Cell order:
# ╠═71998b4b-1bfc-4c3d-8c7e-5741de700162
# ╠═a0a26a0d-a7e9-4168-996b-51f44fd96dad
# ╠═0e6831be-f2e1-11ee-1ecc-1f97692f0c45
# ╠═aa29b69a-5750-4a5f-8e8c-fdc49d9a45b6
# ╠═1511b58f-0ec9-4dc3-90ec-d3169cf9189e
# ╟─61ccdc43-c12f-4b59-b746-316c82b342a5
# ╠═f2f7d603-4836-4010-95c7-10def9254590
# ╠═4cbdb8af-68ca-421f-853c-5bd8c4f4b464
# ╠═b2a4fb51-90a0-49af-9a3e-f0d21a54495d
# ╠═d3c346a0-8970-4e12-aaf8-5884a0ef435c
# ╠═d74cdb2f-7bfe-452c-b96b-1b9bf33c29f2
# ╠═69042648-7393-453d-bd59-ad2b345e1ada
# ╠═c948a893-a13c-45c2-a66b-8eff8ab85d68
# ╠═0d836d4e-d9d2-4fa1-9983-e541fb58ab04
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
# ╠═f35be9d2-efea-4f0a-869d-25bd24f946a0
# ╠═e0003431-ffef-429e-80c8-ba1c7f1774d2
# ╠═d333670f-c9c6-461f-81d8-84f7bb974d85
# ╠═768ff383-57d4-4203-8c25-4cba99c993bf
# ╠═3781cba0-2f81-471d-accc-de4d8546c11c
# ╠═20d272ea-c701-4edf-a393-a48f514465fc
# ╠═261776eb-e90e-4471-8276-cb1b934cde28
# ╠═f6cb782c-934b-4342-b008-10be1fda2021
# ╠═379bba48-ebc0-415d-bd54-25707b97bab1
# ╠═51281685-dea0-4c19-b6c6-893b265b051c
# ╠═2df97f42-ec7f-4c6f-9431-3b69a8cb8024
# ╠═daa284d3-f9bc-45b7-9c53-df70ab60c204
