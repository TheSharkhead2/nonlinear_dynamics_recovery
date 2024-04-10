using Random
using Distributions 

function add_noise_to_sample(data::Vector{Vector{Float64}}, η)
    # generate noise out of normal centered at 0 with maginutde η
    noise::Vector{Vector{Float64}} = η .* eachrow(rand(Normal(0,1), (size(data)[1], size(data[1])[1])))

    data + noise
end # function add_noise_to_sample
