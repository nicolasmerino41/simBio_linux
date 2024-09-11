using Pkg

# Packages
using CSV, DataFrames
using Distributions, NamedArrays, StaticArrays
using DynamicGrids, Dispersal
using DimensionalData, Rasters, Serialization, ArchGDAL
using OrderedCollections, StatsBase

# Setup code
include("/workdir/HerpsVsBirmmals.jl")
include("/workdir/efficient_setup.jl")
include("human_footprint.jl")
include("/workdir/New_metrics.jl")

# println("Break 1.jl")
self_regulation = 0.001

caca = deepcopy(iberian_interact_NA)
# println("Break 2.jl")
full_IM = turn_adj_into_inter(caca, sigma, epsilon)
full_IM = Matrix(full_IM)
# println("Break 3.jl")
outdisp = OutwardsDispersal{:birmmals, :birmmals}(;
    formulation=CustomKernel(alpha),
    distancemethod=AreaToArea(30),
    maskbehavior = Dispersal.CheckMaskEdges()
)

DA_birmmals_with_abundances = deepcopy(DA_birmmals)
# Iterate over rows and columns
for row in axes(DA_birmmals, 1), col in axes(DA_birmmals, 2)
    current = DA_birmmals[row, col]
    empty_birmmals = MyBirmmals(SVector{207, Float64}(fill(0.0, 207)))
    
    if current != empty_birmmals
        new_a = SVector{207, Float64}([current.a[i] != 0.0 ? initial_abundance : current.a[i] for i in 1:207])
        DA_birmmals_with_abundances[row, col] = MyBirmmals(new_a)
    end
end

DA_herps_with_abundances = deepcopy(DA_herps)
# Iterate over rows and columns
for row in axes(DA, 1), col in axes(DA, 2)
    current = DA[row, col]
    empty_struct = MyHerps(SVector{49, Float64}(fill(0.0, 49)))
    
    if current != empty_struct
        new_a = SVector{49, Float64}([current.a[i] != 0.0 ? initial_abundance : current.a[i] for i in 1:49])
        DA_herps_with_abundances[row, col] = MyHerps(new_a)
    end
end

# println("Break 8.jl")
#TODO At prevalence 0.277 or higher you get instability
pepe_state = (
    state = Matrix(DA_with_abundances),
    k_DA = Matrix(k_DA_hf),
    npp_DA = Matrix(npp_DA),
    state_richness = Matrix(DA_richness)
)

# println("Break 9.jl")
##### LAX NICHE #####
array_output = ResultOutput(
    pepe; tspan = 1:2000,
    mask = Matrix(DA_sum)
)

# println("Break 10.jl")
# println("full_IM is a ", typeof(full_IM))
# println("birmmals is a ", typeof(Matrix(DA_random_birmmals_with_abundances)))

@time p = sim!(array_output, Ruleset(biotic_rule_k, outdisp; boundary = Reflect()))
# println("Break 11.jl")

richness_eval = 1-richness_evaluation(p)
biomass_distribution = average_bbp(p)
average_shannon = average_shannon_index(p)

# println("Break 12.jl")
####################################################################
println("Initial_abundance: ", round(initial_abundance, digits = 2), ", Epsilon: ", round(epsilon, digits = 2), ", Sigma: ", round(sigma, digits = 2), ", Alpha: ", round(alpha, digits = 2))
println("The richness_eval was: ", richness_eval)
println("The average_shannon was: ", average_shannon)
println("The biomass_distribution was: ", biomass_distribution)
println("The sum of $epsilon + $sigma is: ", epsilon + sigma)
