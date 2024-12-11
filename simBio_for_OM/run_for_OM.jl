num_species = 256
using Pkg
dir = ""
# Packages
using CSV, DataFrames
using Distributions, NamedArrays, StaticArrays
using DynamicGrids, Dispersal
using DimensionalData, Rasters, Serialization, ArchGDAL
using OrderedCollections, StatsBase

# Setup code
include(joinpath(dir, "HerpsVsBirmmals.jl"))
include(joinpath(dir, "kernels_for_OM.jl"))
include(joinpath(dir, "efficient_setup_for_OM.jl"))
include(joinpath(dir, "human_footprint_for_OM.jl"))
include(joinpath(dir, "New_metrics_for_OM.jl"))
include(joinpath(dir, "Implicit_competition_for_herbivores_for_OM.jl"))

initial_abundance = 0.41
DA_with_abundances = deepcopy(DA)
for row in axes(DA, 1), col in axes(DA, 2)
    if !iszero(DA[row, col]) 
        new_a = SVector{256, Float64}([DA[row, col].a[i] != 0.0 ? initial_abundance : DA[row, col].a[i] for i in 1:256])
        DA_with_abundances[row, col] = MyStructs256(new_a)
    end
end

const DA_with_abundances_const = deepcopy(DA_with_abundances)
const iberian_interact_NA_const = deepcopy(iberian_interact_NA)
# Function to save parameters, grid type (k_DA name), and metrics to CSV and append plots to the final PDFs
function run_simulation(sigma, epsilon, alfa, sigma_comp, assymetry)

    k_DA_name = "k_DA_hf_multiplicative"
    position = 1
    
    # Get the corresponding k_DA grid from the list
    k_DA = k_DA_hf_multiplicative

    DA_with_abundances = deepcopy(DA_with_abundances_const)
    pepe_state = (
        state = Matrix(DA_with_abundances),
        k_DA = Matrix(k_DA_hf_additive),
        npp_DA = Matrix(npp_DA),
        state_richness = Matrix(DA_richness)
    )

    ##### Matrices #####
    self_regulation = 1.0
    s = sigma
    e = epsilon
    beta = 3.0
    # Trophic
    caca = deepcopy(iberian_interact_NA_const)
    full_IM = Matrix(turn_adj_into_inter(caca, s, e, self_regulation, beta))
    # Competition
    s_comp = sigma_comp
    ass = assymetry
    competition_NA = deepcopy(iberian_interact_NA)
    competition_NA .= 0.0
    for i in names(competition_NA, 1), j in names(competition_NA, 2)
        if i in herbivore_names && j in herbivore_names
            competition_NA[i, j] = 1.0
        end 
    end
    full_comp = turn_comp_into_inter(competition_NA, s_comp, ass)

    a = alfa

    function GLV(state::MyStructs256, k_DA::MyStructs256)
        return MyStructs256(
            SVector{256, Float64}(
                state.a + (state.a .* (k_DA.a - state.a) + ((full_IM * state.a) .* state.a) + ((full_comp * state.a) .* state.a)) 
            )
        )
    end

    biotic_GLV = Cell{Tuple{:state, :k_DA}, :state}() do data, (state, k_DA), I
        # if any(isinf, state.a) || any(isnan, state.a)
        #     @error "state has NA values"
        #     println(I)
        # end
        return MyStructs256(SVector{256, Float64}(max.(0.0, GLV(state, k_DA).a)))
    end

    outdisp = OutwardsDispersal{:state, :state}(;
        formulation=CustomKernel(a),
        distancemethod=AreaToArea(30),
        maskbehavior = Dispersal.CheckMaskEdges()
    );

    println("sigma  = ", sigma, " epsilon = ", epsilon, " alfa = ", alfa, " sigma_comp = ", sigma_comp, " assymetry = ", assymetry)

    # Run the simulation
    array_output = ResultOutput(
        pepe_state; tspan = 1:500,
        mask = Matrix(DA_sum)
    )
    # println("output done")
    p = sim!(array_output, Ruleset(biotic_GLV, outdisp; boundary = Reflect()))
    # println("simulation done")
    # Step 1: Compute metrics from the last timestep of the simulation (p[end])
    avg_shannon = average_shannon_index(p, 1; modified = true)
    # println("shannon done")
    avg_bbp = average_bbp(p, position; modified = true)
    # println("bbp done")
    richness_sim = richness_similarity(p, position; modified = true)
    # println("richness done")
    alive_preds = alive_predators(p, position; modified = true)
    # println("predators done")
    mean_tl = calculate_mean_tl(p, position; modified = true)
    # println("meantl done")
    mean_sp = mean_n_of_species(p, position; modified = true)
    
    final_state = p[end].state
    NaNs = any(i -> any(isnan, final_state[idx[i][1], idx[i][2]].a), 1:length(idx)) ? 1.0 : 0.0

    return (
        avg_shannon = avg_shannon,
        avg_bbp = avg_bbp,
        richness_similarity = richness_sim,
        alive_predators = alive_preds,
        mean_tl = mean_tl,
        mean_n_of_species = mean_sp,
        # Threads = Threads.nthreads(),
        NaNs = NaNs
    )
end

avg_shannon, avg_bbp, richness_sim, alive_preds, mean_tl, mean_sp, NaNs = run_simulation(sigma, epsilon, alfa, sigma_comp, assymetry)

####################################################################
println("Epsilon: ", round(epsilon, digits = 2), ", Sigma: ", round(sigma, digits = 3), ", Alfa: ", round(alfa, digits = 2), ", Sigma_comp: ", round(sigma_comp, digits = 2), ", Assymetry: ", round(assymetry, digits = 2))
# println("The richness_eval was: ", round(richness_eval, digits = 2))
println("The average_shannon was: ", round(avg_shannon, digits = 2))
println("The biomass_distribution was: ", round(avg_bbp, digits = 2))
# println("MAE was: ", round(mae, digits = 2))
