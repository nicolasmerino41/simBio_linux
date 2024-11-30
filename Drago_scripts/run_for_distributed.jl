# Packages
using Distributed
using ClusterManagers

# **Add worker processes first**
addprocs(SlurmManager(16),exeflags=["--project","--threads=16"])

@everywhere begin

    using Base.Iterators:product  
    using CSV, DataFrames
    using Distributions, NamedArrays, StaticArrays
    using DynamicGrids, Dispersal
    using DimensionalData, Rasters, Serialization, ArchGDAL
    using OrderedCollections, StatsBase

    dir = pwd()
    num_species = 256
end

# Setup code
@everywhere include(joinpath(dir, "Scripts/HerpsVsBirmmals.jl"))

@everywhere include(joinpath(dir, "Scripts/kernels_for_drago.jl"))

@everywhere  include(joinpath(dir, "Scripts/efficient_setup_for_drago.jl"))

@everywhere  include(joinpath(dir, "Scripts/human_footprint_for_drago.jl"))

@everywhere  include(joinpath(dir, "Scripts/New_metrics_for_drago.jl"))

@everywhere include(joinpath(dir, "Scripts/Implicit_competition_for_herbivores.jl"))



# **Define variables and functions after adding workers**
@everywhere begin
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
end

# Function to save parameters, grid type (k_DA name), and metrics to CSV and append plots to the final PDFs
@everywhere function run_simulation(epsilon, alfa, sigma_comp, assymetry)
    sigmas = [0.0001, 0.001, 0.005, 0.01, 0.05, 0.07, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.0, 1.5, 2.0, 3.0]
    
    Threads.@threads for sigma in sigmas          

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
            pepe_state; tspan = 1:1000,
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
        
        # println(NaNs)
        # Step 2: Save the parameters, grid type, and metrics in a CSV
        results_row = DataFrame(
            sigma = sigma,
            epsilon = epsilon,
            alfa = alfa,
            k_DA_name = k_DA_name,
            sigma_comp = sigma_comp,
            assymetry = assymetry,
            avg_shannon = round(avg_shannon, digits = 2),
            avg_bbp = round(avg_bbp, digits = 2),
            richness_similarity = round(richness_sim, digits = 2),
            alive_predators = round(alive_preds, digits = 2),
            mean_tl = round(mean_tl, digits = 2),
            mean_n_of_species = round(mean_sp, digits = 2),
            NaNs = NaNs
        )
                
        # Append or create the CSV file
        csv_filename = "resultados_distribuidos/DirectSamplingResults.csv"
        if isfile(csv_filename)
            CSV.write(csv_filename, results_row, append = true)
        else
            CSV.write(csv_filename, results_row)
        end
        if NaNs == 0.0
            serialize("resultados_distribuidos/outputs/$sigma-$epsilon-$alfa-$sigma_comp-$assymetry.jls", p[end].state)
        end
        return 1.0	
    end   
end

# Simulation parameters
@everywhere begin
    #sigmas = [0.0001, 0.001, 0.005, 0.01, 0.05, 0.07, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.0, 1.5, 2.0, 3.0]
    epsilons = [0.01, 0.1, 0.5, 1.0, 2.0, 3.0, 5.0]
    alfa_values = [0.001, 0.01, 0.05, 0.1, 0.3, 0.6, 0.9, 1.1]
    sigma_comp_values = [0.001, 0.1, 0.5, 1.0, 1.5, 2.0]
    assymetry_values = [0.0, 0.33, 0.66, 1.0]

    k_DA_list = [k_DA.DA_multiplicative, k_DA.DA_additive, k_DA.DA_min, k_DA.DA_harmonic, k_DA.DA_geometric]
    k_DA_names = ["multiplicative", "additive", "min", "harmonic", "geometric"]
    positions = [1, 2, 3, 4, 5]
    conjunto = collect(product(epsilons, alfa_values, sigma_comp_values, assymetry_values))
    #iteracciones = collect(product(sigmas, epsilons, alfa_values, sigma_comp_values, assymetry_values))
end

# Now you can use pmap
result = pmap(args -> run_simulation(args...), conjunto)

