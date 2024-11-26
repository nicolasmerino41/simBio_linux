# Packages
using Distributed
using ClusterManagers

# Add worker processes
addprocs(SlurmManager(16), exeflags=["--project", "--threads=16"])

@everywhere begin
    using CSV, DataFrames
    using DifferentialEquations
    using Distributions
    using Random
    dir = pwd()
    include(joinpath(dir, "Ecosystem_dynamics!.jl"))
end

# Define the iteration function
@everywhere function run_simulation(num_predators, mu, NPP, m_mean_h, m_mean_p, H0_mean, connectivity, c_mean_p)
    
    num_herbivoress = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
        
    Threads.@threads for num_herbivores in num_herbivoress
    
    herbivore_list = create_herbivores_list(num_herbivores; m_mean=m_mean_h, H0_mean=H0_mean, H0_sd=H0_mean/10)
    calculate_growth_rates(herbivore_list, NPP, mu)

    S_star = length(herbivore_list)
    beta_matrix = fill(mu, S_star, S_star)
    for i in 1:S_star
        beta_matrix[i, i] = 1.0
    end

    predator_list = create_predator_list(num_predators; m_mean=m_mean_p, m_sd=0.02, a_mean=0.01, a_sd=0.0001,
        h_mean=0.1, h_sd=0.01, e_mean=0.1, e_sd=0.01, P_init_mean=5.0, P_init_sd=1.0, c_mean=c_mean_p, c_sd=0.01)

    IM = generate_interaction_matrix(num_predators, S_star, connectivity)
    H_init_values = [sp.H_init for sp in herbivore_list]
    P_init_values = [pred.P_init for pred in predator_list]
    u_init = vcat(H_init_values, P_init_values)
    tspan = (0.0, 200.0)

    p = (herbivore_list, beta_matrix, predator_list, IM)
    prob = ODEProblem(ecosystem_dynamics!, u_init, tspan, p)
    sol = solve(prob, Tsit5(); reltol=1e-6, abstol=1e-6)

    herbivore_biomass = sum(sol[1:S_star, end])
    predator_biomass = sum(sol[S_star+1:end, end])
    total_biomass = herbivore_biomass + predator_biomass
    total_biomass_vs_npp = total_biomass/NPP
    num_extinct_herbivores = count(sol[1:S_star, end] .<= 1.0)
    num_extinct_predators = count(sol[S_star+1:end, end] .<= 1.0)

    prop_of_sp_extinct = (num_extinct_herbivores + num_extinct_predators) / (num_herbivores + num_predators)

    # Save results
        results_row = DataFrame(
            num_herbivores=num_herbivores,
            num_predators=num_predators,
            mu=mu,
            NPP=NPP,
            m_mean_h=m_mean_h,
            m_mean_p=m_mean_p,
            H0_mean=H0_mean,
            connectivity=connectivity,
            c_mean_p=c_mean_p,
            total_biomass=total_biomass,
            num_extinct_herbivores=num_extinct_herbivores,
            num_extinct_predators=num_extinct_predators,
            herbivore_biomass=herbivore_biomass,
            predator_biomass=predator_biomass,
            are_there_extinctions=are_there_extinctions,
            prop_of_sp_extinct=prop_of_sp_extinct,
            total_biomass_vs_npp=total_biomass_vs_npp
        )

        csv_filename = "resultados_distribuidos/EcosystemResults.csv"
        if isfile(csv_filename)
            CSV.write(csv_filename, results_row, append=true)
        else
            CSV.write(csv_filename, results_row)
        end
        return 1.0
    end
end

# Define parameter ranges
@everywhere begin
    # num_herbivores_range = 1:16
    num_predators_range = 1:10
    mu_values = range(0.0, 1.0, length=10)
    NPP_values = range(10.0, 10000.0, length=10)
    m_mean_h_values = range(0.001, 0.6, length=10)
    m_mean_p_values = range(0.001, 0.6, length=10)
    H0_mean_values = range(1.0, 1000.0, length=10)
    connectivity_values = range(0.0, 1.0, length=10)
    c_mean_p_values = range(0.0, 1.0, length=10)

    conjunto = collect(product(num_predators_range, mu_values, NPP_values, m_mean_h_values, m_mean_p_values, H0_mean_values, connectivity_values, c_mean_p_values))   
end

# Now you can use pmap
result = pmap(args -> run_simulation(args...), conjunto)


