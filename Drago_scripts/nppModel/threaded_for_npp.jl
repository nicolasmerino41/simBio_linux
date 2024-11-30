PC = "nicol"
using Pkg
using GLM
using CSV, DataFrames
using DifferentialEquations
using Distributions
using Random

# Define the Herbivores struct
mutable struct Herbivores
    m::Float64        # Mortality rate
    H0::Float64       # Characteristic density (H_i^0)
    H_init::Float64   # Initial abundance (H_i(0))
    g::Float64        # Growth rate (to be calculated)
    p_i::Float64      # Portion
end

# Outer constructor to accept keyword arguments
Herbivores(; m::Float64, H0::Float64, H_init::Float64, g::Float64=0.0, p_i::Float64) = Herbivores(m, H0, H_init, g, p_i)

# Function to create herbivores_list
function create_herbivores_list(num_herbivores::Int; m_mean::Float64, m_sd::Float64,
                                H0_mean::Float64, H0_sd::Float64,
                                H_init::Float64,
                                p_i_mean::Float64, p_i_sd::Float64)
    herbivores_list = Herbivores[]
    for i in 1:num_herbivores
        m = abs(rand(Normal(m_mean, m_sd)))          # Mortality rate
        H0 = abs(rand(Normal(H0_mean, H0_sd)))        # Characteristic density
        p_i = abs(rand(Normal(p_i_mean, p_i_sd)))
        push!(herbivores_list, Herbivores(m=m, H0=H0, H_init=H_init, p_i=p_i))
    end
    return herbivores_list
end

# Function to calculate growth rates based on NPP
function calculate_growth_rates(herbivores_list, NPP, mu)
    S_star = length(herbivores_list)
    # Calculate Fi for each herbivore
    F_list = [sp.H0 * sp.m for sp in herbivores_list]
    # Calculate the numerator of the competition term
    competition_numerator = 1 + mu * (S_star - 1)
    # Calculate gi for each herbivore
    for (i, sp) in enumerate(herbivores_list)
        Fi = F_list[i]
        # if approximate == true
        #     sp.g = sp.m * sqrt((competition_numerator / S_star) * (NPP / Fi)) # The approximation
        # elseif approximate == false  
        #     sp.g = sp.m * ((1+sqrt(1+((4*competition_numerator*NPP)/(S_star*Fi))))/2) # The exact
        # elseif isnothing(approximate)
            sp.g = sp.m * ((1+sqrt(1+((4*competition_numerator*sp.p_i*NPP)/(Fi))))/2) # The new
        # end
    end
end

# Predator struct definition
mutable struct Predator
    m::Float64        # Mortality rate
    a::Float64        # Attack rate
    h::Float64        # Handling time
    e::Float64        # Conversion efficiency
    P_init::Float64   # Initial abundance
    c::Float64        # Self-regulation coefficient
end

# Predator constructor
Predator(; m::Float64, a::Float64, h::Float64, e::Float64, P_init::Float64, c::Float64) = Predator(m, a, h, e, P_init, c)

# Function to create a list of predators with adjusted parameters
function create_predator_list(num_predators::Int; m_mean::Float64=0.1, m_sd::Float64=0.02,
                              a_mean::Float64=0.01, a_sd::Float64=0.0001,
                              h_mean::Float64=0.1, h_sd::Float64=0.01,
                              e_mean::Float64=0.1, e_sd::Float64=0.01,
                              P_init_mean::Float64=P_init_mean, P_init_sd::Float64=P_init_mean/10,
                              c_mean::Float64=0.1, c_sd::Float64=0.01)
    predator_list = Predator[]
    for _ in 1:num_predators
        m = rand(Normal(m_mean, m_sd))              # Mortality rate
        a = rand(Normal(a_mean, a_sd))            # Attack rate
        h = rand(Normal(h_mean, h_sd))              # Handling time
        e = rand(Normal(e_mean, e_sd))              # Conversion efficiency
        P_init = rand(Normal(P_init_mean, P_init_sd))  # Initial abundance
        c = rand(Normal(c_mean, c_sd))              # Self-regulation coefficient
        push!(predator_list, Predator(m=m, a=a, h=h, e=e, P_init=P_init, c=c))
    end
    return predator_list
end

# Function to generate the interaction matrix
function generate_interaction_matrix(num_predators::Int, num_prey::Int, connectivity::Float64)
    IM = zeros(Bool, num_predators, num_prey)
    for k in 1:num_predators
        for i in 1:num_prey
            IM[k, i] = rand() < connectivity  # Assign 1 with probability equal to connectivity
        end
    end
    return IM
end

# Ecosystem dynamics function with IM and self-regulation
function ecosystem_dynamics!(du, u, p, t)
    herbivores_list, beta_matrix, predator_list, IM = p
    S_star = length(herbivores_list)
    num_predators = length(predator_list)
    H = u[1:S_star]  # Herbivore densities
    P = u[S_star+1:end]  # Predator densities
    du_H = zeros(S_star)
    du_P = zeros(num_predators)

    # Herbivore dynamics
    for i in 1:S_star
        if H[i] > 0  # Only update if species is alive
            sp = herbivores_list[i]
            # Compute competition term
            competition = sum(beta_matrix[i, :] .* H) / sp.H0
            # Compute predation term
            predation = 0.0
            for k in 1:num_predators
                if IM[k, i] && P[k] > 0
                    pred = predator_list[k]
                    f_ki = (pred.a * H[i]) / (1 + pred.a * pred.h * H[i])  # Functional response
                    predation += P[k] * f_ki
                end
            end
            # Compute derivative for herbivores
            du_H[i] = H[i] * sp.m * ((sp.g / sp.m) - 1 - competition) - predation
        else
            du_H[i] = 0.0  # Keep derivative at zero if extinct
        end
    end

    # Predator dynamics with self-regulation
    for k in 1:num_predators
        if P[k] > 0  # Only update if species is alive
            pred = predator_list[k]
            ingestion = 0.0
            for i in 1:S_star
                if IM[k, i] && H[i] > 0
                    f_ki = (pred.a * H[i]) / (1 + pred.a * pred.h * H[i])  # Functional response
                    ingestion += f_ki
                end
            end
            # Compute derivative for predators with self-regulation
            du_P[k] = P[k] * (pred.e * ingestion - pred.m - pred.c * P[k])
        else
            du_P[k] = 0.0  # Keep derivative at zero if extinct
        end
    end

    # Assign derivatives
    du[1:S_star] = du_H
    du[S_star+1:end] = du_P
end

# Import SpinLock for thread-safety
using Base.Threads: SpinLock

# Initialize a SpinLock for file access
const file_lock = SpinLock()

# Function to write results safely
function write_results_safely(results_row, csv_filename)
    lock(file_lock)  # Lock the file for thread-safe access
    try
        if isfile(csv_filename)
            CSV.write(csv_filename, results_row, append=true)
        else
            CSV.write(csv_filename, results_row)
        end
    finally
        unlock(file_lock)  # Unlock the file
    end
end

# Updated run_simulation function with thread-safe CSV writing
function run_simulation(num_herbivores, num_predators, mu, NPP, m_mean_h, m_mean_p, connectivity, c_mean_p)
    
    H0_mean = NPP/num_herbivores
    H0_mean_aprox = H0_mean
    p_i_sd = 0.0
    m_sd_h = 0.01

    # Your simulation logic remains the same
    herbivores_list = create_herbivores_list(num_herbivores; m_mean=m_mean_h, m_sd=m_sd_h,
        H0_mean=NPP/num_herbivores, H0_sd=NPP/num_herbivores/0.1,
        H_init=NPP/num_herbivores, 
        p_i_mean=NPP/num_herbivores, p_i_sd=0.0
    )
    calculate_growth_rates(herbivores_list, NPP, mu)

    S_star = length(herbivores_list)
    beta_matrix = fill(mu, S_star, S_star)
    for i in 1:S_star
        beta_matrix[i, i] = 1.0
    end

    predator_list = create_predator_list(num_predators; 
        m_mean=m_mean_p, m_sd=0.02, 
        a_mean=0.01, a_sd=0.0001,
        h_mean=0.1, h_sd=0.01,
        e_mean=0.1, e_sd=0.01,
        P_init_mean= NPP/num_herbivores*0.1, P_init_sd=1.0, 
        c_mean=c_mean_p, c_sd=0.01
    )

    IM = generate_interaction_matrix(num_predators, S_star, connectivity)
    H_init_values = [sp.H_init for sp in herbivores_list]
    P_init_values = [pred.P_init for pred in predator_list]
    u_init = vcat(H_init_values, P_init_values)
    tspan = (0.0, 1000.0)

    p = (herbivores_list, beta_matrix, predator_list, IM)
    prob = ODEProblem(ecosystem_dynamics!, u_init, tspan, p)
    sol = solve(prob, Tsit5(); reltol=1e-6, abstol=1e-6)

    herbivore_biomass = sum(sol[1:S_star, end])
    predator_biomass = sum(sol[S_star+1:end, end])
    herbivore_data = sol[1:length(herbivores_list), :]  # Herbivore dynamics
    predator_data = sol[length(herbivores_list)+1:end, :]  # Predator dynamics
    total_biomass = herbivore_biomass + predator_biomass
    total_biomass_vs_npp = total_biomass / NPP
    num_extinct_herbivores = count(sol[1:S_star, end] .<= 1.0)
    num_extinct_predators = count(sol[S_star+1:end, end] .<= 1.0)
    are_there_extinctions = num_extinct_herbivores > 0 || num_extinct_predators > 0
    prop_of_sp_extinct = (num_extinct_herbivores + num_extinct_predators) / (num_herbivores + num_predators)

    fi_over_4 = sum([sp.H0 for sp in herbivores_list].*[sp.m for sp in herbivores_list]./4)
    sum_pi_NPP_larger_than_fi_over_4 = NPP > fi_over_4

    # Equation holding true?
    growth_rates = [sp.g for sp in herbivores_list]
    holding = (NPP/sum(growth_rates .* herbivore_data[:, end]) > 0.99) && (NPP/sum(growth_rates .* herbivore_data[:, end]) < 1.001)

    # Create results row
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
        p_i = 1/num_herbivores,
        num_extinct_herbivores=num_extinct_herbivores,
        num_extinct_predators=num_extinct_predators,
        herbivore_biomass=herbivore_biomass,
        predator_biomass=predator_biomass,
        are_there_extinctions=are_there_extinctions,
        prop_of_sp_extinct=prop_of_sp_extinct,
        total_biomass_vs_npp=total_biomass_vs_npp,
        holding=holding,
        fi_over_4=fi_over_4,
        NPP_vS_fi4=sum_pi_NPP_larger_than_fi_over_4
    )

    # Write results using the thread-safe function
    csv_filename = "resultados_threaded/simulation_results_30_11.csv"
    write_results_safely(results_row, csv_filename)

    # println("go")
    return 1.0
end

using IterTools

# Simulation parameters
num_herbivores_range = [1, 5, 10, 15]
num_predators_range = [0, 1, 5, 10]
mu_values = range(0.0, 1.0, length=16)
NPP_values = [1.0, 10.0, 100.0, 1000.0, 10000.0]
m_mean_h_values = [0.1, 0.2, 0.3]
m_mean_p_values = [0.1, 0.2, 0.3]
connectivity_values = [0.1, 0.25, 0.5, 0.75, 1.0]
c_mean_p_values = [0.0, 0.5, 1.0]

# Generate all parameter combinations
parameter_combinations = collect(IterTools.product(
    num_herbivores_range,
    num_predators_range,
    mu_values,
    NPP_values,
    m_mean_h_values,
    m_mean_p_values,
    connectivity_values,
    c_mean_p_values
))

# Parallelized loop over all combinations
Threads.@threads for params in parameter_combinations
    num_herbivores, num_predators, mu, NPP, m_mean_h, m_mean_p, connectivity, c_mean_p = params
    run_simulation(num_herbivores, num_predators, mu, NPP, m_mean_h, m_mean_p, connectivity, c_mean_p)
end


