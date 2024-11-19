# Required Packages
using DifferentialEquations
using Distributions
using Random

# OpenMOLE-provided variables (assumed to be defined)
# mu = 0.2                 # Competition coefficient (0 to 1)
# NPP = 100.0               # Net Primary Productivity (10 to 10000)
num_predators = Int(round(num_pred, digits = 0))      # Number of predator species
num_herbivores = Int(round(num_herb, digits = 0))     # Number of herbivore species
# H0_mean = 10.0            # Mean characteristic density of herbivores
# # Optionally, you can define 'seed' for reproducibility
# c_mean_p = 0.1
# # Derived parameters
H0_mean_aprox = H0_mean  # H0_mean provided by OpenMOLE
# connectivity = 1.0       # You can also make this an input if needed

# Define the Herbivore struct
mutable struct Herbivore
    m::Float64        # Mortality rate
    H0::Float64       # Characteristic density (H_i^0)
    H_init::Float64   # Initial abundance (H_i(0))
    g::Float64        # Growth rate (to be calculated)
end

# Constructor for Herbivore
Herbivore(; m::Float64, H0::Float64, H_init::Float64, g::Float64=0.0) = Herbivore(; m, H0, H_init, g)

# Function to create herbivores_list
function create_herbivores_list(num_herbivores::Int; m_mean::Float64=0.1, m_sd::Float64=0.02,
    H0_mean::Float64=H0_mean_aprox, H0_sd::Float64=H0_mean_aprox/10,
    H_init_mean::Float64=5.0, H_init_sd::Float64=1.0)
    herbivores_list = []
    for i in 1:num_herbivores
    m = rand(Normal(m_mean, m_sd))          # Mortality rate
    H0 = rand(Normal(H0_mean, H0_sd))        # Characteristic density
    H_init = H0  # Initial abundance set to characteristic density
    push!(herbivores_list, Herbivore(m=m, H0=H0, H_init=H_init))
    end
    return herbivores_list
end

# Function to calculate growth rates based on NPP
function calculate_growth_rates(herbivore_list, NPP, mu)
    S_star = length(herbivore_list)
    # Calculate Fi for each herbivore
    F_list = [sp.H0 * sp.m for sp in herbivore_list]
    # Calculate the numerator of the competition term
    competition_numerator = 1 + mu * (S_star - 1)
    # Calculate gi for each herbivore
    for (i, sp) in enumerate(herbivore_list)
        Fi = F_list[i]
        sp.g = sp.m * sqrt((competition_numerator / S_star) * (NPP / Fi))
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
    P_init_mean::Float64=5.0, P_init_sd::Float64=1.0,
    c_mean::Float64=0.1, c_sd::Float64=0.01)
predator_list = []
for _ in 1:num_predators
m = rand(Normal(m_mean, m_sd))              # Mortality rate
a = rand(Normal(a_mean, a_sd))            # Attack rate
h = rand(Normal(h_mean, h_sd))              # Handling time
e = rand(Normal(e_mean, e_sd))              # Conversion efficiency
P_init = rand(Normal(P_init_mean, P_init_mean/10))  # Initial abundance
c = rand(Normal(c_mean, c_sd))  # Define mean and sd for c
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

# Ecosystem dynamics function with IM
function ecosystem_dynamics!(du, u, p, t)
    herbivore_list, beta_matrix, predator_list, IM = p
    S_star = length(herbivore_list)
    num_predators = length(predator_list)
    H = u[1:S_star]  # Herbivore densities
    P = u[S_star+1:end]  # Predator densities
    du_H = zeros(S_star)
    du_P = zeros(num_predators)

    # Herbivore dynamics
    for i in 1:S_star
        sp = herbivore_list[i]
        # Compute competition term
        competition = sum(beta_matrix[i, :] .* H) / sp.H0
        # Compute predation term
        predation = 0.0
        for k in 1:num_predators
            if IM[k, i]
                pred = predator_list[k]
                f_ki = (pred.a * H[i]) / (1 + pred.a * pred.h * H[i])  # Functional response
                predation += P[k] * f_ki
            end
        end
        # Compute derivative for herbivores
        du_H[i] = H[i] * sp.m * ((sp.g / sp.m) - 1 - competition) - predation
    end

     # Predator dynamics with self-regulation
     for k in 1:num_predators
        pred = predator_list[k]
        ingestion = 0.0
        for i in 1:S_star
           if IM[k, i]
              f_ki = (pred.a * H[i]) / (1 + pred.a * pred.h * H[i])  # Functional response
              ingestion += f_ki
          end
        end
        # Self-regulation term
        du_P[k] = P[k] * (pred.e * ingestion - pred.m - pred.c * P[k])
    end

    # Assign derivatives
    du[1:S_star] = du_H
    du[S_star+1:end] = du_P
end

# Main simulation function
function run_simulation()
    # Create herbivore list and calculate growth rates
    herbivore_list = create_herbivores_list(num_herbivores; m_mean=m_mean_h, H0_mean=H0_mean_aprox)
    calculate_growth_rates(herbivore_list, NPP, mu)

    # Create beta_matrix
    S_star = length(herbivore_list)
    beta_matrix = fill(mu, S_star, S_star)
    for i in 1:S_star
        beta_matrix[i, i] = 1.0  # Self-competition
    end

    # Create predator list
    predator_list = create_predator_list(
        num_predators; 
        m_mean=m_mean_p, a_mean=a_mean_p, 
        h_mean=h_mean_p, e_mean=e_mean_p,
        c_mean=c_mean_p
    )

    # Generate the interaction matrix
    IM = generate_interaction_matrix(num_predators, S_star, connectivity)

    # Initial conditions for herbivores and predators
    H_init_values = [sp.H_init for sp in herbivore_list]
    P_init_values = [pred.P_init for pred in predator_list]

    # Combined initial conditions for the system
    u_init = vcat(H_init_values, P_init_values)

    # Define the time span for simulation
    tspan = (0.0, 200.0)

    # Define the ODE problem
    p = (herbivore_list, beta_matrix, predator_list, IM)
    prob = ODEProblem(ecosystem_dynamics!, u_init, tspan, p)

    # Solve the ODE
    sol = solve(prob, DifferentialEquations.Tsit5(); reltol=1e-6, abstol=1e-6)

    # Extract herbivore data
    H_array = sol[1:S_star, :]

    # Calculate total biomass at the end
    total_biomass = sum(H_array[:, end])

    herbivore_data = sol[1:length(herbivore_list), :]  # Herbivore dynamics
    predator_data = sol[length(herbivore_list)+1:end, :]  # Predator dynamics
    herbivore_biomass = sum(herbivore_data[:, end])
    predator_biomass = sum(predator_data[:, end])

    if true #any(herbivore_data[:, end] .<= 1.0) 
        num_extinct_herbivores = count(herbivore_data[:, end] .<= 1.0)
        # println(num_extinct_herbivores, " herbivore(s) went extinct.")
    else
        num_extinct_herbivores = 0
    end
    if true #any(predator_data[:, end] .<= 1.0)
        num_extinct_predators = count(predator_data[:, end] .<= 1.0)
        # println(num_extinct_predators, " predator(s) went extinct.")
    else
        num_extinct_predators = 0
    end
    are_there_extinctions = num_extinct_herbivores > 0 || num_extinct_predators > 0
    are_there_extinctions = Int(are_there_extinctions)

    prop_of_sp_extinct = (num_extinct_herbivores + num_extinct_predators) / (num_herbivores + num_predators)
    # Return the total biomass
    return total_biomass/NPP, num_extinct_herbivores, num_extinct_predators, herbivore_biomass, predator_biomass, are_there_extinctions, prop_of_sp_extinct
end

# Run the simulation and obtain the output
total_biomass, num_extinct_herbivores, num_extinct_predators, herbivore_biomass, predator_biomass, are_there_extinctions, prop_of_sp_extinct = run_simulation()
# Output the result (OpenMOLE will capture the variable 'total_biomass')
println("total_biomass = ", total_biomass)
println("herbivore_biomass = ", herbivore_biomass)
println("predator_biomass = ", predator_biomass)
println("pred/herb ratio = ", predator_biomass/herbivore_biomass)
println("num_extinct_herbivores = $num_extinct_herbivores/$num_herbivores")
println("num_extinct_predators = $num_extinct_predators/$num_predators")
println("are_there_extinctions = ", are_there_extinctions)
