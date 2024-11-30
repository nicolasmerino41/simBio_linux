# Define the Herbivore struct
mutable struct Herbivore
    m::Float64        # Mortality rate
    H0::Float64       # Characteristic density (H_i^0)
    H_init::Float64   # Initial abundance (H_i(0))
    g::Float64        # Growth rate (to be calculated)
end

# Constructor for Herbivore
Herbivore(; m::Float64, H0::Float64, H_init::Float64, g::Float64=0.0) = Herbivore(m, H0, H_init, g)

# Function to create herbivores list
function create_herbivores_list(num_herbivores::Int; m_mean::Float64=0.1, m_sd::Float64=0.02,
    H0_mean::Float64, H0_sd::Float64, H_init_mean::Float64=5.0, H_init_sd::Float64=1.0)
    herbivores_list = Herbivore[]
    for i in 1:num_herbivores
        m = abs(rand(Normal(m_mean, m_sd)))           # Mortality rate
        H0 = abs(rand(Normal(H0_mean, H0_sd)))        # Characteristic density
        H_init = H0  # Initial abundance set to characteristic density
        push!(herbivores_list, Herbivore(m=m, H0=H0, H_init=H_init))
    end
    return herbivores_list
end

# Function to calculate growth rates based on NPP
function calculate_growth_rates(herbivore_list, NPP, mu)
    S_star = length(herbivore_list)
    F_list = [sp.H0 * sp.m for sp in herbivore_list]
    competition_numerator = 1 + mu * (S_star - 1)
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

# Function to create a list of predators
function create_predator_list(num_predators::Int; m_mean::Float64, m_sd::Float64, a_mean::Float64, a_sd::Float64,
    h_mean::Float64, h_sd::Float64, e_mean::Float64, e_sd::Float64, P_init_mean::Float64, P_init_sd::Float64,
    c_mean::Float64, c_sd::Float64)
    predator_list = Predator[]
    for _ in 1:num_predators
        m = abs(rand(Normal(m_mean, m_sd)))             # Mortality rate
        a = abs(rand(Normal(a_mean, a_sd)))             # Attack rate
        h = abs(rand(Normal(h_mean, h_sd)))             # Handling time
        e = abs(rand(Normal(e_mean, e_sd)))             # Conversion efficiency
        P_init = abs(rand(Normal(P_init_mean, P_init_sd)))  # Initial abundance
        c = abs(rand(Normal(c_mean, c_sd)))             # Self-regulation coefficient
        push!(predator_list, Predator(m=m, a=a, h=h, e=e, P_init=P_init, c=c))
    end
    return predator_list
end

# Function to generate the interaction matrix
function generate_interaction_matrix(num_predators::Int, num_prey::Int, connectivity::Float64)
    IM = zeros(Bool, num_predators, num_prey)
    for k in 1:num_predators
        for i in 1:num_prey
            IM[k, i] = rand() < connectivity
        end
    end
    return IM
end

# Ecosystem dynamics function
function ecosystem_dynamics!(du, u, p, t)
    herbivore_list, beta_matrix, predator_list, IM = p
    S_star = length(herbivore_list)
    num_predators = length(predator_list)
    H = u[1:S_star]         # Herbivore densities
    P = u[S_star+1:end]     # Predator densities
    du_H = zeros(S_star)
    du_P = zeros(num_predators)

    # Herbivore dynamics
    for i in 1:S_star
        sp = herbivore_list[i]
        competition = sum(beta_matrix[i, :] .* H) / sp.H0
        predation = 0.0
        for k in 1:num_predators
            if IM[k, i]
                pred = predator_list[k]
                f_ki = (pred.a * H[i]) / (1 + pred.a * pred.h * H[i])
                predation += P[k] * f_ki
            end
        end
        du_H[i] = H[i] * sp.m * ((sp.g / sp.m) - 1 - competition) - predation
    end

    # Predator dynamics with self-regulation
    for k in 1:num_predators
        pred = predator_list[k]
        ingestion = 0.0
        for i in 1:S_star
            if IM[k, i]
                f_ki = (pred.a * H[i]) / (1 + pred.a * pred.h * H[i])
                ingestion += f_ki
            end
        end
        du_P[k] = P[k] * (pred.e * ingestion - pred.m - pred.c * P[k])
    end

    du[1:S_star] = du_H
    du[S_star+1:end] = du_P
end
