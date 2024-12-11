using Statistics, Distributions, LinearAlgebra
using DifferentialEquations, DiffEqCallbacks

#############################
# Parameters
#############################
begin
    # S = 5   # Number of herbivore species
    # R = 2   # Number of predator species
    
    exponent_abundance = 0.0  # exponent for SAD (power law)
    # H0_mean = 100.0
    # M_mean = 0.1
    H0_sd = H0_mean/10
    M_sd = M_mean/10
    NPP = 1000.0
    
    # mu = 0.0                 # average herbivore-herbivore interaction strength
    mu_pred = 0.0           # average predator-predator interaction strength
    mu_predation = 0.1       # average herbivore-predator interaction strength
    asymmetry_competition = 1.0
    asymmetry_predators = 0.7
    asymmetry_predation = 1.0
    # epsilon = 0.01           # assimilation efficiency of predators
    
    # connectivity_hp = 1.0    # Herbivore-Predator connectivity
    connectivity_pp = 0.0    # Predator-Predator connectivity
    
    EXTINCTION_THRESHOLD = 1e-6
    final_time = 500.0
    
    #############################
    # Functions and Computations
    #############################
    function generate_hat_Hi(S::Int, exponent::Float64)
        ranks = collect(1:S)
        abundances = ranks .^ (-exponent)
        abundances ./= sum(abundances)
        return abundances
    end
    
    hat_H = generate_hat_Hi(S, exponent_abundance)
    
    H_i0 = [abs(rand(Normal(H0_mean, H0_sd))) for _ in 1:S]
    m_i = [abs(rand(Normal(M_mean, M_sd))) for _ in 1:S]
    barH = mean(H_i0)
    barM = mean(m_i)
    
    h_i = [H0/(S*barH) for H0 in H_i0]
    q_i = [mm/(S*barM) for mm in m_i]
    
    h_sum = sum(h_i)
    q_sum = sum(q_i)
    h_i = h_i ./ h_sum
    q_i = q_i ./ q_sum
    
    function generate_competition_matrix(S, mu, asymmetry)
        V_symmetric = ones(S, S)
        for i in 1:S, j in 1:S
            if i == j
                V_symmetric[i, j] = 1.0
            else
                V_symmetric[i, j] = mu
            end
        end
    
        V_random = ones(S, S)
        for i in 1:S, j in 1:S
            if i != j
                V_random[i, j] = rand()
            else
                V_random[i, j] = 1.0
            end
        end
    
        off_diag_values = [V_random[i,j] for i in 1:S, j in 1:S if i != j]
        current_mean = mean(off_diag_values)
        desired_mean = 1/(1+mu)
        scaling_factor = desired_mean / current_mean
    
        for i in 1:S, j in 1:S
            if i != j
                V_random[i,j] *= scaling_factor
            end
        end
    
        V = asymmetry * V_symmetric + (1.0 - asymmetry)*V_random
        return V
    end
    
    V = generate_competition_matrix(S, mu, asymmetry_competition)
    
    function generate_interaction_matrix_pred(R, mu_pred, asymmetry, connectivity)
        A_symmetric = fill(mu_pred, R, R)
        for α in 1:R
            A_symmetric[α, α] = -1.0
        end
    
        A_random = Matrix{Float64}(undef, R, R)
        for α in 1:R, β in 1:R
            if α == β
                A_random[α, β] = -1.0
            else
                if rand() < connectivity
                    A_random[α, β] = rand()
                else
                    A_random[α, β] = 0.0
                end
            end
        end
    
        off_vals = [A_random[α,β] for α in 1:R, β in 1:R if α != β]
        if !isempty(off_vals)
            current_mean = mean(off_vals)
            if current_mean != 0.0
                scaling_factor = mu_pred / current_mean
                for α in 1:R, β in 1:R
                    if α != β
                        A_random[α, β] *= scaling_factor
                    end
                end
            end
        end
    
        A = asymmetry * A_symmetric .+ (1.0 - asymmetry)*A_random
        return A
    end
    
    A = generate_interaction_matrix_pred(R, mu_pred, asymmetry_predators, connectivity_pp)
    
    function generate_predator_herbivore_matrix(rows::Int, cols::Int, connectivity::Float64, mu::Float64, asymmetry::Float64)
        @assert 0.0 ≤ connectivity ≤ 1.0
        @assert 0.0 ≤ asymmetry ≤ 1.0
        @assert mu ≥ 0
    
        mat = zeros(rows, cols)
        dist = Uniform(0, 2*mu)
        for i in 1:rows
            for j in 1:cols
                if rand() < connectivity
                    random_part = rand(dist)
                    value = asymmetry*mu + (1-asymmetry)*random_part
                    mat[i, j] = value
                end
            end
        end
        return mat
    end
    
    a_matrix = generate_predator_herbivore_matrix(S, R, connectivity_hp, mu_predation, asymmetry_predation)
    
    A_inv = inv(A)
    
    C = zeros(S,S)
    G = zeros(S)
    
    m_alpha = [M_mean for _ in 1:R]
    
    for i in 1:S
        for j in 1:S
            val = 0.0
            for α in 1:R
                for β in 1:R
                    val += epsilon * a_matrix[i, α]*A_inv[α, β]*a_matrix[j, β]
                end
            end
            C[i,j] = val
        end
    end
    
    for i in 1:S
        val = 0.0
        for β in 1:R
            for α in 1:R
                val += a_matrix[i, α]*A_inv[α, β]*m_alpha[β]
            end
        end
        G[i] = val
    end
    
    V_inv = inv(V)
    mu_matrix = V_inv - I
    
    M_modified = copy(mu_matrix)
    for i in 1:S, j in 1:S
        M_modified[i,j] = mu_matrix[i,j] + (C[i,j]*H_i0[i]/m_i[i])
    end
    
    hat_p = similar(hat_H)
    for i in 1:S
        interaction_sum = 0.0
        for j in 1:S
            interaction_sum += M_modified[i,j]*hat_H[j]
        end
        hat_p[i] = (hat_H[i] + interaction_sum)/h_i[i]
    end
    
    p = hat_p ./ sum(hat_p)
    
    IplusM = I + M_modified
    IM_inv = inv(IplusM)
    
    A_vec = zeros(S)
    B_vec = zeros(S)
    for i in 1:S
        A_val = 0.0
        B_val = 0.0
        for j in 1:S
            A_val += IM_inv[i,j]*p[j]
            B_val += IM_inv[i,j]*(G[j]/m_i[j] -1)
        end
        A_vec[i] = A_val
        B_vec[i] = B_val
    end
    
    A_coef = sum(p[i]*m_i[i]*H_i0[i]*A_vec[i] for i in 1:S)
    B_coef = sum(p[i]*m_i[i]*H_i0[i]*B_vec[i] + G[i]*H_i0[i]*A_vec[i] for i in 1:S)
    C_coef = sum(G[i]*H_i0[i]*B_vec[i] for i in 1:S) - NPP
    
    discriminant = B_coef^2 - 4*A_coef*C_coef
    
    ###################################
    # Modified behavior for x solution
    ###################################
    if discriminant < 0
        # No real solution for x: convergent=0, assign arbitrary values and skip dynamics
        convergent = 0
        # Arbitrary values for outputs:
        all_survived = false
        proportion_of_survived = 0.0
        herb_pred_bodymass_ratio = 0.0
        total_biomass = 0.0
        system_true = false
    else
        x_candidates = [(-B_coef + sqrt(discriminant))/(2*A_coef), (-B_coef - sqrt(discriminant))/(2*A_coef)]
        x_candidates = filter(x->x>0, x_candidates)
        if isempty(x_candidates)
            # No positive x: convergent=0
            convergent = 0
            # Arbitrary values:
            all_survived = false
            proportion_of_survived = 0.0
            herb_pred_bodymass_ratio = 0.0
            total_biomass = 0.0
            system_true = false
        else
            # Valid x found, proceed
            convergent = 1
            x_final = maximum(x_candidates)
    
            g_i = [x_final * p[i] * m_i[i] for i in 1:S]
    
            # Initial conditions
            H_init = [H_i0[i] for i in 1:S]
            P_init = fill(10.0, R)
            u0 = vcat(H_init, P_init)
    
            tspan = (0.0, final_time)
    
            params = (S, R, H_i0, m_i, g_i, G, M_modified, a_matrix, A, epsilon, m_alpha)
    
            function ecosystem_dynamics!(du, u, p, t)
                S, R, H_i0, m_i, g_i, G, M, a_matrix, A, epsilon, m_alpha = p
    
                H = @view u[1:S]
                P = @view u[S+1:S+R]
    
                duH = zeros(S)
                duP = zeros(R)
    
                for i in 1:S
                    if H[i] > 0.0
                        m_ii = m_i[i]
                        numerator = (g_i[i] + G[i])/m_ii - 1.0
                        interaction_sum = H[i]
                        for j in 1:S
                            interaction_sum += M[i,j]*H[j]
                        end
                        duH[i] = H[i]*m_ii*(numerator - interaction_sum/H_i0[i])
                    else
                        duH[i] = 0.0
                    end
                end
    
                for α in 1:R
                    if P[α] > 0.0
                        predation_sum = 0.0
                        for j in 1:S
                            predation_sum += a_matrix[j, α]*H[j]
                        end
                        predator_interactions = 0.0
                        for β in 1:R
                            predator_interactions += A[α, β]*P[β]
                        end
                        duP[α] = P[α]*(epsilon*predation_sum - m_alpha[α] + predator_interactions)
                    else
                        duP[α] = 0.0
                    end
                end
    
                du[1:S] = duH
                du[S+1:S+R] = duP
            end
    
            callbacks = []
            push!(callbacks, PositiveDomain())
    
            for i in 1:S
                condition(u, t, integrator) = u[i] - EXTINCTION_THRESHOLD
                affect!(integrator) = (integrator.u[i] = 0.0)
                push!(callbacks, ContinuousCallback(condition, affect!))
            end
    
            offset = S
            for α in 1:R
                idx = offset + α
                condition(u, t, integrator) = u[idx] - EXTINCTION_THRESHOLD
                affect!(integrator) = (integrator.u[idx] = 0.0)
                push!(callbacks, ContinuousCallback(condition, affect!))
            end
    
            cb = CallbackSet(callbacks...)
    
            sol = solve(ODEProblem(ecosystem_dynamics!, u0, tspan, params), Tsit5(); callback=cb, reltol=1e-6, abstol=1e-6)
    
            H_data = sol[1:S, :]
            P_data = sol[S+1:S+R, :]
    
            herb_survivors = count(H_data[:, end] .> EXTINCTION_THRESHOLD)
            pred_survivors = count(P_data[:, end] .> EXTINCTION_THRESHOLD)
    
            total_species = S + R
            surviving_species = herb_survivors + pred_survivors
    
            all_survived = (surviving_species == total_species)
            proportion_of_survived = surviving_species / total_species
    
            H_biomass = sum(H_data[:, end][H_data[:, end] .> EXTINCTION_THRESHOLD])
            P_biomass = sum(P_data[:, end][P_data[:, end] .> EXTINCTION_THRESHOLD])
    
            herb_pred_bodymass_ratio = H_biomass / (P_biomass != 0 ? P_biomass : Inf)
            total_biomass = H_biomass + P_biomass
    
            # Check if sum((g_i+G_i)*H_i) ≈ NPP
            H_final = H_data[:, end]
            system_true = isapprox(sum((g_i .+ G) .* H_final), NPP; atol=1e-3)
        end
    end
    
    # No printing, just output variables:
    all_survived = all_survived
    proportion_of_survived = proportion_of_survived
    herb_pred_bodymass_ratio = herb_pred_bodymass_ratio
    total_biomass = total_biomass
    system_true = system_true
    convergent = convergent

end
