##### METRICS FOR THREADING QUICKLY #####
# Function to compute Shannon index for a given abundance vector
function shannon_index(abundance_vector)
    total_abundance = sum(abundance_vector)
    if total_abundance == 0
        return 0.0  # Handle case where total abundance is zero
    end
    proportions = abundance_vector / total_abundance
    return -sum(p * log(p) for p in proportions if p > 0)
end

# Function to compute the average Shannon index over specified indices
function average_shannon_index(array_output, position; modified = false, caca = false)
    if !modified && !caca
        # Merge birmmals and herps
        combined_abundances = (deepcopy(array_output[end].herps) + deepcopy(array_output[end].birmmals)).*lambda_DA[position]
    elseif modified && !caca
        combined_abundances = array_output[end].state.*lambda_DA[position]
    end
    if caca
        combined_abundances = array_output.*lambda_DA[position]
    end
    
    # Initialize a list to store Shannon indices
    shannon_indices = Float64[]

    # Iterate over the specified indices
    for index in idx
        cell_abundance_vector = combined_abundances[index].a
        if !any(isnan, cell_abundance_vector)
            push!(shannon_indices, shannon_index(cell_abundance_vector))
        end
    end

    # Compute and return the average Shannon index
    return mean(shannon_indices)
end
############### Mean Trophic Level #################
####################################################
TrophInd = CSV.File(joinpath(dir, "DFs/TLs.csv")) |> DataFrame
TrophInd = TrophInd[1:256, 1:2]
TrophInd[findall(x -> x < 1.05, TrophInd[:, 2]), 2] .= 1.0
# TrophInd[:, 2] = TrophInd[:, 2].-1
# TrophInd[256, 2] = 1.0 # For some reason last line had floating point error
rename!(TrophInd, Symbol("Column1") => :Species, Symbol("TL") => :TL)
TrophInd[findall(x -> 1.98 < x < 2.05, TrophInd[:, 2]), 2] .= 2.0
order_indices = indexin(spain_names, TrophInd[:, :Species])
TrophInd = TrophInd[order_indices, :]
TrophInd_vector = TrophInd[:, :TL]

# Function to calculate mean trophic level
function calculate_mean_tl(array_output, position; modified = false, caca = false)
    if !modified && !caca
        # Merge birmmals and herps
        combined_abundances = (deepcopy(array_output[end].herps) + deepcopy(array_output[end].birmmals)).*lambda_DA[position]
    elseif modified && !caca
        combined_abundances = array_output[end].state.*lambda_DA[position]
    end
    if caca
        combined_abundances = array_output.*lambda_DA[position]
    end
    
    meanTL_matrix = DimArray(reshape([NaN32 for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))
    
    for cell in idx
        abundances = combined_abundances[cell].a
        presence = abundances .> body_mass_vector
        # Calculate mean trophic level for the present species
        if sum(presence) > 0 && !any(isnan, abundances)
            meanTL_matrix[cell] = mean(TrophInd_vector[presence])
        else
            meanTL_matrix[cell] = NaN
        end
    end
    # Exclude NaN values from the mean
    return mean(filter(!isnan, meanTL_matrix))
end

############## BIOMASS DISTRIBUTION ################
####################################################
function average_bbp(array_output, position; modified = false, caca = false)
    if !modified && !caca
        # Merge birmmals and herps
        combined_abundances = (deepcopy(array_output[end].herps) + deepcopy(array_output[end].birmmals)).*lambda_DA[position]
    elseif modified && !caca
        combined_abundances = array_output[end].state.*lambda_DA[position]
    end
    if caca
        combined_abundances = array_output.*lambda_DA[position]
    end
    bbp_vector = Float64[]
    
    for cell in idx
        abundances = combined_abundances[cell].a
        presence = abundances .> body_mass_vector
        
        if sum(presence) > 0 && !any(isnan, abundances)
            present_abundances = abundances[presence]
            present_trophic_levels = TrophInd_vector[presence]
            total_biomass = sum(present_abundances)
            weighted_trophic_sum = sum(present_abundances .* present_trophic_levels)
            bbp_vector = push!(bbp_vector, weighted_trophic_sum / total_biomass)
        else
            bbp_vector = push!(bbp_vector, NaN)
        end
    end
    return mean(filter(!isnan, bbp_vector))
end

######### RICHNESS SIMILARITY ################
function richness_similarity(array_output, position; modified = false, caca = false)
    if !modified && !caca
        # Merge birmmals and herps
        combined_abundances = (deepcopy(array_output[end].herps) + deepcopy(array_output[end].birmmals)).*lambda_DA[position]
    elseif modified && !caca
        combined_abundances = array_output[end].state.*lambda_DA[position]
    end
    if caca
        combined_abundances = array_output.*lambda_DA[position]
    end

    # Create a matrix to store simulated species richness
    simulated_richness = DimArray(reshape([0.0 for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))
    indices_to_remove = [
    CartesianIndex(34, 9), CartesianIndex(34, 10),
    CartesianIndex(34, 11), CartesianIndex(34, 12), 
    CartesianIndex(34, 13)
    ]

    # Use filter! to remove unwanted indices from idx
    idx_removed = filter!(x -> !(x in indices_to_remove), idx)

    # Calculate presence/absence and simulated richness
    for cell in idx_removed
        if !any(isnan, combined_abundances[cell].a)
            abundances = combined_abundances[cell].a
            presence = abundances .> body_mass_vector
            simulated_richness[cell] = sum(presence)
            # if simulated_richness[cell] != DA_richness[cell]
            #     print("cell is: ", cell, "\n")
            # end
            # println(simulated_richness[cell])
            # println(DA_richness[cell])
        elseif any(isnan, combined_abundances[cell].a)
            simulated_richness[cell] = 0.0
        end
    end

    # Calculate Mean Absolute Error (MAE) between real and simulated richness
    differences = [abs(DA_richness[cell] - simulated_richness[cell]) for cell in idx_removed]

    return mean(differences)
end

# # Iterate over the grid and create new SVector replacing 10.0 with 500.0
# triall = deserialize("Objects/DA_with_abundances_all10.jls")::DimArray{MyStructs256{Float64},2}
# for row in axes(triall, 1), col in axes(triall, 2)
#     old_vector = triall[row, col].a
#     new_vector = SVector{256, Float64}(replace(old_vector, 10.0 => 500.0))  # Create a new SVector
#     triall[row, col] = MyStructs256(new_vector)  # Assign the new MyStructs256 with updated vector
# end
# richness_similarity(Matrix(triall), modified = true, caca = true)

carnivores_vector = deepcopy(herb_carv_vector)
carnivores_vector[carnivores_vector .== 1.0] .= 0.0
carnivores_vector[carnivores_vector .== 0.00000001] .= 1.0

function alive_predators(array_output, position; modified = false, caca = false)
    if !modified && !caca
        # Merge birmmals and herps
        combined_abundances = (deepcopy(array_output[end].herps) + deepcopy(array_output[end].birmmals)).*lambda_DA[position]
    elseif modified && !caca
        combined_abundances = array_output[end].state.*lambda_DA[position]
    end
    if caca
        combined_abundances = array_output.*lambda_DA[position]
    end

    carn_alive = Float64[]
    # Calculate presence/absence and simulated richness
    for cell in idx
        if !any(isnan, combined_abundances[cell].a)
            abundances = combined_abundances[cell].a
            presence = abundances .> body_mass_vector
            
            simulated_predator_richness = sum(presence.*carnivores_vector)/106
            carn_alive = push!(carn_alive, simulated_predator_richness)
        end
    end

    return mean(carn_alive)
end
# alive_predators(triall, modified = true, caca = true)

function total_biomass(array_output, position; modified = false, caca = false)
    if !modified && !caca
        # Merge birmmals and herps
        combined_abundances = (deepcopy(array_output[end].herps) + deepcopy(array_output[end].birmmals)).*lambda_DA[position]
    elseif modified && !caca
        combined_abundances = array_output[end].state.*lambda_DA[position]
    end
    if caca
        combined_abundances = array_output.*lambda_DA[position]
    end
end
######### MEAN_N_OF_SPECIES ################
function mean_n_of_species(array_output, position; modified = false, caca = false)
    if !modified && !caca
        # Merge birmmals and herps
        combined_abundances = (deepcopy(array_output[end].herps) + deepcopy(array_output[end].birmmals)).*lambda_DA[position]
    elseif modified && !caca
        combined_abundances = array_output[end].state.*lambda_DA[position]
    end
    if caca
        combined_abundances = array_output.*lambda_DA[position]
    end
    vector = Float64[]
    # Calculate presence/absence and simulated richness
    for cell in idx
        if !any(isnan, combined_abundances[cell].a)
            abundances = combined_abundances[cell].a
            presence = abundances .> body_mass_vector
            simulated_richness = sum(presence)
            vector = push!(vector, simulated_richness)
        end
    end
    return mean(vector)
end
