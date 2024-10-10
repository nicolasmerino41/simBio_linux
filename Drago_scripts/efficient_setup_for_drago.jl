################################################################################################
####################### REAL SIMULATION ############################
####################################################################
####################################################################
####################################################################
######################## Size_selection kernel ########################
#######################################################################
function size_selection_kernel(predator_mass, prey_mass, sd, beta)
    intensity = exp(-((log(float(predator_mass)) / (beta * float(prey_mass)))^2.0) / (2.0 * float(sd)^2.0))
    return float(intensity)
end
beta = float(3)

gbif_sizes = CSV.read(joinpath(dir, "DFs/gbif_sizes.csv"), DataFrame)[:, 2:end]
###################### ZERO DIAGONAL ##################################
######################################################################
function zero_out_diagonal!(matrix)
    n = size(matrix, 1)
    for i in 1:n
        matrix[i, i] = 0
    end
    return matrix
end
#################### FILL DIAGONAL ###################################
######################################################################
# Assuming there exists a function with name `fill_diagonal!` to fill the diagonal of a matrix.
# If not, it needs to be defined as follows:
function fill_diagonal!(mat, val)
    for i in 1:min(size(mat)...)
        mat[i, i] = val
    end
end
#################### get_neighbors ###################################
######################################################################
# Helper function to get the neighbors
function get_neighbors(matrix, row, col)
    neighbors = []
    rows, cols = size(matrix)
    
    for r in max(row-1, 1):min(row+1, rows)
        for c in max(col-1, 1):min(col+1, cols)
            if (r != row || c != col) && !isnan(matrix[r, c])
                push!(neighbors, matrix[r, c])
            end
        end
    end
    return neighbors
end

web = CSV.read(joinpath(dir, "DFs/TetraEU_pairwise_interactions.csv"), DataFrame)

web = DataFrame(predator = web.sourceTaxonName, prey = web.targetTaxonName)

web.predator = string.(web.predator)
web.prey = string.(web.prey)
web = web[:, [:predator, :prey]]

unique_predators = unique(web.predator)
unique_preys = unique(web.prey)

x = vcat(unique_predators, unique_preys)
unique_species = unique(x)

# Read the CSV file
diets = CSV.File(joinpath(dir, "DFs/TetraEU_generic_diet.csv")) |> DataFrame
diets = hcat(diets.sourceTaxonName, diets.targetGenericItemName)

Amph = CSV.read(joinpath(dir, "DFs/DB_Amphibians_IP.txt"), delim='\t', DataFrame)
Bird = CSV.read(joinpath(dir, "DFs/DB_Birds_IP.txt"), delim='\t', DataFrame)
Mamm = CSV.read(joinpath(dir, "DFs/DB_Mammals_IP.txt"), delim='\t', DataFrame)
Rept = CSV.read(joinpath(dir, "DFs/DB_Reptiles_IP.txt"), delim='\t', DataFrame)

amphibian_names = names(Amph)[2:end]
reptile_names = names(Rept)[2:end]
mammal_names = names(Mamm)[2:end]
bird_names = names(Bird)[2:end]
herps_names = append!(deepcopy(amphibian_names), deepcopy( reptile_names))
birmmals_names = append!(deepcopy(mammal_names), deepcopy(bird_names))
spain_fauna = append!(deepcopy(herps_names), deepcopy(birmmals_names)) 

# Refactored code to merge columns of DataFrames for Spanish fauna, keeping UTMCODE just once
spanish_fauna = hcat(Amph, Rept[:, Not(:UTMCODE)], Mamm[:, Not(:UTMCODE)], Bird[:, Not(:UTMCODE)])

# Merge the `web` DataFrame with `spain_fauna` using inner join on 'predator' from `web` and 'species' from `spain_fauna`
merged_web = innerjoin(web, DataFrame(species=spain_fauna), on=(:predator => :species))
# Filter the merged DataFrame based on the prey column to include only species found in Spain
merged_web = merged_web[in.(merged_web.prey, Ref(spain_fauna)), :]
# Obtaining the species names that are at least predator/prey
unique_species_in_web = unique(vcat(merged_web.predator, merged_web.prey))
# println("There are ", length(unique_species_in_web), " unique species in the food web")

# Initializing an empty matrix with zeros for the Iberian species interaction
n = length(unique_species_in_web)
iberian_interact_matrix = zeros(Int, n, n)
iberian_interact_matrix = NamedArray(iberian_interact_matrix, (unique_species_in_web, unique_species_in_web))
# Ordering the matrix by Apmhibians, Reptiles, Mammals, Birds
spain_names = filter(name -> name in unique_species_in_web, names(hcat(Amph[:, Not(:UTMCODE)], Rept[:, Not(:UTMCODE)], Mamm[:, Not(:UTMCODE)], Bird[:, Not(:UTMCODE)])))
iberian_interact_matrix = iberian_interact_matrix[:, spain_names]
iberian_interact_matrix = iberian_interact_matrix[spain_names, :]

## Creating a mapping from species names to matrix indices
species_to_index = Dict(zip(spain_names, 1:n))
species_names = collect(keys(species_to_index))
#Filling the matrix with 1s where there are predator-prey interactions
for i in 1:nrow(merged_web)
    iberian_interact_matrix[merged_web.predator[i], merged_web.prey[i]] = 1
end
iberian_interact_matrix = float(iberian_interact_matrix)
# Count the amount of 1s in the iberian_interact_matrix
interaction_count = sum(iberian_interact_matrix .== 1)
# println("The number of 1s in the iberian_interact_matrix is $interaction_count")

#### STABLISHING TROPHIC LEVELS AND HERBIVORE NAMES
TrophInd = CSV.File(joinpath(dir, "DFs/TLs.csv")) |> DataFrame
TrophInd = TrophInd[1:256, 2:3]
TrophInd[:, 1] = round.(TrophInd[:, 1], digits = 2)
# TrophInd[:, 2] = TrophInd[:, 2].-1
# TrophInd[256, 2] = 1.0 # For some reason last line had floating point error
rename!(TrophInd, Symbol("species") => :Species, Symbol("TrophInd") => :TL)
TrophInd[:, :TL] = round.(TrophInd[:, :TL].-1, digits = 2)
order_indices = indexin(spain_names, TrophInd[:, :Species])
TrophInd = TrophInd[order_indices, :]
TrophInd_vector = TrophInd[:, :TL]

herbivore_names = TrophInd[TrophInd[:, :TL] .== 1, :Species]

# Turn iberian_interact_matrix into a DataFrame
# Convert the iberian_interact_matrix into a DataFrame with appropriate column and row names
self_regulation = 0.001
iberian_interact_df = DataFrame(iberian_interact_matrix, species_names)
function turn_adj_into_inter(adjacencyy, sigma, epsilon, self_regulation, beta)
    adjacency = deepcopy(adjacencyy)
    epsilon = float(epsilon)
    u = adjacency
    for i in names(adjacency, 1)
        for j in names(adjacency, 2)
            if adjacency[i, j] != 0.0 && i != j && adjacency[i, j] > 0.0 && adjacency[j, i] == 0.0
                predator_mass = Float64(gbif_sizes[gbif_sizes.species .== i, :bodyMass][1])
                prey_mass = Float64(gbif_sizes[gbif_sizes.species .== j, :bodyMass][1])
                sd = Float64(gbif_sizes[gbif_sizes.species .== i, :sigma][1])  # Use the sigma value for the predator
                # println(length(predator_mass))
                # Calculate interaction strength based on size-selection kernel
                kernel = size_selection_kernel(predator_mass, prey_mass, sd, beta)
                intensity = max(0.001*sd, kernel)
                # We need to create a Normal distribution first
                normal_dist = Normal(0, sigma*intensity)
                x = round(abs(rand(normal_dist)), digits = 20)
                u[i, j] = x / epsilon
                u[j, i] = -x
            elseif adjacency[i, j] != 0.0 && i != j && adjacency[i, j] > 0.0 && adjacency[j, i] > 0.0
                predator_mass = Float64(gbif_sizes[gbif_sizes.species .== i, :bodyMass][1])
                prey_mass = Float64(gbif_sizes[gbif_sizes.species .== j, :bodyMass][1])
                sd = Float64(gbif_sizes[gbif_sizes.species .== i, :sigma][1])  # Use the sigma value for the predator
                # println(length(predator_mass))
                # Calculate interaction strength based on size-selection kernel
                kernel = size_selection_kernel(predator_mass, prey_mass, sd, beta)
                intensity = max(0.001 * sigma, kernel)
                
                # Draw from a semi-Gaussian distribution
                normal_dist = Normal(0, sigma * intensity)
                x = round(abs(rand(normal_dist)), digits = 20)
                u[i, j] = x / 4.0
            elseif i == j
                u[i, j] = -self_regulation
            end
        end
    end
    adjacency = u
    return adjacency
end

##################### UTMtoRASTER ##################################
####################################################################
####################################################################
####################################################################
species_df = CSV.File(joinpath(dir, "DFs/Species_spain_df.csv")) |> DataFrame

variables = species_df[!, 2:5]
rename!(variables, [:ID, :Value, :sum, :UTMCODE])

species = species_df[!, unique_species_in_web]
species = species[:, spain_names]

species_df = hcat(variables, species, makeunique=true)
species_df_matrix = Matrix(species_df)

utmraster = Raster(joinpath(dir, "Rasters/updated_utmraster.tif"))
utmraster_DA = DimArray(utmraster)
utmraster_da = map(x -> isnothing(x) || isnan(x) ? false : true, utmraster_DA)

DA = deserialize(joinpath(dir, "Objects/DA.jls"))
DA_herps = deserialize(joinpath(dir, "Objects/DA_herps.jls"))
DA_birmmals = deserialize(joinpath(dir, "Objects/DA_birmmals.jls"))
# @load "Objects1_9/DA.jld2" DA
# @load "Objects1_9/DA_herps.jld2" DA_herps
# @load "Objects1_9/DA_birmmals.jld2" DA_birmmals

initial_abundance = 0.41
DA_with_abundances = deepcopy(DA)
for row in axes(DA, 1), col in axes(DA, 2)
    if DA[row, col] != MyStructs256(SVector{256, Float64}(fill(0.0, 256)))
        new_a = SVector{256, Float64}([DA[row, col].a[i] != 0.0 ? initial_abundance : DA[row, col].a[i] for i in 1:256])
        DA_with_abundances[row, col] = MyStructs256(new_a)
    end
end

DA_sum = deserialize(joinpath(dir, "Objects/DA_sum.jls"))
# Turnin DA_sum(float) into boolean
# DA_sum = DA_sums .== 1.0

# @load "Objects1_9/DA_sum.jld2" DA_sum
DA_sum_r = reverse(DA_sum, dims=1)
DA_sum_p = permutedims(DA_sum, (2, 1))

DA_with_abundances_r = reverse(DA_with_abundances, dims=1)
DA_with_abundances_p = permutedims(DA_with_abundances, (2, 1))
DA_with_abundances_p_masked = deepcopy(DA_with_abundances_p)

DA_richness = deserialize(joinpath(dir, "Objects/DA_richness.jls"))::DimArray{Float64,2}
# @load "Objects1_9/DA_richness.jld2" DA_richness
# @load "Objects1_9/DA_richness_birmmals.jld2" DA_richness_birmmals
# @load "Objects1_9/DA_richness_herps.jld2" DA_richness_herps
########################## IDX #####################################
idx = findall(x -> x == 1, DA_sum)
DA_with_presences = DimArray([fill(0.0, 256) for _ in 1:125, _ in 1:76], (Dim{:a}(1:125), Dim{:b}(1:76)))

for row in axes(DA_with_abundances, 1), col in axes(DA_with_abundances, 2)
    if DA_with_abundances[row, col].b != 0.0
        for i in findall(!iszero, DA_with_abundances[row, col].a)
            DA_with_presences[row, col][i] = 1.0
        end
    end
end

######################### NPP ####################################
npp_absolute = CSV.File(joinpath(dir, "DFs/npp_absolute_df.csv")) |> DataFrame
rename!(npp_absolute, [:ID, :UTMCODE, :npp, :X, :Y]) 
npp_absolute_in_kg = deepcopy(npp_absolute)
npp_absolute_in_kg.npp = npp_absolute.npp .* 1000
npp_absolute_in_kg = npp_absolute_in_kg[:, [2, 3]]

species_df = leftjoin(species_df, npp_absolute_in_kg, on = :UTMCODE, makeunique = true)
species_df_matrix = Matrix(species_df)

npp_DA = deserialize(joinpath(dir, "Objects/npp_DA.jls"))
# @load "Objects1_9/npp_DA.jld2" npp_DA
# npp_DA = npp_DA./10000
################### EFFICIENT MATRIX FRAMEWORK #####################
####################################################################
####################################################################
####################################################################
# Load a DataFrame from a serialized file ('.jls' format).
# iberian_interact_df = deserialize(joinpath(dir, "Objects/iberian_interact_df.jls"))
iberian_interact_df = CSV.File("Objects/iberian_interact_df.csv") |> DataFrame
# Convert the DataFrame to a matrix for easier manipulation.
iberian_interact_matrix = iberian_interact_df |> Matrix
# Convert the modified matrix back to a DataFrame, preserving the original column names.
iberian_interact_df = DataFrame(iberian_interact_matrix, names(iberian_interact_df))
# Create a NamedArray from the matrix, using the DataFrame's column names for both dimensions.
iberian_interact_NA = NamedArray(
    iberian_interact_matrix, 
    (names(iberian_interact_df), names(iberian_interact_df)),
    ("Species", "Species")
)
iberian_interact_NA = iberian_interact_NA[spain_names, spain_names]

##################### NEW NICHES ###########################
######## bio rasters  ##############
bio5_DA = deserialize(joinpath(dir, "Objects/bio5.jls"))
bio6_DA = deserialize(joinpath(dir, "Objects/bio6.jls"))
bio12_DA = deserialize(joinpath(dir, "Objects/bio12.jls"))
futurebio5_DA = bio5_DA .+ 1.0
futurebio6_DA = bio6_DA .+ 1.0
futurebio12_DA = bio12_DA .+ rand(Normal(0, 100), 125, 76)

######## body_mass_vector  ##############
# Initialize an empty vector
body_mass_vector = Float64[]
# Loop through each species in spain_names and push the bodyMass into the vector
for i in spain_names
    # Find the bodyMass for the species and push it to the vector
    body_mass = gbif_sizes[gbif_sizes.species .== i, :bodyMass]
    if !isempty(body_mass)
        push!(body_mass_vector, body_mass[1])
    end
end
body_mass_vector_herps = body_mass_vector[1:49]
body_mass_vector_birds = body_mass_vector[50:256]

######## niches_df  ##############
species_niches = CSV.File(joinpath(dir, "DFs/iberian_species_niches_withbinned_TH.csv"), decimal = ',') |> DataFrame
order_indices = indexin(spain_names, species_niches[:, :Species])
species_niches = species_niches[order_indices, :]

lax_species_niches = CSV.File(joinpath(dir, "DFs/iberian_species_niches_withLaxNiche.csv"), decimal = ',') |> DataFrame
order_indices = indexin(spain_names, lax_species_niches[:, :Species])
lax_species_niches = lax_species_niches[order_indices, :]

strict_species_niches = CSV.File(joinpath(dir, "DFs/iberian_species_niches_withVeryStrictNiche.csv"), decimal = ',') |> DataFrame
order_indices = indexin(spain_names, strict_species_niches[:, :Species])
strict_species_niches = strict_species_niches[order_indices, :]

# herbivore_names = CSV.File(joinpath(dir, "DFs/herbivore_names.csv")) |> DataFrame
# herbivore_names = convert(Vector{String}, herbivore_names[:, 2])
# binary_vector = [name in herbivore_names ? 1 : 0 for name in names(iberian_interact_df)]
# opposite_binary_vector = [name in herbivore_names ? 0 : 1 for name in names(iberian_interact_df)]

###################### GENERATE Ki(z) DimArray ################
###############################################################
# Initialize the empty DA arrays for each suitability method (same shape as k_DA)
DA_multiplicative = DimArray(reshape([MyStructs256(SVector{256, Float64}(fill(0.0, 256))) for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))
DA_additive = DimArray(reshape([MyStructs256(SVector{256, Float64}(fill(0.0, 256))) for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))
DA_geometric = DimArray(reshape([MyStructs256(SVector{256, Float64}(fill(0.0, 256))) for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))
DA_min = DimArray(reshape([MyStructs256(SVector{256, Float64}(fill(0.0, 256))) for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))
DA_harmonic = DimArray(reshape([MyStructs256(SVector{256, Float64}(fill(0.0, 256))) for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))

# Define the herbivore carnivore vector
herb_carv_svector = SVector{256, Float64}([name in herbivore_names ? 1.0 : 0.00000001 for name in spain_names])
herb_carv_vector = [name in herbivore_names ? 1.0 : 0.00000001 for name in spain_names]
@everywhere herb_carv_svector = deepcopy(herb_carv_svector)
@everywhere herb_carv_vector = deepcopy(herb_carv_vector)

# Loop through the axes of the DA arrays
for row in axes(DA_multiplicative, 1), col in axes(DA_multiplicative, 2)
    if isone(DA_sum[row, col])
        # Calculate the scaled deviations for bio5, bio6, and bio12
        S_bio5 = 1 ./ (1 .+ abs.(bio5_DA[row, col] .- lax_species_niches.mean_bio5) ./ lax_species_niches.sd_bio5)
        S_bio6 = 1 ./ (1 .+ abs.(bio6_DA[row, col] .- lax_species_niches.mean_bio6) ./ lax_species_niches.sd_bio6)
        S_bio12 = 1 ./ (1 .+ abs.(bio12_DA[row, col] .- lax_species_niches.mean_bio12) ./ lax_species_niches.sd_bio12)

        # 1. Multiplicative Approach (Original)
        multiplicative_suitability = S_bio5 .* S_bio6 .* S_bio12
        DA_multiplicative[row, col] = MyStructs256(SVector{256, Float64}(multiplicative_suitability .* herb_carv_svector))

        # 2. Additive Approach
        additive_suitability = (S_bio5 .+ S_bio6 .+ S_bio12) / 3
        DA_additive[row, col] = MyStructs256(SVector{256, Float64}(additive_suitability .* herb_carv_svector))

        # 3. Geometric Mean Approach
        geometric_suitability = (S_bio5 .* S_bio6 .* S_bio12).^(1/3)
        DA_geometric[row, col] = MyStructs256(SVector{256, Float64}(geometric_suitability .* herb_carv_svector))

        # 4. Minimum Suitability Approach
        min_suitability = min(S_bio5, S_bio6, S_bio12)
        DA_min[row, col] = MyStructs256(SVector{256, Float64}(min_suitability .* herb_carv_svector))

        # 5. Harmonic Mean Approach
        harmonic_suitability = 3 ./ (1 ./ S_bio5 .+ 1 ./ S_bio6 .+ 1 ./ S_bio12)
        DA_harmonic[row, col] = MyStructs256(SVector{256, Float64}(harmonic_suitability .* herb_carv_svector))
    end
end

k_DA = (DA_multiplicative = DA_multiplicative, DA_additive = DA_additive, DA_geometric = DA_geometric, DA_min = DA_min, DA_harmonic = DA_harmonic)

###################### GENERATE lambda_DA DimArray ################
###############################################################
# Function to calculate the lambda scalar for a given k_hat and NPP
function calculate_lambda_scalar(k_hat, NPP)
    lambda = NPP / sum(k_hat)
    return lambda
end

# Initialize a NamedTuple to store lambda_DA for each k_DA
lambda_DA = NamedTuple((
    multiplicative = DimArray(reshape([NaN64 for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76))),
    additive = DimArray(reshape([NaN64 for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76))),
    geometric = DimArray(reshape([NaN64 for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76))),
    minimum = DimArray(reshape([NaN64 for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76))),
    harmonic = DimArray(reshape([NaN64 for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))
))
# lambda_DA[1]

# Loop through the axes to calculate lambda_DA for each suitability method
for row in axes(lambda_DA.multiplicative, 1), col in axes(lambda_DA.multiplicative, 2)
    if isone(DA_sum[row, col])
        # Calculate lambda for multiplicative suitability
        lambda_DA.multiplicative[row, col] = calculate_lambda_scalar(k_DA.DA_multiplicative[row, col].a, npp_DA[row, col])
        # Calculate lambda for additive suitability
        lambda_DA.additive[row, col] = calculate_lambda_scalar(k_DA.DA_additive[row, col].a, npp_DA[row, col])
        # Calculate lambda for geometric suitability
        lambda_DA.geometric[row, col] = calculate_lambda_scalar(k_DA.DA_geometric[row, col].a, npp_DA[row, col])
        # Calculate lambda for minimum suitability
        lambda_DA.minimum[row, col] = calculate_lambda_scalar(k_DA.DA_min[row, col].a, npp_DA[row, col])
        # Calculate lambda for harmonic suitability
        lambda_DA.harmonic[row, col] = calculate_lambda_scalar(k_DA.DA_harmonic[row, col].a, npp_DA[row, col])
    end
end

dimensions = (125, 76)
idx_tupled = [(i, j) for i in 1:dimensions[1], j in 1:dimensions[2]]
