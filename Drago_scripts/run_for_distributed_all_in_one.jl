# Packages
using Distributed
using Base.Iterators:product


addprocs(4,exeflags=["--project","--threads=2"])

@everywhere begin
    using CSV, DataFrames
    using Distributions, NamedArrays, StaticArrays
    using DynamicGrids, Dispersal
    using DimensionalData, Rasters, Serialization, ArchGDAL
    using OrderedCollections, StatsBase

    dir = pwd()
    num_species = 256
end

@everywhere begin
    ######################## COMPLEX RULES #############################
####################################################################
####################################################################
####################################################################
######################## DEFINING BASIC MYSTRUCTS256 METHODS ####################################
#################################################################################################
struct MyStructs256{T <: AbstractFloat} <: FieldVector{2, T}
    a::SVector{num_species, T}
    b::T

    # Custom constructor for automatic sum calculation
    function MyStructs256(a::SVector{num_species, T}) where T <: AbstractFloat
        new{T}(a, sum(a))
    end

    # Explicit constructor allowing manual setting of both `a` and `b`
    MyStructs256(a::SVector{num_species, T}, b::T) where T <: AbstractFloat = new{T}(a, b)
end

# Define zero and oneunit for MyStructs256
Base.zero(::Type{MyStructs256{T}}) where {T <: AbstractFloat} = MyStructs256(SVector{num_species, T}(fill(zero(T), num_species)), zero(T))
Base.zero(x::MyStructs256{T}) where {T <: AbstractFloat} = MyStructs256(SVector{num_species, T}(fill(zero(T), num_species)), zero(T))
Base.oneunit(::Type{MyStructs256{T}}) where {T <: AbstractFloat} = MyStructsnum_species(fill(oneunit(T), num_species), oneunit(T))
Base.iszero(x::MyStructs256{T}) where {T <: AbstractFloat} = all(iszero, x.a) && iszero(x.b)
# Comparison based on 'b' field
Base.isless(x::MyStructs256, y::MyStructs256) = isless(x.b, y.b)
Base.isless(x::MyStructs256, y::AbstractFloat) = isless(x.b, y)

# Element-wise arithmetic operations ensuring 'b' is recalculated correctly
Base.:+(x::MyStructs256, y::MyStructs256) = MyStructs256(x.a .+ y.a, sum(x.a .+ y.a))
Base.:-(x::MyStructs256, y::MyStructs256) = MyStructs256(x.a .- y.a, sum(x.a .- y.a))
Base.:*(x::MyStructs256, scalar::Real) = MyStructs256(x.a .* scalar, sum(x.a .* scalar))
Base.:/(x::MyStructs256, scalar::Real) = MyStructs256(x.a ./ scalar, sum(x.a ./ scalar))
Base.:-(x::MyStructs256, scalar::Real) = MyStructs256(x.a .- scalar, x.b - scalar * num_species)
Base.:+(x::MyStructs256, scalar::Real) = MyStructs256(x.a .+ scalar, x.b + scalar * num_species)
Base.:*(x::MyStructs256, y::AbstractVector) = MyStructs256(x.a .* SVector{num_species, Float64}(y), sum(x.a .* SVector{num_species, Float64}(y)))
Base.:/(x::MyStructs256, y::AbstractVector) = MyStructs256(x.a ./ SVector{num_species, Float64}(y), sum(x.a ./ SVector{num_species, Float64}(y)))
Base.:*(y::AbstractVector, x::MyStructs256) = MyStructs256(SVector{num_species, Float64}(y) .* x.a, sum(SVector{num_species, Float64}(y) .* x.a))
Base.:/(y::AbstractVector, x::MyStructs256) = MyStructs256(SVector{num_species, Float64}(y) ./ x.a, sum(SVector{num_species, Float64}(y) ./ x.a))


# Define what a NaN is for MyStructs256
Base.isnan(x::MyStructs256) = isnan(x.b) || any(isnan, x.a)

# Adding a method in the sum function for MyStructs256
function Base.sum(structs::MyStructs256...)
    # Sum the 'a' vectors
    summed_a = sum([s.a for s in structs])

    # Sum the 'b' values
    summed_b = sum([s.b for s in structs])

    # Create a new MyStructs256 instance with the summed results
    return MyStructs256(summed_a, summed_b)
end

# Adding a method to maximum
# Define maximum for MyStructs256
function Base.maximum(a::MyStructs256, b::MyStructs256)
    return MyStructs256(max.(a.a, b.a))
end

# Define maximum for MyStructs256 with a scalar
function Base.maximum(a::MyStructs256, b::AbstractFloat)
    return MyStructs256(max.(a.a, b))
end

# Define maximum for a scalar with MyStructs256
function Base.maximum(a::AbstractFloat, b::MyStructs256)
    return MyStructs256(max.(a, b.a))
end

# Define maximum for MyStructs256
function Base.maximum(a::MyStructs256)
    return maximum(a.a)
end

# Define maximum for a matrix of MyStructs256
function Base.maximum(a::Matrix{MyStructs256{AbstractFloat}})
    # Extract all `b` values from each MyStructs256 element in the matrix and find the maximum
    return maximum(map(x -> x.b, a))
end

################################## MYHERPS METHODS ###########################################
#################################################################################################
#################################################################################################
#################################################################################################
struct MyHerps{T <: AbstractFloat} <: FieldVector{2, T}
    a::SVector{49, T}
    b::T

    # Custom constructor for automatic sum calculation
    function MyHerps(a::SVector{49, T}) where T <: AbstractFloat
        new{T}(a, sum(a))
    end

    # Explicit constructor allowing manual setting of both `a` and `b`
    MyHerps(a::SVector{49, T}, b::T) where T <: AbstractFloat = new{T}(a, b)
end

# Define zero and oneunit for MyHerps
Base.zero(::Type{MyHerps{T}}) where {T <: AbstractFloat} = MyHerps(SVector{49, T}(fill(zero(T), 49)), zero(T))
Base.zero(x::MyHerps{T}) where {T <: AbstractFloat} = MyHerps(SVector{49, T}(fill(zero(T), 49)), zero(T))
Base.oneunit(::Type{MyHerps{T}}) where {T <: AbstractFloat} = MyHerps(fill(oneunit(T), 49), oneunit(T))

# Comparison based on 'b' field
Base.isless(x::MyHerps, y::MyHerps) = isless(x.b, y.b)
Base.isless(x::MyHerps, y::AbstractFloat) = isless(x.b, y)

# Element-wise arithmetic operations ensuring 'b' is recalculated correctly
Base.:+(x::MyHerps, y::MyHerps) = MyHerps(x.a .+ y.a, sum(x.a .+ y.a))
Base.:-(x::MyHerps, y::MyHerps) = MyHerps(x.a .- y.a, sum(x.a .- y.a))
Base.:*(x::MyHerps, scalar::Real) = MyHerps(x.a .* scalar, sum(x.a .* scalar))
Base.:/(x::MyHerps, scalar::Real) = MyHerps(x.a ./ scalar, sum(x.a ./ scalar))
Base.:-(x::MyHerps, scalar::Real) = MyHerps(x.a .- scalar, x.b - scalar * 49)
Base.:+(x::MyHerps, scalar::Real) = MyHerps(x.a .+ scalar, x.b + scalar * 49)

# Define what a NaN is for MyHerps
Base.isnan(x::MyHerps) = isnan(x.b) || any(isnan, x.a)

# Adding a method in the sum function for MyHerps
function Base.sum(structs::MyHerps...)
    # Sum the 'a' vectors
    summed_a = sum([s.a for s in structs])

    # Sum the 'b' values
    summed_b = sum([s.b for s in structs])

    # Create a new MyHerps instance with the summed results
    return MyHerps(summed_a, summed_b)
end

# Adding a method to maximum
# Define maximum for MyHerps
function Base.maximum(a::MyHerps, b::MyHerps)
    return MyHerps(max.(a.a, b.a))
end

# Define maximum for MyHerps with a scalar
function Base.maximum(a::MyHerps, b::AbstractFloat)
    return MyHerps(max.(a.a, b))
end

# Define maximum for a scalar with MyHerps
function Base.maximum(a::AbstractFloat, b::MyHerps)
    return MyHerps(max.(a, b.a))
end

# Define maximum for MyHerps
function Base.maximum(a::MyHerps)
    return maximum(a.a)
end

# Define maximum for a matrix of MyHerps
function Base.maximum(a::Matrix{MyHerps{AbstractFloat}})
    # Extract all `b` values from each MyHerps element in the matrix and find the maximum
    return maximum(map(x -> x.b, a))
end

################################## MYBIRMMALS METHODS ###########################################
#################################################################################################
#################################################################################################
#################################################################################################
struct MyBirmmals{T <: AbstractFloat} <: FieldVector{2, T}
    a::SVector{207, T}
    b::T

    # Custom constructor for automatic sum calculation
    function MyBirmmals(a::SVector{207, T}) where T <: AbstractFloat
        new{T}(a, sum(a))
    end

    # Explicit constructor allowing manual setting of both `a` and `b`
    MyBirmmals(a::SVector{207, T}, b::T) where T <: AbstractFloat = new{T}(a, b)
end

# Define zero and oneunit for MyBirmmals
Base.zero(::Type{MyBirmmals{T}}) where {T <: AbstractFloat} = MyBirmmals(SVector{207, T}(fill(zero(T), 207)), zero(T))
Base.zero(x::MyBirmmals{T}) where {T <: AbstractFloat} = MyBirmmals(SVector{207, T}(fill(zero(T), 207)), zero(T))
Base.oneunit(::Type{MyBirmmals{T}}) where {T <: AbstractFloat} = MyBirmmals(fill(oneunit(T), 207), oneunit(T))

# Comparison based on 'b' field
Base.isless(x::MyBirmmals, y::MyBirmmals) = isless(x.b, y.b)
Base.isless(x::MyBirmmals, y::AbstractFloat) = isless(x.b, y)

# Element-wise arithmetic operations ensuring 'b' is recalculated correctly
Base.:+(x::MyBirmmals, y::MyBirmmals) = MyBirmmals(x.a .+ y.a, sum(x.a .+ y.a))
Base.:-(x::MyBirmmals, y::MyBirmmals) = MyBirmmals(x.a .- y.a, sum(x.a .- y.a))
Base.:*(x::MyBirmmals, scalar::Real) = MyBirmmals(x.a .* scalar, sum(x.a .* scalar))
Base.:/(x::MyBirmmals, scalar::Real) = MyBirmmals(x.a ./ scalar, sum(x.a ./ scalar))
Base.:-(x::MyBirmmals, scalar::Real) = MyBirmmals(x.a .- scalar, x.b - scalar * 207)
Base.:+(x::MyBirmmals, scalar::Real) = MyBirmmals(x.a .+ scalar, x.b + scalar * 207)

# Define what a NaN is for MyBirmmals
Base.isnan(x::MyBirmmals) = isnan(x.b) || any(isnan, x.a)

# Adding a method in the sum function for MyBirmmals
function Base.sum(structs::MyBirmmals...)
    # Sum the 'a' vectors
    summed_a = sum([s.a for s in structs])

    # Sum the 'b' values
    summed_b = sum([s.b for s in structs])

    # Create a new MyBirmmals instance with the summed results
    return MyBirmmals(summed_a, summed_b)
end

# Adding a method to maximum
# Define maximum for MyBirmmals
function Base.maximum(a::MyBirmmals, b::MyBirmmals)
    return MyBirmmals(max.(a.a, b.a))
end

# Define maximum for MyBirmmals with a scalar
function Base.maximum(a::MyBirmmals, b::AbstractFloat)
    return MyBirmmals(max.(a.a, b))
end

# Define maximum for a scalar with MyBirmmals
function Base.maximum(a::AbstractFloat, b::MyBirmmals)
    return MyBirmmals(max.(a, b.a))
end

# Define maximum for MyBirmmals
function Base.maximum(a::MyBirmmals)
    return maximum(a.a)
end

# Define maximum for a matrix of MyBirmmals
function Base.maximum(a::Matrix{MyBirmmals{AbstractFloat}})
    # Extract all `b` values from each MyBirmmals element in the matrix and find the maximum
    return maximum(map(x -> x.b, a))
end

# Define zeros for all three types
function Base.zeros(dims::NTuple{2, Int}, type = nothing)
    if type == MyBirmmals{AbstractFloat}
        return [MyBirmmals(fill(0.0, 207)) for _ in 1:dims[1], _ in 1:dims[2]]
    elseif type == MyStructs256{AbstractFloat}
        return [MyStructs256(fill(0.0, 256)) for _ in 1:dims[1], _ in 1:dims[2]]
    elseif type == MyHerps{AbstractFloat}
        return [MyHerps(fill(0.0, 49)) for _ in 1:dims[1], _ in 1:dims[2]]
    else
        return [0.0 for _ in 1:dims[1], _ in 1:dims[2]]
    end
end

################# MYBIRMMALS KERNEL METHODS ################
###############################################################
###############################################################
# Define kernel product for MyStructs256
function Dispersal.kernelproduct(hood::Window{1, 2, 9, MyBirmmals{AbstractFloat}}, kernel::SVector{9, AbstractFloat})
    
    result_a = SVector{207, AbstractFloat}(fill(0.0f0, 207))
    
    for (i, k) in enumerate(kernel)
        result_a += hood[i].a .* k
    end
    return MyBirmmals(result_a)
end

function Dispersal.kernelproduct(hood::Window{2, 2, 25, MyBirmmals{AbstractFloat}}, kernel::SVector{25, AbstractFloat})
    
    result_a = SVector{207, AbstractFloat}(fill(0.0f0, 207))
    
    for (i, k) in enumerate(kernel)
        result_a += hood[i].a .* k
    end

    return MyBirmmals(result_a)
end

import Base.Broadcast: broadcastable

broadcastable(x::MyBirmmals) = Ref(x)
broadcastable(x::MyHerps) = Ref(x)

Base.:+(x::MyHerps, y::MyBirmmals) = MyStructs256(SVector{256, typeof(x.b)}(vcat(x.a, y.a)), x.b + y.b)
Base.:+(x::MyBirmmals, y::MyHerps) = MyStructs256(SVector{256, typeof(x.b)}(vcat(y.a, x.a)), y.b + x.b)
end

@everywhere begin
     ################# MYSTRUCTS256 KERNEL METHODS ################
###############################################################
###############################################################
struct CustomKernel <: KernelFormulation
    a::AbstractFloat
end

abstract type AbstractKernelNeighborhood end

struct CustomDispersalKernel{N<:DynamicGrids.Neighborhood, F<:KernelFormulation} <: AbstractKernelNeighborhood
    neighborhood::N
    formulation::F
end

function CustomDispersalKernel(; 
    neighborhood::DynamicGrids.Neighborhood=Moore(1), 
    formulation::KernelFormulation=CustomKernel(1.0)
)
    CustomDispersalKernel{typeof(neighborhood), typeof(formulation)}(neighborhood, formulation)
end

function (kernel::CustomKernel)(distance)
    return exp(-(distance^2) / (2*(kernel.a^2)))
end
end

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

hf = Rasters.Raster(joinpath(dir, "Rasters/wildareas-v3-2009-human-footprint.tif"))

hf_reprojected = resample(hf; to=utmraster)

hf_spain = map(x -> Float32(x), hf_reprojected)
hf_spain = map(x -> x == 128.0 ? 0.0 : x, hf_spain)

inverted_hf = 1 ./ deepcopy(hf_spain)
inverted_hf = map(x -> x == Inf ? 0.0 : x, inverted_hf)

######## FOR DA'S ##########
inverted_hf_DA = DimArray(Matrix(inverted_hf), (Dim{:a}(1:125), Dim{:b}(1:76)))
# Ensure values below 0.1 (but greater than 0.0) in inverted_hf_DA are set to 0.1
inverted_hf_DA .= ifelse.((inverted_hf_DA .> 0.0) .& (inverted_hf_DA .< 0.1), 0.1, inverted_hf_DA)

function adjust_inverted_hf_DA(lambda_DA, inverted_hf_DA)
    inverted_hf_DA1 = deepcopy(inverted_hf_DA)
    rows, cols = size(lambda_DA.multiplicative) # Any lambda is ok, it's just for size matching
    adjusted_hf_DA = deepcopy(inverted_hf_DA1)  # Start with a copy of the original matrix

    # First pass: Calculate averages
    for row in 1:rows
        for col in 1:cols
            if !isnan(lambda_DA.multiplicative[row, col])
                neighbors = get_neighbors(inverted_hf_DA1, row, col)
                non_zero_neighbors = filter(x -> x != 0.0, neighbors)

                if !isempty(non_zero_neighbors)
                    adjusted_hf_DA[row, col] = mean(non_zero_neighbors)
                else
                    adjusted_hf_DA[row, col] = 0.1  # or some other default value if all neighbors are zero
                end
            end
        end
    end

    # Second pass: Set NA cells in lambda_DA to 0.0 in adjusted_hf_DA
    for row in 1:rows
        for col in 1:cols
            if isnan(lambda_DA.multiplicative[row, col])
                adjusted_hf_DA[row, col] = 0.0
            end
        end
    end

    return adjusted_hf_DA
end

# Adjust the inverted_hf_DA raster based on lambda_DA
adjusted_inverted_hf_DA = adjust_inverted_hf_DA(lambda_DA, inverted_hf_DA)

k_DA_hf_multiplicative = k_DA.DA_multiplicative .* adjusted_inverted_hf_DA
k_DA_hf_additive = k_DA.DA_additive .* adjusted_inverted_hf_DA
k_DA_hf_geometric = k_DA.DA_geometric .* adjusted_inverted_hf_DA
k_DA_hf_minimum = k_DA.DA_min .* adjusted_inverted_hf_DA
k_DA_hf_harmonic = k_DA.DA_harmonic .* adjusted_inverted_hf_DA

@everywhere begin
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
end

@everywhere function turn_comp_into_inter(adjacencyy, sigma, assymetry, self_regulation = 0.0)
    adjacency = deepcopy(adjacencyy)
    u = adjacency  # Initialize the interaction matrix
    
    # Loop through the matrix
    for i in names(adjacency, 1)
        for j in names(adjacency, 2)
            # Ensure we only process species pairs where there is competition (non-zero entries)
            if i != j && adjacency[i, j] != 0.0 && adjacency[j, i] != 0.0
                # Skip if the reverse interaction is already calculated
                if u[j, i] != 0.0 && u[j, i] != 1.0
                    continue
                end

                # Draw interaction strength from normal distribution
                normal_dist = Normal(0, sigma)
                x = round(abs(rand(normal_dist)), digits = 20)  # Interaction for i, j

                # For j, i: adjust based on assymetry (correlation factor)
                y = assymetry * x + (1 - assymetry) * round(abs(rand(normal_dist)), digits = 10)

                # Assign interactions with asymmetry
                u[i, j] = - x  # Species i competes with j
                u[j, i] = - y  # Species j competes with i (asymmetrically based on assymetry)

            elseif i == j
                u[i, j] = -self_regulation  # Self-regulation on diagonal
            end
        end
    end
    return u
end

########## DISTRIBUTING ##########
@everywhere DA = DA
@everywhere k_DA_hf_multiplicative = k_DA_hf_multiplicative
@everywhere npp_DA = npp_DA
@everywhere DA_richness = DA_richness

initial_abundance = 0.41
DA_with_abundances = deepcopy(DA)
for row in axes(DA, 1), col in axes(DA, 2)
    if !iszero(DA[row, col]) 
        new_a = SVector{256, Float64}([DA[row, col].a[i] != 0.0 ? initial_abundance : DA[row, col].a[i] for i in 1:256])
        DA_with_abundances[row, col] = MyStructs256(new_a)
    end
end

@everywhere const DA_with_abundances_const = deepcopy(DA_with_abundances)
@everywhere const iberian_interact_NA_const = deepcopy(iberian_interact_NA)
# Function to save parameters, grid type (k_DA name), and metrics to CSV and append plots to the final PDFs
@everywhere function run_simulation(sigma, epsilon, alfa, sigma_comp, assymetry)

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
        pepe_state; tspan = 1:10,
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
    serialize("resultados_distribuidos/outputs/$sigma-$epsilon-$alfa-$sigma_comp-$assymetry.jls", p[end].state)
    # Append or create the CSV file
    csv_filename = "resultados_distribuidos/DirectSamplingResults.csv"
    if isfile(csv_filename)
        CSV.write(csv_filename, results_row, append = true)
    else
        CSV.write(csv_filename, results_row)
    end
end

# Simulation parameters
sigmas = [0.0001, 0.001, 0.005, 0.01, 0.05, 0.07, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.0, 1.5, 2.0, 3.0]
epsilons = [0.01, 0.1, 0.5, 1.0, 2.0, 3.0, 5.0]
alfa_values = [0.001, 0.01, 0.05, 0.1, 0.3, 0.6, 0.9, 1.1]
sigma_comp_values = [0.001, 0.1, 0.5, 1.0, 1.5, 2.0]
assymetry_values = [0.0, 0.33, 0.66, 1.0]

k_DA_list = [k_DA.DA_multiplicative, k_DA.DA_additive, k_DA.DA_min, k_DA.DA_harmonic, k_DA.DA_geometric]
k_DA_names = ["multiplicative", "additive", "min", "harmonic", "geometric"]
positions = [1, 2, 3, 4, 5]

iteracciones = collect(product(sigmas, epsilons, alfa_values, sigma_comp_values, assymetry_values))
result = pmap(args->run_simulation(args...),iteracciones)