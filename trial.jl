using Distributions
using Random

# OpenMOLE-provided variables (assumed to be defined)
# total_biomass

# Set random seed for reproducibility (optional)
# Random.seed!(seed)

result = total_biomass * rand()
# Output the result (OpenMOLE will capture the variable 'result')
println("result = ", result)
