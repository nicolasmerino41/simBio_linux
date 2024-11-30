 ################# MYSTRUCTS256 KERNEL METHODS ################
###############################################################
###############################################################
struct CustomKernel <: KernelFormulation
    α::AbstractFloat
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
    return exp(-(distance^2) / (2*(kernel.α^2)))
end

