using NamedGraphs.NamedGraphGenerators: named_grid
using Statistics: mean

# Load the Plots library for plotting results
using Plots: Plots, plot

include("ising_tensornetwork.jl")
include("belief_propagation.jl")

"""
    main(; kwargs...)

Creates a grid graph of size `Lx` by `Ly`, constructs the Ising tensor network on it, and
computes the Bethe-Peierls free energy density using belief propagation.

# Keywords
- `Lx::Int = 3`: The number of columns in the grid.
- `Ly::Int = 3`: The number of rows in the grid.
- `beta::Number = 0.2`: The inverse temperature parameter.
- `periodic::Bool = false`: Whether to use periodic boundary conditions.
- `outputlevel::Int = 1`: Controls how much information will be printed by the script.

# Returns
A named tuple containing:
- `bp_phi_tn::Number`: The Bethe-Peierls free energy density computed via belief propagation.
- `exact_phi_onsager::Number`: The exact free energy density from Onsager's solution in the thermodynamic limit.
- `niters::Int`: The number of iterations taken for convergence in belief propagation.
"""
function main(; Lx::Int = 3, Ly::Int = 3, beta::Number = 0.2, periodic = false, outputlevel::Int = 1)
    g = named_grid((Lx, Ly); periodic)

    tensornetwork = ising_tensornetwork(g, beta)
    messages, niters = belief_propagation(tensornetwork, g; niters = 1000, outputlevel)

    bp_phi_tn = bp_phi(tensornetwork, messages, g)
    exact_phi_onsager = ising_phi(beta)
    return (; bp_phi_tn, exact_phi_onsager, niters)
end
