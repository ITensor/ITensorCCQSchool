using NamedGraphs.NamedGraphGenerators: named_grid
using Statistics: mean

# Load the Plots library for plotting results
using Plots: Plots, plot

include("ising_tensornetwork.jl")
include("belief_propagation.jl")

"""
    main(; kwargs...)  

Creates a periodic 5x5 grid graph, constructs the Ising tensor network on it, and computes
both the Bethe-Peierls free energy density using belief propagation and a loop-corrected
free energy density.

# Keywords
- `beta::Number = 0.2`: The inverse temperature parameter.
- `outputlevel::Int = 1`: Controls how much information will be printed by the script.

# Returns
A named tuple containing:
- `bp_phi_tn::Number`: The Bethe-Peierls free energy density computed via belief propagation.
- `bp_corrected_phi_tn::Number`: The loop-corrected free energy density.
- `exact_phi_onsager::Number`: The exact free energy density from Onsager's solution in the thermodynamic limit.
- `niters::Int`: The number of iterations taken for convergence in belief propagation.
"""
function main(; beta::Number = 0.2, outputlevel::Int = 1)
    graph = named_grid((5, 5); periodic = true)

    tensornetwork = ising_tensornetwork(graph, beta)
    messages, niters = belief_propagation(tensornetwork, graph; niters = 1000, outputlevel)

    bp_phi_tn = bp_phi(tensornetwork, messages, graph)
    smallest_loop_size = 4
    bp_corrected_phi_tn = bp_corrected_phi(tensornetwork, messages, graph; smallest_loop_size)
    exact_phi_onsager = ising_phi(beta)

    return (; bp_phi_tn, bp_corrected_phi_tn, exact_phi_onsager, niters)
end
