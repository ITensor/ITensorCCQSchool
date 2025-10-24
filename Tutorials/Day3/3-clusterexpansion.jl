using NamedGraphs.NamedGraphGenerators: named_grid
using Statistics: mean

# Load the Plots library for plotting results
using Plots: Plots, plot

include("isingtensornetwork.jl")
include("beliefpropagationfunctions.jl")

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
- `bp_phi_g::Number`: The Bethe-Peierls free energy density computed via belief propagation.
- `bp_corrected_phi_g::Number`: The loop-corrected free energy density.
- `exact_phi_onsager::Number`: The exact free energy density from Onsager's solution in the thermodynamic limit.
- `niterations::Int`: The number of iterations taken for convergence in belief propagation.
"""
function main(; beta::Number = 0.2, outputlevel::Int=1)
    g = named_grid((5,5); periodic = true)

    tensornetwork = ising_tensornetwork(g, beta)
    messages, niterations = belief_propagation(tensornetwork, g, 1000; outputlevel)

    bp_phi_g = bp_phi(tensornetwork, messages, g)
    smallest_loop_size = 4
    bp_corrected_phi_g = bp_corrected_phi(tensornetwork, messages, g, smallest_loop_size)
    exact_phi_onsager = ising_phi(beta)

    return (; bp_phi_g, bp_corrected_phi_g, exact_phi_onsager, niterations)
end
