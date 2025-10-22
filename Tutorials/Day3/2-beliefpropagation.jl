using NamedGraphs.NamedGraphGenerators: named_grid
using Statistics: mean

# Load the Plots library for plotting results
using Plots: Plots, plot

include("isingtensornetwork.jl")
include("beliefpropagationfunctions.jl")

"""
    main(; Lx::Int = 3, Ly::Int = 3, beta::Number = 0.2, periodic = false, outputlevel::Int = 1)
    Creates a grid graph of size Lx by Ly, constructs the Ising tensor network on it, and computes the Bethe-Peierls free energy density using belief propagation.
    # Arguments
    - `Lx::Int`: The number of vertices in the x-direction (default is 3).
    - `Ly::Int`: The number of vertices in the y-direction (default is 3).
    - `beta::Number`: The inverse temperature parameter (default is 0.2).
    - `periodic::Bool`: Whether to use periodic boundary conditions (default is false).
    - `outputlevel::Int`: Level of output detail (default is 1).
    # Returns
    - `bp_phi_g::Number`: The Bethe-Peierls free energy density computed via belief propagation.
    - `exact_phi_onsager::Number`: The exact free energy density from Onsager's solution in the thermodynamic limit.
    - `niterations::Int`: The number of iterations taken for convergence in belief propagation.
"""
function main(; Lx::Int = 3, Ly::Int = 3, beta::Number = 0.2, periodic = false, outputlevel::Int = 1)
    g = named_grid((Lx,Ly); periodic)

    tensornetwork = ising_tensornetwork(g, beta)
    messages, niterations = belief_propagation(tensornetwork, g, 1000; outputlevel)

    bp_phi_g = bp_phi(tensornetwork, messages, g)
    exact_phi_onsager = ising_phi(beta)
    return (; bp_phi_g, exact_phi_onsager, niterations)
end
