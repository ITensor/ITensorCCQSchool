using NamedGraphs.NamedGraphGenerators: named_grid
using Statistics: mean

# Load the Plots library for plotting results
using Plots: Plots, plot

include("isingtensornetwork.jl")
include("beliefpropagationfunctions.jl")

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
