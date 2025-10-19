using NamedGraphs.NamedGraphGenerators: named_grid
using Statistics: mean
using QuadGK: quadgk

# Load the Plots library for plotting results
using Plots: Plots, plot

include("isingtensornetwork.jl")
include("beliefpropagationfunctions.jl")

# -β f as a function of β
function ising_phi(β)
    g(θ1, θ2) = log(
        cosh(2β)*cosh(2β) -
        sinh(2β)*cos(θ1) -
        sinh(2β)*cos(θ2)
    )
    inner(θ2) = quadgk(θ1 -> g(θ1, θ2), 0, 2π)[1]
    return -log(2) + (1/(8π^2)) * quadgk(inner, 0, 2π)[1]
end

function main(; beta::Number = 0.2, outputlevel::Int=1)
    g = named_grid((5,5); periodic = true)

    tensornetwork = ising_tensornetwork(g, beta)
    messages, niterations = belief_propagation(tensornetwork, g, 100; outputlevel)

    bp_phi_g = bp_phi(tensornetwork, messages, g)
    smallest_loop_size = 4
    bp_corrected_phi_g = bp_corrected_phi(tensornetwork, messages, g, smallest_loop_size)
    exact_phi_onsager = ising_phi(beta)

    return (; bp_phi_g, bp_corrected_phi_g, exact_phi_onsager, niterations)
end
