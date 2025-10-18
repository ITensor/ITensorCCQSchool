using NamedGraphs: NamedEdge, NamedGraph, vertices, edges, src, dst, neighbors, simplecycles_limited_length
using LinearAlgebra: normalize, dot
using NamedGraphs.NamedGraphGenerators: named_grid
using ITensors: Index, ITensor, prime, commonind, onehot
using Statistics: mean
using QuadGK: quadgk

# Load the Plots library for plotting results
using Plots: Plots, plot

include("isingtensornetwork.jl")

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

function updated_message(tensornetwork::Dict, messages::Dict, g::NamedGraph, e::NamedEdge)
    incoming_vs = filter(v -> v ≠ dst(e), neighbors(g, src(e)))
    local_tensor = tensornetwork[src(e)]
    messages = copy(messages)
    if isempty(incoming_vs) 
        m = normalize(local_tensor)
    else
        incoming_messages = [messages[NamedEdge(vn, src(e))] for vn in incoming_vs]
        m = normalize(local_tensor * prod(incoming_messages))
    end
    return m
end

function update_messages(tensornetwork::Dict, messages::Dict, g::NamedGraph)
    updated_messages = copy(messages)
    for e in edges(g)
        me = updated_message(tensornetwork, messages, g, e)
        mer = updated_message(tensornetwork, messages, g, reverse(e))
        updated_messages[e], updated_messages[reverse(e)] = me, mer
    end
    return updated_messages
end

function belief_propagation(tn::Dict, g::NamedGraph, niters::Int; tol::Float64=1e-10)
    all_edges = vcat(edges(g), reverse.(edges(g)))
    messages = Dict(e => normalize(onehot(commonind(tn[src(e)], tn[dst(e)]) => 1)) for e in all_edges)
    for i in 1:niters   
        old_messages = copy(messages)
        messages = update_messages(tn, messages, g)
        if mean([1 - dot(messages[e], old_messages[e])^2 for e in all_edges]) < tol
            println("BP Algorithm Converged after $i iterations")
            return messages, i
            break
        end

    end
    return messages, Inf
end

function local_factor(tn::Dict, messages::Dict, g::NamedGraph, v)
    incoming_vs = neighbors(g, v)
    incoming_messages = [messages[NamedEdge(vn => v)] for vn in incoming_vs]
    m = prod([[tn[v]]; incoming_messages])
    return m
end

function _bp_phi(tn::Dict, messages::Dict, g::NamedGraph)
    return sum([log(local_factor(tn, messages, g, v)[]) for v in vertices(g)]) / length(vertices(g))
end

function first_order_cluster_expansion_phi(tn::Dict, messages::Dict, g::NamedGraph, smallest_loop_size::Int)
    bp_phi = _bp_phi(tn, messages, g)
    rescaled_tn = Dict(v => tn[v] / local_factor(tn, messages, g, v) for v in vertices(g))
    cycles = filter(c -> length(c) >  2, simplecycles_limited_length(g, smallest_loop_size))
    cycles = unique(Set.(cycles))
    isempty(cycles) && error("No cycles found with length $smallest_loop_size")
    cycle_weights = []
    for cycle in cycles
        incoming_es = [vn => v for v in cycle for vn in neighbors(g, v) if vn ∉ cycle]
        incoming_messages = [messages[NamedEdge(e)] for e in incoming_es]
        local_tensors = [rescaled_tn[v] for v in cycle]
        weight = prod([local_tensors; incoming_messages])[]
        push!(cycle_weights, weight)
    end
    return sum(log.(cycle_weights)) / length(vertices(g)) + bp_phi
end

function main(; beta::Number = 0.2)
    g = named_grid((5,5); periodic = true)

    tensornetwork = ising_tensornetwork(g, beta)
    messages, niterations = belief_propagation(tensornetwork, g, 100)

    bp_phi = _bp_phi(tensornetwork, messages, g)
    smallest_loop_size = 4
    bp_corrected_phi = first_order_cluster_expansion_phi(tensornetwork, messages, g, smallest_loop_size)
    exact_phi_onsager = ising_phi(beta)

    bp_err = abs(bp_phi - exact_phi_onsager)
    bp_corrected_err = abs(bp_corrected_phi - exact_phi_onsager)

    return (; bp_phi, bp_err, bp_corrected_err, exact_phi_onsager, bp_corrected_phi, niterations)
end
