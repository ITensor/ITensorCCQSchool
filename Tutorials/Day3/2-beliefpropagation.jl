using NamedGraphs: NamedEdge, NamedGraph, vertices, edges, src, dst, neighbors
using LinearAlgebra: normalize, dot
using NamedGraphs.NamedGraphGenerators: named_grid
using ITensors: Index, ITensor, prime, commonind, onehot
using Statistics: mean
using QuadGK: quadgk

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

function calculate_bp_phi(tn::Dict, messages::Dict, g::NamedGraph)
    f_node = 0.0
    for v in vertices(g)
        incoming_messages = [messages[NamedEdge(vn => v)] for vn in neighbors(g, v)]
        m = prod([[tn[v]]; incoming_messages])
        f_node += log(m[])
    end
    f_edge = 0.0
    for e in edges(g)
        m1 = messages[e]
        m2 = messages[reverse(e)]
        f_edge += log((m1 * m2)[])
    end
    return (f_node - f_edge) / length(vertices(g))
end

function main(; Lx::Int, Ly::Int, beta::Number = 0.2, periodic = false)
    g = named_grid((Lx,Ly); periodic)

    tensornetwork = ising_tensornetwork(g, beta)
    messages, niterations = belief_propagation(tensornetwork, g, 100)

    bp_phi = calculate_bp_phi(tensornetwork, messages, g)
    exact_phi_onsager = ising_phi(beta)
    return (; bp_phi, exact_phi_onsager, niterations)
end
