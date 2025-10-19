using NamedGraphs: NamedEdge, NamedGraph, vertices, edges, src, dst, neighbors, simplecycles_limited_length
using NamedGraphs.GraphsExtensions: all_edges, boundary_edges
using LinearAlgebra: normalize, dot
using ITensors: Index, ITensor, commonind, onehot

function updated_message(tensornetwork::Dict, messages::Dict, g::NamedGraph, e::NamedEdge)
    incoming_es = setdiff(boundary_edges(g, [src(e)]), [reverse(e)])
    local_tensor = tensornetwork[src(e)]
    messages = copy(messages)
    if isempty(incoming_vs) 
        m = normalize(local_tensor)
    else
        incoming_messages = [messages[e] for e in incoming_es]
        m = normalize(local_tensor * prod(incoming_messages))
    end
    return m
end

function update_messages(tensornetwork::Dict, messages::Dict, g::NamedGraph)
    updated_messages = copy(messages)
    for e in all_edges(g)
        updated_messages[e] = updated_message(tensornetwork, messages, g, e)
    end
    return updated_messages
end

function initial_messages(tn::Dict, g::NamedGraph)
    edges = all_edges(g)
    return Dict(e => normalize(onehot(commonind(tn[src(e)], tn[dst(e)]) => 1)) for e in edges)
end

function binormalized_messages(messages::Dict, g::NamedGraph)
    binorm_messages = copy(messages)
    for e in edges(g)
        n = dot(messages[e], messages[reverse(e)])
        if sign(n) < 0
            n = abs(n)
            s = -1
        else
            s = 1
        end
        binorm_messages[e] = s * messages[e] / sqrt(n)
        binorm_messages[reverse(e)] = messages[reverse(e)] / sqrt(n)
    end
    return binorm_messages
end

"""
    bp_phi(tn::Dict, messages::Dict, g::NamedGraph)
    Computes the Bethe-Peierls free energy density.
    # Arguments
    - `tn::Dict`: A dictionary representing the tensor network, where keys are vertices and values are tensors.
    - `messages::Dict`: A dictionary containing the messages for each edge in the graph.
    - `g::NamedGraph`: The named graph representing the structure of the tensor network.
    # Returns
    - `Number`: The free energy estimate per vertex.
"""
function bp_phi(tn::Dict, messages::Dict, g::NamedGraph)
    messages = binormalized_messages(messages, g)
    return sum([log(local_factor(tn, messages, g, v)[]) for v in vertices(g)]) / length(vertices(g))
end

"""
    belief_propagation(tn::Dict, g::NamedGraph, niters::Int; tol::Float64=1e-10, outputlevel::Int = 1)
    Performs the Belief Propagation algorithm on a given tensor network defined over a named graph.
    # Arguments
    - `tn::Dict`: A dictionary representing the tensor network, where keys are vertices and values are tensors.
    - `g::NamedGraph`: The named graph representing the structure of the tensor network.
    - `niters::Int`: The maximum number of iterations to perform.
    # Keyword Arguments
    - `tol::Float64=1e-10`: The tolerance for convergence.
    - `outputlevel::Int = 1`: The verbosity level of the output.
    # Returns
    - `messages::Dict`: A dictionary containing the converged messages for each edge in the graph.
    - `niterations::Int`: The number of iterations taken to converge, or `Inf` if not converged within `niters`.
"""
function belief_propagation(tn::Dict, g::NamedGraph, niters::Int; tol::Float64=1e-10, outputlevel::Int = 1)
    edges = all_edges(g)
    messages = initial_messages(tn, g)
    for i in 1:niters   
        old_messages = copy(messages)
        messages = update_messages(tn, messages, g)
        if mean([1 - dot(messages[e], old_messages[e])^2 for e in edges]) < tol
            outputlevel >= 1 && println("BP Algorithm Converged after $i iterations")
            return messages, i
            break
        end

    end
    return messages, Inf
end

function local_factor(tn::Dict, messages::Dict, g::NamedGraph, v)
    incoming_messages = [messages[e] for e in boundary_edges(g, [v])]
    m = prod([[tn[v]]; incoming_messages])
    return m
end

"""
    bp_corrected_phi(tn::Dict, messages::Dict, g::NamedGraph, smallest_loop_size::Int)
    Computes the first order corrected Bethe-Peierls free energy.
    # Arguments
    - `tn::Dict`: A dictionary representing the tensor network, where keys are vertices and values are tensors.
    - `messages::Dict`: A dictionary containing the messages for each edge in the graph.
    - `g::NamedGraph`: The named graph representing the structure of the tensor network.
    - `smallest_loop_size::Int`: The size of the smallest loops to consider for the correction.
    # Returns
    - `Number`: The corrected free energy estimate per vertex.
"""
function bp_corrected_phi(tn::Dict, messages::Dict, g::NamedGraph, smallest_loop_size::Int)
    messages = binormalized_messages(messages, g)
    bp_phi_g = bp_phi(tn, messages, g)
    rescaled_tn = Dict(v => tn[v] / local_factor(tn, messages, g, v) for v in vertices(g))
    cycles = filter(c -> length(c) >  2, simplecycles_limited_length(g, smallest_loop_size))
    cycles = unique(Set.(cycles))
    isempty(cycles) && error("No cycles found with length $smallest_loop_size")
    cycle_weights = []
    for cycle in cycles
        incoming_messages = [messages[e] for e in boundary_edges(g, cycle)]
        local_tensors = [rescaled_tn[v] for v in cycle]
        weight = prod([local_tensors; incoming_messages])[]
        push!(cycle_weights, weight)
    end
    return sum(log.(cycle_weights)) / length(vertices(g)) + bp_phi_g
end