using NamedGraphs: NamedEdge, NamedGraph, add_edges
using NamedGraphs.NamedGraphGenerators: named_path_graph

include("ising_tensornetwork.jl")

"""
    main(; kwargs...)
    
Creates a simple path graph, constructs the Ising tensor network on it, and computes the
partition function Z by contracting the tensor network.

# Keywords
- `beta::Number = 0.2`: The inverse temperature parameter.
- `outputlevel::Int = 1`: Controls how much information will be printed by the script.

# Returns
A named tuple containing:
- `graph::NamedGraph`: The created graph.
- `tensornetwork::Dict{Any, ITensor}`: The tensor network representation of the Ising model on the graph.
- `z::Number`: The partition function Z computed by contracting the tensor network.
"""
function main(; beta::Number = 0.2, outputlevel::Int = 1)
    # Create a simple graph
    graph = NamedGraph([1, 2, 3])
    edges = [1 => 2, 2 => 3]
    graph = add_edges(graph, edges)

    tensornetwork = ising_tensornetwork(graph, beta)

    z = prod(values(tensornetwork))[]

    return (; graph, tensornetwork, z)
end
