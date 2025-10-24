using NamedGraphs: NamedEdge, NamedGraph, add_edges

using NamedGraphs.NamedGraphGenerators: named_path_graph

include("isingtensornetwork.jl")

"""
    main(; kwargs...)
    
Creates a simple path graph, constructs the Ising tensor network on it, and computes the
partition function Z by contracting the tensor network.

# Keywords
- `beta::Number = 0.2`: The inverse temperature parameter.
- `outputlevel::Int = 1`: Controls how much information will be printed by the script.

# Returns
A named tuple containing:
- `g::NamedGraph`: The created graph.
- `tensornetwork::Dict{Any, ITensor}`: The tensor network representation of the Ising model on the graph.
- `z::Number`: The partition function Z computed by contracting the tensor network.
"""
function main(; beta::Number = 0.2, outputlevel::Int=1)
    # Create a simple graph
    g = NamedGraph([1,2,3])
    edges = [1 => 2, 2 => 3]
    g = add_edges(g, edges)

    tensornetwork = ising_tensornetwork(g, beta)

    z = prod(values(tensornetwork))[]

    return (; g, tensornetwork, z)
end
