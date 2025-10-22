using NamedGraphs: NamedEdge, NamedGraph, add_edges

using NamedGraphs.NamedGraphGenerators: named_path_graph

include("isingtensornetwork.jl")

"""
    main(; beta::Number = 0.2, outputlevel::Int=1)
    Creates a simple path graph with 3 vertices, constructs the Ising tensor network on it, and computes the partition function Z.
    # Arguments
    - `beta::Number`: The inverse temperature parameter (default is 0.2).
    - `outputlevel::Int`: Level of output detail (default is 1).
    # Returns
    - `g::NamedGraph`: The created graph.
    - `tensornetwork::Dict{Any, ITensor}`: The tensor network representation of the Ising model on the graph.
    - `z::Number`: The computed partition function Z.
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
