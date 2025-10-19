using NamedGraphs: NamedEdge, NamedGraph, add_edges

using NamedGraphs.NamedGraphGenerators: named_path_graph

include("isingtensornetwork.jl")

function main(; beta::Number = 0.2, outputlevel::Int)
    # Create a simple graph
    g = NamedGraph([1,2,3])
    edges = [1 => 2, 2 => 3]
    g = add_edges(g, edges)

    tensornetwork = ising_tensornetwork(g, beta)

    z = prod(values(tensornetwork))[]

    return (; g, tensornetwork, z)
end
