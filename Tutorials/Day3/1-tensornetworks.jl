using NamedGraphs: NamedEdge, NamedGraph, edges, vertices, add_edges, dst, src, neighbors
using LinearAlgebra: normalize
using ITensors: ITensor, Index, delta, prime, apply

using NamedGraphs.NamedGraphGenerators: named_grid, named_hexagonal_lattice_graph, named_path_graph

"""
    ising_tn(g::NamedGraph, β::Real)
    Construct the tensor network representation of the Ising model on graph `g` at inverse temperature `β`.
    Args:
        g: A `NamedGraph` representing the lattice or graph structure.
        β: Inverse temperature (1/kT).
    Returns:
        A dictionary mapping each vertex in `g` to its corresponding tensor in the tensor network.
    The tensor network is constructed such that contracting all tensors yields the partition function of the Ising model on the graph.
"""
function ising_tn(g::NamedGraph, β::Real)
    nv = length(vertices(g))
    link = Dict(e => Index(2, "e$(src(e))_$(dst(e))") for e in edges(g))

    # symmetric sqrt of Boltzmann matrix W = exp(β σσ')
    λ1, λ2 = cosh(β), sinh(β)
    α = 0.5*(sqrt(λ1) + sqrt(λ2))
    ϕ = 0.5*(sqrt(λ1) - sqrt(λ2))
    sqrt_W = [α ϕ; ϕ α]

    T = Dict()
    for v in vertices(g)
        inds = [link[e] for e in edges(g) if src(e)==v || dst(e)==v]
        T[v] = delta(inds) *(2^(1/length(vertices(g))))
        for vn in neighbors(g, v)
            e = NamedEdge(v, vn) ∈ edges(g) ? NamedEdge(v, vn) : NamedEdge(vn, v)
            T[v] = apply(T[v], ITensor(sqrt_W, link[e], prime(link[e]))) 
        end
    end
    return T
end

function contract_tensornetwork(tensornetwork::Dict)
    return prod(ITensor[tensornetwork[v] for v in keys(tensornetwork)])[]
end

function main(; β::Number = 0.2)
    # Create a simple graph
    g = NamedGraph([1,2,3])
    edges = [1 => 2, 2 => 3]
    g = add_edges(g, edges)

    tensornetwork = ising_tn(g, β)

    z = contract_tensornetwork(tensornetwork)

    return (; g, tensornetwork, z)
end