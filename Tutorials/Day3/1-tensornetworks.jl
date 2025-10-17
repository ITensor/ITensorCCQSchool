using NamedGraphs: NamedEdge, NamedGraph, edges, vertices, add_edges, dst, src, neighbors
using LinearAlgebra: normalize
using ITensors: ITensor, Index, delta, prime, apply

using NamedGraphs.NamedGraphGenerators: named_grid, named_hexagonal_lattice_graph, named_path_graph

"""
    ising_tensornetwork(g::NamedGraph, β::Real)
    Construct bond dimension 2 tensor network representation of the partition function of the Ising model on graph `g` at inverse temperature `β`.
    Args:
        g: A `NamedGraph` representing the lattice or graph structure.
        β: Inverse temperature (1/kT).
    Returns:
        A dictionary mapping each vertex in `g` to its corresponding tensor in the tensor network.
"""
function ising_tensornetwork(g::NamedGraph, β::Real)
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
        T[v] = delta(inds)
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

function main(; beta::Number = 0.2)
    # Create a simple graph
    g = NamedGraph([1,2,3])
    edges = [1 => 2, 2 => 3]
    g = add_edges(g, edges)

    tensornetwork = ising_tensornetwork(g, beta)

    z = contract_tensornetwork(tensornetwork)

    return (; g, tensornetwork, z)
end
