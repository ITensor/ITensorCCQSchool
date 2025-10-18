
using NamedGraphs: NamedGraph, edges, vertices, dst, src, neighbors
using LinearAlgebra
using ITensors: ITensor, Index, delta, prime, apply

function ising_tensornetwork(g::NamedGraph, β::Real)
    nv = length(vertices(g))
    links = Dict(e => Index(2, "e$(src(e))_$(dst(e))") for e in edges(g))
    links = merge(links, Dict(reverse(e) => links[e] for e in edges(g)))

    # symmetric sqrt of Boltzmann matrix W = exp(β σσ')
    λ1, λ2 = cosh(β), sinh(β)
    α = 0.5*(sqrt(λ1) + sqrt(λ2))
    ϕ = 0.5*(sqrt(λ1) - sqrt(λ2))
    sqrt_W = [α ϕ; ϕ α]

    T = Dict()
    for v in vertices(g)
        inds = [links[e] for e in edges(g) if src(e)==v || dst(e)==v]
        T[v] = delta(inds)
        for vn in neighbors(g, v)
            e = NamedEdge(v, vn)
            T[v] = apply(T[v], ITensor(sqrt_W, links[e], prime(links[e]))) 
        end
    end
    return T
end
