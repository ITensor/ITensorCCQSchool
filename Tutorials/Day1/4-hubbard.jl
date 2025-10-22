using ITensorMPS: MPO, OpSum, dmrg, maxlinkdim, random_mps, siteinds
# Functions for performing 2D DMRG
using ITensorMPS: square_lattice
# Functions for performing measurements of MPS
using ITensorMPS: ITensorMPS, AbstractObserver, expect, inner
# Use to set the RNG seed for reproducibility
using StableRNGs: StableRNG
# Load the Plots package for plotting
using Plots: Plots, plot, scatter, quiver!

include("../src/animate.jl")

function plot_hubbard(res, i::Int = length(res.ns))
    points = vec(reverse.(Tuple.(CartesianIndices((res.ny, res.nx)))))
    xs = first.(points)
    ys = last.(points)
    n = vec(res.ns[i])
    markersize = 10 .* n ./ maximum(n)
    plt = scatter(
        xs, ys; xlims = (0, res.nx + 1), ylims = (0, res.ny + 1), markersize, markershape = :circle,
        markercolor = :lightgreen, markeralpha = 0.9, label = "⟨n⟩"
    )
    s = vec(res.szs[i])
    quiver = (zeros(length(xs)), s)
    quiver!(plt, xs, ys; quiver, color = :blue, arrow = :closed)
    return plt
end

function animate_hubbard(res; fps = res.nsite)
    return animate(i -> plot_hubbard(res, i); nframes = length(res.ns), fps)
end

@kwdef struct SzObserver <: AbstractObserver
    szs::Vector{Vector{Float64}} = Vector{Float64}[]
end
function ITensorMPS.measure!(obs::SzObserver; psi, kwargs...)
    push!(obs.szs, expect(psi, "Sz"))
    return nothing
end

@kwdef struct SzAndDensityObserver <: AbstractObserver
    nx::Int
    ny::Int
    szs::Vector{Matrix{Float64}} = Vector{Float64}[]
    ns::Vector{Matrix{Float64}} = Vector{Float64}[]
end
function ITensorMPS.measure!(obs::SzAndDensityObserver; psi, kwargs...)
    push!(obs.szs, reshape(expect(psi, "Sz"), (obs.ny, obs.nx)))
    push!(obs.ns, reshape(expect(psi, "n↑↓"), (obs.ny, obs.nx)))
    return nothing
end

"""
    main(; kwargs...)
    
Perform DMRG on a 2D Hubbard model and measure ⟨Sz⟩ and site densities.
    
# Keyword Arguments
- `nx::Int=4`: Number of sites in x direction.
- `ny::Int=2`: Number of sites in y direction.
- `U::Float64=8.0`: On-site interaction strength.
- `t::Float64=1.0`: Hopping parameter.
- `nsweeps::Int=5`: Number of DMRG sweeps.
- `maxdim::Vector{Int}=[100,200,400,800,1600]`: Maximum bond dimensions for each sweep.
- `cutoff::Vector{Float64}=[1.0e-6]`: Cutoff values for each sweep.
- `noise::Vector{Float64}=[1.0e-6,1.0e-7,1.0e-8,0.0]`: Noise values for each sweep.
- `outputlevel::Int=1`: Level of output detail.

# Outputs
A named tuple containing:
- `energy` the optimized ground state energy.
- `H` the Hamiltonian as an MPO.
- `psi` the ground state wavefunction as an MPS.
- `nx::Int`: Same as above.
- `ny::Int`: Same as above.
- `nsite::Int`: Total number of sites.
- `U::Float64`: Same as above.
- `szs`: Vector of ⟨Sz⟩ measurements at each DMRG step.
- `ns`: Vector of site density measurements at each DMRG step.
- `nsweeps::Int`: Same as above.
- `maxdim::Vector{Int}`: Same as above.
- `cutoff::Vector{Float64}`: Same as above.
- `noise::Vector{Float64}`: Same as above.
"""
function main(;
        # Number of sites in x and y
        nx = 4,
        ny = 2,
        # Hubbard parameters
        U = 8.0,
        t = 1.0,
        # DMRG parameters
        nsweeps = 5,
        maxdim = [100, 200, 400, 800, 1600],
        cutoff = [1.0e-6],
        noise = [1.0e-6, 1.0e-7, 1.0e-8, 0.0],
        outputlevel = 1,
    )
    nsite = nx * ny

    # Build the physical indices
    sites = siteinds("Electron", nsite; conserve_qns = true)

    # Build the lattice
    lattice = square_lattice(nx, ny; yperiodic = true)

    # Build the Hubbard Hamiltonian as an MPO
    terms = OpSum()
    for b in lattice
        i, j = b.s1, b.s2
        terms -= t, "c†↑", i, "c↑", j
        terms -= t, "c†↑", j, "c↑", i
        terms -= t, "c†↓", i, "c↓", j
        terms -= t, "c†↓", j, "c↓", i
    end
    for j in 1:nsite
        terms += U, "n↑↓", j
    end
    H = MPO(terms, sites)

    if outputlevel > 0
        println("MPO bond dimension: ", maxlinkdim(H))
    end

    # Initial (half-filled) state for DMRG
    state = [isodd(j) ? "↑" : "↓" for j in 1:nsite]
    rng = StableRNG(123)
    psi0 = random_mps(rng, sites, state; linkdims = 10)

    # It starts with a bond dimension 10
    if outputlevel > 0
        println("Initial MPS bond dimension: ", maxlinkdim(psi0))
    end

    # Run DMRG, measuring Sz and site densities during the sweep
    observer = SzAndDensityObserver(; nx, ny)
    energy, psi = dmrg(
        H, psi0; nsweeps, maxdim, cutoff, observer, outputlevel = min(outputlevel, 1)
    )
    szs = observer.szs
    ns = observer.ns

    res = (; energy, H, psi, nx, ny, nsite, U, szs, ns, nsweeps, maxdim, cutoff, noise)
    if outputlevel > 1
        animate_hubbard(res)
    end
    return res
end
