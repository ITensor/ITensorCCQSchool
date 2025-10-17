using ITensorMPS: MPO, OpSum, dmrg, maxlinkdim, random_mps, siteinds
# Functions for performing 2D DMRG
using ITensorMPS: square_lattice
# Functions for performing measurements of MPS
using ITensorMPS: ITensorMPS, AbstractObserver, expect, inner
# Use to set the RNG seed for reproducibility
using StableRNGs: StableRNG
# Load the Plots package for plotting
using Plots: Plots, plot

include("../src/animate.jl")

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
    psi0 = random_mps(sites, state; linkdims = 10)

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

    return (; energy, H, psi, nx, ny, nsite, U, szs, ns, nsweeps, maxdim, cutoff, noise)
end
