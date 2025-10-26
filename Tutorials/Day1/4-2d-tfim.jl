using ITensors: flux, pause
using ITensorMPS: MPO, OpSum, dmrg, maxlinkdim, random_mps, siteinds
# Functions for performing 2D DMRG
using ITensorMPS: square_lattice
# Functions for performing measurements of MPS
using ITensorMPS: ITensorMPS, AbstractObserver, expect, inner
# Use to set the RNG seed for reproducibility
using StableRNGs: StableRNG
# Load the Plots package for plotting
using Plots: Plots, plot, scatter, quiver, quiver!
using Random: shuffle

include("../src/animate.jl")

function plot_spins(res, i::Int = length(res.szs))
    points = vec(reverse.(Tuple.(CartesianIndices((res.ny, res.nx)))))
    xs = first.(points)
    ys = last.(points)
    s = vec(res.szs[i])
    quiver_data = (zeros(length(xs)), s)
    plt = quiver(xs, ys, quiver=quiver_data; color = :blue, arrow = :closed)
    return plt
end

function animate_spins(res; fps = res.nsite)
    return animate(i -> plot_spins(res, i); nframes = length(res.szs), fps)
end

@kwdef struct SzObserver <: AbstractObserver
    nx::Int
    ny::Int
    szs::Vector{Matrix{Float64}} = Vector{Float64}[]
end
function ITensorMPS.measure!(obs::SzObserver; psi, kwargs...)
    #println("Pushing to obs.szs"); pause()
    push!(obs.szs, reshape(expect(psi, "Sz"), (obs.ny, obs.nx)))
    return nothing
end


"""
    main(; kwargs...)
    
Perform DMRG on a 2D transverse-field Ising model (TFIM) and measure ⟨Sz⟩ and site densities.
    
# Keywords
- `nx::Int = 4`: Number of sites in x direction.
- `ny::Int = 2`: Number of sites in y direction.
- `h_max::Float64 = 5.0`: maximum on-site transverse field (reached at right-hand side of system)
- `ramp_width::Float64 = 1.0`: width of ramp function for magnetic field
- `nsweeps::Int = 5`: Number of DMRG sweeps.
- `maxdim::Vector{Int} = [100, 200, 400, 800, 1600]`: Maximum bond dimensions for each sweep.
- `cutoff::Vector{Float64} = [1.0e-6]`: Cutoff values for each sweep.
- `noise::Vector{Float64} = [1.0e-6, 1.0e-7, 1.0e-8, 0.0]`: Noise values for each sweep.
- `outputlevel::Int = 1`: Controls how much information will be printed by the script.

# Returns
A named tuple containing:
- `energy::Float64`: The optimized ground state energy.
- `H::MPO`: The Hamiltonian as an MPO.
- `psi::MPS`: The ground state wavefunction as an MPS.
- `nx::Int`: Same as above.
- `ny::Int`: Same as above.
- `nsite::Int`: Total number of sites.
- `szs::Vector{Vector{Float64}}`: Vector of ⟨Sz⟩ measurements at each DMRG step.
- `nsweeps::Int`: Same as above.
- `maxdim::Vector{Int}`: Same as above.
- `cutoff::Vector{Float64}`: Same as above.
- `noise::Vector{Float64}`: Same as above.
"""
function main(;
        # Number of sites in x and y
        nx = 4,
        ny = 2,
        # Model parameters
        h_max = 5.0,
        ramp_width = 1.0,
        # DMRG parameters
        nsweeps = 5,
        maxdim = [100, 200, 400, 800, 1600],
        cutoff = [1.0e-6],
        noise = [1.0e-6, 1.0e-7, 1.0e-8, 0.0],
        outputlevel = 1,
    )
    nsite = nx * ny

    field(j) = h_max*(1/2+1/2*tanh((j-nx÷2)/ramp_width))

    # Build the physical indices
    sites = siteinds("S=1/2", nsite)

    # Build the lattice
    lattice = square_lattice(nx, ny; yperiodic = true)

    # Build the transverse-field Ising model as an MPO
    terms = OpSum()
    for b in lattice
        i, j = b.s1, b.s2
        xi, xj = Int(b.x1), Int(b.x2)
        terms += -1,"Z", i, "Z", j
        terms += field(xi)/2, "X", i
        terms += field(xj)/2, "X", j
    end
    H = MPO(terms, sites)

    if outputlevel > 0
        println("MPO bond dimension: ", maxlinkdim(H))
    end

    # Initial state for DMRG
    state = [iseven(j) ? "↑" : "↓" for j=1:nsite]
    rng = StableRNG(123)
    psi0 = random_mps(rng, sites, state; linkdims = 10)

    # It starts with a bond dimension 10
    if outputlevel > 0
        println("Initial MPS bond dimension: ", maxlinkdim(psi0))
    end

    # Run DMRG, measuring Sz and site densities during the sweep
    observer = SzObserver(; nx, ny)
    energy, psi = dmrg(
        H, psi0; nsweeps, maxdim, cutoff, observer, outputlevel = min(outputlevel, 1)
    )
    szs = observer.szs

    res = (; energy, H, psi, nx, ny, nsite, szs, nsweeps, maxdim, cutoff, noise)
    if outputlevel > 1
        animate_spins(res)
    end
    return res
end
