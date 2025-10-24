using ITensorMPS: MPO, OpSum, dmrg, maxlinkdim, random_mps, siteinds
# Functions for performing measurements of MPS
using ITensorMPS: ITensorMPS, AbstractObserver, correlation_matrix, expect, inner
# Use to set the RNG seed for reproducibility
using StableRNGs: StableRNG
# Load the Plots package for plotting
using Plots: Plots, plot

include("../src/animate.jl")

function animate_dmrg_sz(res; fps = res.nsite)
    return animate(; nframes = length(res.szs), fps) do i
        return plot(
            res.szs[i]; xlim = (1, res.nsite), ylim = (-0.25, 0.25), xlabel = "Site j",
            ylabel = "⟨Szⱼ⟩", legend = false, title = "Sweep = $(i ÷ (2 * res.nsite) + 1)"
        )
    end
end

@kwdef struct SzObserver <: AbstractObserver
    szs::Vector{Vector{Float64}} = Vector{Float64}[]
end
function ITensorMPS.measure!(obs::SzObserver; psi, kwargs...)
    push!(obs.szs, expect(psi, "Sz"))
    return nothing
end

"""
    main(; kwargs...)

Perform DMRG on a Heisenberg spin-1/2 chain and measure ⟨Sz⟩ and ⟨SzⱼSz⟩.

# Keywords
- `nsite::Int = 30`: Number of sites in the spin chain.
- `nsweeps::Int = 5`: Number of DMRG sweeps.
- `maxdim::Vector{Int} = [10, 20, 100, 100, 200]`: Maximum bond dimensions for each sweep.
- `cutoff::Vector{Float64} = [1.0e-10]`: Cutoff values for each sweep.
- `outputlevel::Int = 1`: Controls how much information will be printed by the script.

# Results
A named tuple containing:
- `energy::Float64`: The optimized ground state energy.
- `H::MPO`: The Hamiltonian as an MPO.
- `psi::MPS`: The optimized ground state wavefunction as an MPS.
- `sz::Vector{Float64}`: Vector of ⟨Sz⟩ measurements.
- `szsz::Matrix{Float64}`: Correlation matrix of ⟨SzⱼSz⟩.
- `szs::Vector{Vector{Float64}}`: Vector of ⟨Sz⟩ measurements at each DMRG step.
- `nsite::Int`: Same as above.
- `nsweeps::Int`: Same as above.
- `maxdim::Vector{Int}`: Same as above.
- `cutoff::Vector{Float64}`: Same as above.
"""
function main(;
        # Number of sites
        nsite = 30,
        # DMRG parameters
        nsweeps = 5,
        maxdim = [10, 20, 100, 100, 200],
        cutoff = [1.0e-10],
        outputlevel = 1,
    )
    # Build the physical indices for nsite spins (spin 1/2)
    sites = siteinds("S=1/2", nsite)

    # Build the Heisenberg Hamiltonian as an MPO
    terms = OpSum()
    for j in 1:(nsite - 1)
        terms += 1 / 2, "S+", j, "S-", j + 1
        terms += 1 / 2, "S-", j, "S+", j + 1
        terms += "Sz", j, "Sz", j + 1
    end
    H = MPO(terms, sites)

    # It has bond dimension 5
    if outputlevel > 0
        println("MPO bond dimension: ", maxlinkdim(H))
    end

    # Initial state for DMRG
    psi0 = random_mps(sites; linkdims = 10)

    # It starts with a bond dimension 10
    if outputlevel > 0
        println("Initial MPS bond dimension: ", maxlinkdim(psi0))
    end

    # Run DMRG
    observer = SzObserver()
    energy, psi = dmrg(
        H, psi0; nsweeps, maxdim, cutoff, observer, outputlevel = min(outputlevel, 1)
    )
    szs = observer.szs

    if outputlevel > 0
        println("Optimized MPS bond dimension: ", maxlinkdim(psi))
        println("Energy: ", energy)
        println("⟨ψ|ψ⟩: ", inner(psi, psi))
        println("⟨ψ|H|ψ⟩: ", inner(psi', H, psi))
    end

    sz = expect(psi, "Sz")
    if outputlevel > 0
        display(
            plot(
                sz; xlim = (1, nsite), ylim = (-0.25, 0.25),
                xlabel = "Site", ylabel = "⟨Sz⟩", legend = false
            )
        )
    end

    szsz = correlation_matrix(psi, "Sz", "Sz")
    if outputlevel > 0
        display(
            plot(
                szsz[nsite ÷ 2, :]; xlim = (1, nsite), ylim = (-0.25, 0.25),
                xlabel = "Site", ylabel = "⟨SzⱼSz⟩", legend = false
            )
        )
    end

    res = (; energy, H, psi, sz, szsz, szs, nsite, nsweeps, maxdim, cutoff)
    if outputlevel > 1
        animate_dmrg_sz(res)
    end
    return res
end
