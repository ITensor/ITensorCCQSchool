using ITensorMPS: MPO, OpSum, dmrg, maxlinkdim, random_mps, siteinds
# Functions for performing measurements of MPS
using ITensorMPS: expect, inner
# Functions for time evolution
using ITensorMPS: apply, op
using LinearAlgebra: normalize
# Use to set the RNG seed for reproducibility
using StableRNGs: StableRNG
# Load the Plots package for plotting
using Plots: Plots, plot

include("../src/animate.jl")

function plot_tebd_sz(res; step::Int)
    return plot(
        res.szs[step]; xlim = (1, res.nsite), ylim = (-0.5, 0.5), xlabel = "Site j",
        ylabel = "⟨Szⱼ(t=$(res.times[step]))⟩", legend = false
    )
end

function animate_tebd_sz(res; fps = res.nsite)
    return animate(i -> plot_tebd_sz(res; step = i); nframes = length(res.szs), fps)
end

"""
   main(; kwargs...)

Perform time-evolving block decimation (TEBD) on a 1D Heisenberg spin-1/2 chain to perform
imaginary time evolution to find the ground state.

# Keywords
- `nsite::Int = 30`: Number of sites in the spin chain.
- `time::Float64 = 20.0`: Total time for evolution.
- `timestep::Float64 = 0.2`: Time step for each TEBD application.
- `cutoff::Float64 = 1.0e-10`: Cutoff for truncation during TEBD.
- `outputlevel::Int = 1`: Controls how much information will be printed by the script.

# Returns
A named tuple containing:
- `energy::Float64`: The final energy after time evolution.
- `H::MPO`: The Hamiltonian as an MPO.
- `psi::MPS`: The final wavefunction after time evolution as an MPS.
- `times::Vector{Float64}`: Vector of time points at which measurements were taken.
- `szs::Vector{Vector{Float64}}`: Vector of ⟨Sz⟩ measurements at each time point.
- `energies::Vector{Float64}`: Vector of energy measurements at each time point.
- `nsite::Int`: Same as above.
- `time::Float64`: Same as above.
- `timestep::Float64`: Same as above.
- `cutoff::Float64`: Same as above.
"""
function main(;
        # Number of sites
        nsite = 30,
        # TEBD parameters
        time = 20.0,
        timestep = 0.2,
        cutoff = 1.0e-10,
        outputlevel = 1,
    )
    # Build the physical indices for nsite spins (spin 1/2)
    sites = siteinds("S=1/2", nsite)

    # Run DMRG to get starting state for time evolution
    terms = OpSum()
    for j in 1:(nsite - 1)
        terms += 1 / 2, "S+", j, "S-", j + 1
        terms += 1 / 2, "S-", j, "S+", j + 1
        terms += "Sz", j, "Sz", j + 1
    end
    H = MPO(terms, sites)
    psi0 = random_mps(sites; linkdims = 10)
    energy, psi = dmrg(
        H, psi0; nsweeps = 5, maxdim = [10, 20, 100, 100, 200],
        cutoff = [1.0e-10], outputlevel = min(outputlevel, 1)
    )

    # Make gates (1, 2), (2, 3), (3, 4), ...
    gates = map(1:(nsite - 1)) do j
        si, sj = sites[j], sites[j + 1]
        hj = 1 / 2 * op("S+", si) * op("S-", sj) +
            1 / 2 * op("S-", si) * op("S+", sj) +
            op("Sz", si) * op("Sz", sj)
        return exp(-timestep / 2 * hj)
    end
    # Include gates in reverse order too
    # (N, N - 1), (N - 1, N - 2), ...
    append!(gates, reverse(gates))

    # Make starting state
    rng = StableRNG(123)
    psit = random_mps(rng, sites)
    psit = normalize(psit)

    szs = [expect(psit, "Sz")]
    energies = [inner(psit', H, psit)]
    times = 0.0:timestep:time
    print_every = 5
    for current_time in times[2:end]
        psit = apply(gates, psit; cutoff)
        psit = normalize(psit)
        energy_t = inner(psit', H, psit)
        sz_t = expect(psit, "Sz")
        push!(szs, sz_t)
        push!(energies, energy_t)
        if floor(current_time - timestep + 10eps()) ≠ floor(current_time) &&
                floor(current_time) % print_every == 0
            if outputlevel > 0
                println("time: ", current_time)
                println("Bond dimension: ", maxlinkdim(psit))
                println("⟨ψₜ|Szⱼ|ψₜ⟩: ", sz_t[nsite ÷ 2])
                println("∑ⱼ⟨ψₜ|Szⱼ|ψₜ⟩: ", sum(sz_t))
                println("⟨ψₜ|H|ψₜ⟩: ", energy_t)
                println()
            end
        end
    end

    res = (; energy, H, psi, times, szs, energies, nsite, time, timestep, cutoff)
    if outputlevel > 1
        animate_tebd_sz(res)
    end
    return res
end
