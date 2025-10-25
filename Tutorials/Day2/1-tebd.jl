using ITensorMPS: MPO, OpSum, dmrg, maxlinkdim, random_mps, siteinds, linkinds
# Functions for performing measurements of MPS
using ITensorMPS: expect, inner, MPS, orthogonalize, svd
# Functions for time evolution
using ITensorMPS: apply, op
using LinearAlgebra: normalize, diag
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

function entanglement_entropy(ψ::MPS, bond::Int = round(Int, length(ψ) / 2); cutoff = 1e-12, kwargs...)
    ψ = normalize(ψ)
    ψ = orthogonalize(ψ, bond)
    U, S, V = svd(ψ[bond], (linkinds(ψ, bond - 1)..., siteinds(ψ, bond)...))
    Sd = Array(diag(S))
    Sd = Sd .^ 2
    return sum(d -> (d > cutoff) ? -d * log(d) : 0.0, Sd)
end

"""
   main(; kwargs...)

Perform time-evolving block decimation (TEBD) on a 1D Heisenberg spin-1/2 chain to perform
real time evolution of the Heisenberg ground state with a spin flip in the center of the
chain.

# Keywords
- `nsite::Int = 30`: Number of sites in the spin chain.
- `time::Float64 = 5.0`: Total time for evolution.
- `timestep::Float64 = 0.1`: Time step for each TEBD application.
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
        time = 6.0,
        timestep = 0.1,
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
    initial_energy, psi = dmrg(
        H, psi0; nsweeps = 5, maxdim = [10, 20, 100, 100, 200],
        cutoff = [1.0e-10], outputlevel = min(outputlevel, 1)
    )

    # Make gates (1, 2), (2, 3), (3, 4), ...
    gates = map(1:(nsite - 1)) do j
        si, sj = sites[j], sites[j + 1]
        hj = 1 / 2 * op("S+", si) * op("S-", sj) +
            1 / 2 * op("S-", si) * op("S+", sj) +
            op("Sz", si) * op("Sz", sj)
        return exp(-im * timestep / 2 * hj)
    end
    # Include gates in reverse order too
    # (N, N - 1), (N - 1, N - 2), ...
    append!(gates, reverse(gates))

    # Make starting state
    j = nsite ÷ 2
    psit = apply(op("S+", sites[j]), psi)
    psit = normalize(psit)

    szs = [expect(psit, "Sz")]
    energies = ComplexF64[inner(psit', H, psit)]
    entanglements = [entanglement_entropy(psit, nsite ÷ 2)]
    times = 0.0:timestep:time
    print_every = 1
    for current_time in times[2:end]
        psit = apply(gates, psit; cutoff)
        psit = normalize(psit)
        energy_t = inner(psit', H, psit)
        sz_t = expect(psit, "Sz")
        push!(szs, sz_t)
        push!(energies, energy_t)
        push!(entanglements, entanglement_entropy(psit, nsite ÷ 2))
        if floor(current_time - timestep + 10eps()) ≠ floor(current_time) &&
                floor(current_time) % print_every == 0
            if outputlevel > 0
                println("time: ", current_time)
                println("Bond dimension: ", maxlinkdim(psit))
                println("⟨ψₜ|Szⱼ|ψₜ⟩: ", sz_t[j])
                println("∑ⱼ⟨ψₜ|Szⱼ|ψₜ⟩: ", sum(sz_t))
                println("⟨ψₜ|H|ψₜ⟩: ", energy_t)
                println()
            end
        end
    end

    res = (; inital_energy, H, psi, times, szs, energies, entanglements, nsite, time, timestep, cutoff)
    if outputlevel > 1
        animate_tebd_sz(res)
    end
    return res
end
