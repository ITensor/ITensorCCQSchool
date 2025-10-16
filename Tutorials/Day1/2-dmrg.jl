using ITensorMPS: MPO, OpSum, dmrg, maxlinkdim, random_mps, siteinds
# Use to set the RNG seed for reproducibility
using StableRNGs: StableRNG
# Load the Plots library for plotting results
using Plots: Plots, plot
# Set the plotting backend to UnicodePlots to print the plots in the terminal
Plots.unicodeplots()

"""
    main(; kwargs...)

Run the Density Matrix Renormalization Group (DMRG) algorithm to compute the ground state
energy and wavefunction of a 1D Heisenberg spin-1/2 chain.

This function constructs the spin chain Hamiltonian as a Matrix Product Operator (MPO),
initializes a random Matrix Product State (MPS) as the starting point, and performs a series
of DMRG sweeps to optimize the ground state energy and wavefunction.

# Keyword Arguments
- `nsite::Int = 30`: Number of sites in the spin chain.
- `nsweeps::Int = 5`: Number of DMRG sweeps to perform.
- `maxdim::Vector{Int} = [10, 20, 100, 100, 200]`: Maximum bond dimensions for each sweep.
- `cutoff::Vector{Float64} = [1.0e-10]`: Cutoff for truncation.
- `outputlevel::Int = 1`: Controls how much information will be printed by the script.

# Outputs
A named tuple containing:
- `energy`: The optimized ground state energy.
- `H`: The Hamiltonian as an MPO.
- `psi`: The optimized ground state wavefunction as an MPS.
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

    outputlevel > 0 && println("nsite: ", nsite)
    outputlevel > 0 && println("nsweeps: ", nsweeps)
    outputlevel > 0 && println("maxdim: ", maxdim)
    outputlevel > 0 && println("cutoff: ", cutoff)

    # Build the physical indices for nsite spins (spin 1/2)
    sites = siteinds("S=1/2", nsite)

    # Build the Heisenberg Hamiltonian as an MPO
    os = OpSum()
    for j in 1:(nsite - 1)
        os += 1 / 2, "S+", j, "S-", j + 1
        os += 1 / 2, "S-", j, "S+", j + 1
        os += "Sz", j, "Sz", j + 1
    end
    H = MPO(os, sites)

    # It has bond dimension 5
    outputlevel > 0 && println("MPO bond dimension: ", maxlinkdim(H))

    # Initial state for DMRG
    rng = StableRNG(123)
    psi0 = random_mps(rng, sites; linkdims = 10)

    # It starts with a bond dimension 10
    outputlevel > 0 && println("Initial MPS bond dimension: ", maxlinkdim(psi0))

    # Run DMRG
    energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff, outputlevel)

    outputlevel > 0 && println("Optimized MPS bond dimension: ", maxlinkdim(psi))
    outputlevel > 0 && println("DMRG energy: ", energy)

    return (; energy, H, psi, nsite, nsweeps, maxdim, cutoff)
end
