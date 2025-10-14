using ITensorMPS: MPO, OpSum, correlation_matrix, dmrg, expect, inner, maxlinkdim,
    random_mps, siteinds
using Plots: Plots, plot

# Load the UnicodePlots backend for plotting to the terminal.
Plots.unicodeplots()

let
    # Build the physical indices for 30 spins (spin 1/2)
    N = 30
    sites = siteinds("S=1/2", N)

    # Build the Heisenberg Hamiltonian as an MPO
    os = OpSum()
    for j in 1:(N - 1)
        os += 1 / 2, "S+", j, "S-", j + 1
        os += 1 / 2, "S-", j, "S+", j + 1
        os += "Sz", j, "Sz", j + 1
    end
    H = MPO(os, sites)

    # It has bond dimension 5
    println("MPO bond dimension is $(maxlinkdim(H))")

    # Initial state for DMRG
    psi0 = random_mps(sites; linkdims = 10)

    # It starts with a bond dimension 10
    println("Initial MPS bond dimension is $(maxlinkdim(psi0))")

    # DMRG Parameters
    nsweeps = 5
    maxdim = [10, 20, 100, 100, 200]
    cutoff = [1.0e-10]

    # Run DMRG
    energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff)

    println("Optimized MPS bond dimension is $(maxlinkdim(psi))")

    @show inner(psi, psi)
    @show inner(psi', H, psi)

    sz = expect(psi, "Sz")
    display(plot(sz; xlabel = "Site", ylabel = "⟨Sz⟩"))

    sz_sz = correlation_matrix(psi, "Sz", "Sz")
    display(plot(sz_sz[15, :]; xlabel = "Site", ylabel = "⟨Sz₁₅Szⱼ⟩"))

    return (; psi, sz, sz_sz)
end
