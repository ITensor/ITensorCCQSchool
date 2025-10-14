using ITensorMPS: MPO, OpSum, dmrg, maxlinkdim, random_mps, siteinds

function main(;
        # Number of sites
        N = 30,
        # DMRG parameters
        nsweeps = 5,
        maxdim = [10, 20, 100, 100, 200],
        cutoff = [1.0e-10],
        outputlevel = 1,
    )
    # Build the physical indices for N spins (spin 1/2)
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
    outputlevel > 0 && println("MPO bond dimension: ", maxlinkdim(H))

    # Initial state for DMRG
    psi0 = random_mps(sites; linkdims = 10)

    # It starts with a bond dimension 10
    outputlevel > 0 && println("Initial MPS bond dimension: ", maxlinkdim(psi0))

    # Run DMRG
    energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff, outputlevel)

    outputlevel > 0 && println("Optimized MPS bond dimension: ", maxlinkdim(psi))
    outputlevel > 0 && println("Energy: ", energy)

    return (; energy, H, psi)
end
