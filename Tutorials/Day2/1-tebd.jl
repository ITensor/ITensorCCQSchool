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
# Load the UnicodePlots backend for plotting to the terminal
Plots.unicodeplots()

"""
    animate(f; nframes::Int, fps::Int = 30)

Call `f(i)` for each frame `i in 1:nframes`, where `f(i)` returns an object to print to the
terminal at each frame. Renders each frame in-place in the terminal at `fps` frames per
second.
"""
function animate(f; nframes::Int, fps::Real = 30)
    io = IOBuffer()
    # Hide cursor
    print(stdout, "\x1b[?25l")
    return try
        # Clear screen
        print(stdout, "\x1b[2J")
        for i in 1:nframes
            seekstart(io)
            show(io, MIME("text/plain"), f(i))
            frame_str = String(take!(io))
            # Move home
            print(stdout, "\x1b[H")
            println(frame_str)
            sleep(1 / fps)
        end
    finally
        # Show cursor again
        print(stdout, "\x1b[?25h")
    end
end

function main(;
        # Number of sites
        nsite = 30,
        # TEBD parameters
        time = 5.0,
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
        return exp(-im * timestep / 2 * hj)
    end
    # Include gates in reverse order too
    # (N, N - 1), (N - 1, N - 2), ...
    append!(gates, reverse(gates))

    # Make starting state
    j = nsite ÷ 2
    psit = apply(op("S+", sites[j]), psi)

    szs = [expect(psit, "Sz")]
    energies = ComplexF64[inner(psit', H, psi)]
    for current_time in 0.0:timestep:time
        psit = apply(gates, psit; cutoff)
        psit = normalize(psit)
        energy = inner(psit', H, psit)
        sz = expect(psit, "Sz")
        push!(szs, sz)
        push!(energies, energy)
        if floor(current_time - timestep + eps()) ≠ floor(current_time)
            if outputlevel > 0
                println("time: ", current_time)
                println("Bond dimension: ", maxlinkdim(psit))
                println("⟨ψₜ|Szⱼ|ψₜ⟩: ", sz[j])
                println("∑ⱼ⟨ψₜ|Szⱼ|ψₜ⟩: ", sum(sz))
                println("⟨ψₜ|H|ψₜ⟩: ", energy)
            end
        end
    end

    return (; energy, H, psi, szs, nsite, time, timestep, cutoff)
end
