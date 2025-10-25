# Day 2 Hands-On Tutorials

## Table of Contents

- [Tutorial 1: Real Time Evolution](#tutorial-1)
- [Tutorial 2: Imaginary Time Evolution](#tutorial-2)
- [Tutorial 3: Finite Temperature](#tutorial-3)

<a id="tutorial-1"></a>
<details>
  <summary><h2>Tutorial 1: Real Time Evolution</h2></summary>
  <hr>

In this tutorial we will simulate the time evolution of several initial states under the 1D spin-1/2 Heisenberg
Hamiltonian using the time evolving block decimation (TEBD) algorithm. See
the [ITensorMPS.jl tutorial on TEBD](https://docs.itensor.org/ITensorMPS/stable/tutorials/MPSTimeEvolution.html)
for more background on the algorithm.

The initial state constructed in `main` is the ground state of the Hamiltonian with the central spin excited. If we run this, the dynamics up until time $t=5.0$ will be simulated
```julia
julia> include("1-tebd.jl")
main

julia> res = main();
After sweep 1 energy=-13.09409636493637  maxlinkdim=10 maxerr=2.21E-03 time=0.033
After sweep 2 energy=-13.111284122964177  maxlinkdim=20 maxerr=2.27E-07 time=0.042
After sweep 3 energy=-13.111355718189827  maxlinkdim=46 maxerr=9.95E-11 time=0.083
After sweep 4 energy=-13.111355752014155  maxlinkdim=47 maxerr=9.93E-11 time=0.116
After sweep 5 energy=-13.111355752019346  maxlinkdim=47 maxerr=9.40E-11 time=0.111
time: 1.0
Bond dimension: 40
⟨ψₜ|Szⱼ|ψₜ⟩: 0.35502158812014994
∑ⱼ⟨ψₜ|Szⱼ|ψₜ⟩: 1.0000000000005
⟨ψₜ|H|ψₜ⟩: -11.929346376724634 - 2.0921203223608671e-16im

time: 2.0
Bond dimension: 47
⟨ψₜ|Szⱼ|ψₜ⟩: 0.08355269964032301
∑ⱼ⟨ψₜ|Szⱼ|ψₜ⟩: 1.0000000000007212
⟨ψₜ|H|ψₜ⟩: -11.929346333250871 - 4.398013998619169e-15im

time: 3.0
Bond dimension: 61
⟨ψₜ|Szⱼ|ψₜ⟩: -0.05027843626159247
∑ⱼ⟨ψₜ|Szⱼ|ψₜ⟩: 1.000000000000663
⟨ψₜ|H|ψₜ⟩: -11.929346353874365 - 1.3576397525559596e-15im

time: 4.0
Bond dimension: 81
⟨ψₜ|Szⱼ|ψₜ⟩: 0.00558010009598968
∑ⱼ⟨ψₜ|Szⱼ|ψₜ⟩: 1.0000000000004898
⟨ψₜ|H|ψₜ⟩: -11.929346384611089 + 4.990890715889457e-15im

time: 5.0
Bond dimension: 99
⟨ψₜ|Szⱼ|ψₜ⟩: 0.06961334092046022
∑ⱼ⟨ψₜ|Szⱼ|ψₜ⟩: 1.000000000002958
⟨ψₜ|H|ψₜ⟩: -11.929346264124122 - 3.920014918768254e-15im


julia> Plots.unicodeplots(); # Plot in the terminal

julia> plot_tebd_sz(res; step = 1) # S⁺|ψ⟩
     ┌────────────────────────────────────────┐
  0.5│⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⢷⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⠃⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡜⠀⠀⠀⡇⠀⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠀⠀⡀⠀⢀⡆⠀⡇⠀⠀⠀⢣⠀⡿⡀⢀⢆⠀⢠⡀⠀⢀⠀⠀⡀⠀⠀⠀⠀⠀│
⟨Szⱼ⟩│⡤⠦⢤⠴⠵⢤⠮⢦⢤⠮⢦⡼⠼⡤⡼⢼⢤⠧⠤⠤⠤⢼⢼⠤⢧⡼⠼⣤⠧⠵⣤⠧⠧⡴⠭⢦⠴⠵⠤⠴│
     │⠀⠀⠀⠀⠀⠀⠀⠀⠁⠀⠀⠁⠀⠱⠁⠈⣾⠀⠀⠀⠀⢸⡎⠀⠸⠁⠀⠈⠀⠀⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀│
     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⠀⠀⠀⠀⠀⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
 -0.5│⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
     └────────────────────────────────────────┘
     ⠀1⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀Site j⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀30⠀

julia> res.energies # Energy is approximately conserved
51-element Vector{ComplexF64}:
 -11.929346423742393 + 0.0im
  -11.92934641725427 - 3.6360627981309294e-15im
 -11.929346413315944 - 6.60465373202584e-15im
 -11.929346409625618 + 2.158849497713054e-15im
 -11.929346405215828 + 6.2303492855615294e-15im
                     ⋮
 -11.929346350879868 + 2.796818504340759e-15im
 -11.929346336090644 - 1.6361784674014638e-16im
  -11.92934631672588 - 2.352668301281432e-15im
 -11.929346291758428 - 4.641265263485175e-15im
 -11.929346259946318 - 8.775550971328244e-16im

julia> sum.(res.szs) # Total spin at each time is approximately conserved
51-element Vector{Float64}:
 0.9999999999864227
 1.0000000000026976
 1.0000000000026796
 1.0000000000026528
 ⋮
 1.000000000003992
 1.000000000005273
 1.00000000000655
 1.0000000000072606

julia> animate_tebd_sz(res) # Animation of ⟨Szⱼ⟩ as a function of time
[...]

```
The animation lets us visualize how the excitation propagates through the system as a function of time.

1. Included in `main()` is a function `vonneumann_ee(ψ::MPS, bond::Int = round(Int, length(ψ) / 2))` to compute the von Neumann entanglement entropy between sites `1..bond` and `bond+1...N` of the MPS. The vector of half-chain entanglement entropies is outputted by `main` as `vn_ees`.
Plot this half chain entanglement entropy as a function of time, how does it behave?

```julia
julia> Plots.unicodeplots()

julia> plot(res.times, res.vn_ees, xlabel = "Time", ylabel = "Entanglement")
            ┌────────────────────────────────────────┐  
     1.28615│⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠤⠔⠒⠲⠤⣀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│y1
            │⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡔⠁⠀⠀⠀⠀⠀⠀⠉⠒⠢⣄⣀⠀⠀⠀⢀⣀⠤⠀│  
            │⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠎⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠉⠉⠁⠀⠀⠀│  
            │⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
            │⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
            │⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠎⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
            │⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠎⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
Entanglement│⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡜⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
            │⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
            │⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡰⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
            │⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⡔⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
            │⠀⡇⠀⠀⠀⠀⠀⠀⠀⢀⠴⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
            │⠀⡇⠀⠀⠀⠀⠀⠀⢠⠊⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
            │⠀⡇⠀⠀⠀⣀⡠⠊⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
     0.43001│⠀⡧⠔⠒⠋⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
            └────────────────────────────────────────┘  
            ⠀-0.18⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀Time⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀6.18⠀ 
```
Is this what you would expect for a quench? Why or why not? What happens around time `t ~ 5.0`? Try increasing the time of the simulation to `time = 8.0` to resolve the long-time behaviour better.

2. We can change the initial state to something different. Let's try a state where all the spins are polarised along the z-axis. This can be done via the line
```julia
julia> psit = ITensorMPS.MPS(sites, ["Z+" for i in 1:nsite])
```
Note that you will have to comment out parts of the code where the initial state was created by DMRG and then excited. It is sufficient to comment out lines 81-85:
```julia
    psi0 = random_mps(sites; linkdims = 10)
    initial_energy, psi = dmrg(
        H, psi0; nsweeps = 5, maxdim = [10, 20, 100, 100, 200],
        cutoff = [1.0e-10], outputlevel = min(outputlevel, 1)
    )
```

and lines 100-102:
```julia
    j = nsite ÷ 2
    psit = apply(op("S+", sites[j]), psi)
    psit = normalize(psit)
```

and replace them with 

```julia
    psit = ITensorMPS.MPS(sites, ["Z+" for i in 1:nsite])
    initial_energy = inner(psit', H, psit)
```

What do you notice about the dynamics of the quench now? Hint: think about the symmetries of the model.

3. Now try initializing the system in an anti-ferromagnetic state instead
```julia
    psit = ITensorMPS.MPS(sites, [iseven(i) ? "Z+" : "Z-" for i in 1:nsite])
```

Plot the entanglement entropy as a function of time.

```julia
julia> plot(res.times, res.vn_ees, xlabel = "Time", ylabel = "Entanglement")
            ┌────────────────────────────────────────┐  
      3.0758│⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⠔⠀│y1
            │⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⠴⠊⠀⠀⠀│  
            │⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠤⠊⠀⠀⠀⠀⠀⠀│  
            │⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠤⠚⠁⠀⠀⠀⠀⠀⠀⠀⠀│  
            │⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡠⠒⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
            │⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣠⠖⠊⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
            │⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠔⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
Entanglement│⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⡠⠊⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
            │⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⠔⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
            │⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⡠⠊⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
            │⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⣀⠔⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
            │⠀⡇⠀⠀⠀⠀⠀⢀⡠⠊⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
            │⠀⡇⠀⠀⠀⢀⠔⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
            │⠀⡇⠀⢀⠜⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
  -0.0895865│⠤⡷⠮⠥⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤│  
            └────────────────────────────────────────┘  
            ⠀-0.18⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀Time⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀6.18⠀ 
``` 

How does this differ to the first initial state we used (the locally excited ground state)? What does this imply for the growth of the bond dimension of the function of time to maintain accuracy? Hint: for an arbitrary MPS of bond dimension $\chi$, the bound $S_{\rm Von Neumann} \leq log_{2}(\chi)$ is true. 

Click [here](#table-of-contents) to return to the table of contents.

</details>

<a id="tutorial-2"></a>
<details>
  <summary><h2>Tutorial 2: Imaginary Time Evolution</h2></summary>
  <hr>

```julia
julia> res = main();
After sweep 1 energy=-13.10580711255933  maxlinkdim=10 maxerr=2.04E-03 time=0.029
After sweep 2 energy=-13.111348929097458  maxlinkdim=20 maxerr=1.41E-07 time=0.040
After sweep 3 energy=-13.11135575001343  maxlinkdim=45 maxerr=9.81E-11 time=0.085
After sweep 4 energy=-13.111355751942149  maxlinkdim=47 maxerr=1.00E-10 time=0.118
After sweep 5 energy=-13.111355751949796  maxlinkdim=47 maxerr=1.00E-10 time=0.112
time: 5.0
Bond dimension: 24
⟨ψₜ|Szⱼ|ψₜ⟩: -0.07015198148930198
∑ⱼ⟨ψₜ|Szⱼ|ψₜ⟩: -0.3554642454935465
⟨ψₜ|H|ψₜ⟩: -12.918726195417213

time: 10.0
Bond dimension: 38
⟨ψₜ|Szⱼ|ψₜ⟩: -0.0007049850288560583
∑ⱼ⟨ψₜ|Szⱼ|ψₜ⟩: -0.12498385296119857
⟨ψₜ|H|ψₜ⟩: -13.082551163963094

time: 15.0
Bond dimension: 40
⟨ψₜ|Szⱼ|ψₜ⟩: 0.007906094151056576
∑ⱼ⟨ψₜ|Szⱼ|ψₜ⟩: -0.035190911135993715
⟨ψₜ|H|ψₜ⟩: -13.105243727446057

time: 20.0
Bond dimension: 40
⟨ψₜ|Szⱼ|ψₜ⟩: 0.005455558345839978
∑ⱼ⟨ψₜ|Szⱼ|ψₜ⟩: -0.010092112111528002
⟨ψₜ|H|ψₜ⟩: -13.109727125683662


julia> res.energies .- res.energy
101-element Vector{Float64}:
 13.344740259203546
 11.607416095173846
  9.788767514056136
  8.027055433275628
  6.478091389576792
  5.219661141730995
  ⋮
  0.00210452773845482
  0.0019989336674051117
  0.001898845407890093
  0.0018039541390795222
  0.0017139718695293737
  0.0016286262661342477

```

Click [here](#table-of-contents) to return to the table of contents.

</details>

<a id="tutorial-3"></a>
<details>
  <summary><h2>Tutorial 3: Finite Temperature</h2></summary>
  <hr>

Click [here](#table-of-contents) to return to the table of contents.

</details>
