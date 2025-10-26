# Day 2 Hands-On Tutorials

## Table of Contents

- [Tutorial 1: Real Time Evolution](#tutorial-1)
- [Tutorial 2: Imaginary Time Evolution](#tutorial-2)
- [Tutorial 3: Finite Temperature](#tutorial-3)
- [Stretch Goals](#stretch-goals)

<a id="tutorial-1"></a>
<details>
  <summary><h2>Tutorial 1: Real Time Evolution</h2></summary>
  <hr>

In this tutorial we will simulate the time evolution of several initial states under the 1D spin-1/2 Heisenberg
Hamiltonian using the time evolving block decimation (TEBD) algorithm. See
the [ITensorMPS.jl tutorial on TEBD](https://docs.itensor.org/ITensorMPS/stable/tutorials/MPSTimeEvolution.html)
for more background on the algorithm. We will work off of the script [1-tebd.jl](./1-tebd.jl).

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

1. Included in `main()` is a function `entanglement_entropy(ψ::MPS, bond::Int = length(ψ) ÷ 2)` to compute the von Neumann entanglement entropy between sites `1..bond` and `bond+1...N` of the MPS. The vector of half-chain entanglement entropies is output by `main` as `entanglements`.
Plot this half chain entanglement entropy as a function of time, how does it behave?

```julia
julia> Plots.unicodeplots()

julia> plot(res.times, res.entanglements, xlabel = "Time", ylabel = "Entanglement")
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
julia> psit = MPS(sites, ["Z+" for i in 1:nsite])
```
Note that you should comment out parts of the code where the initial state was created by DMRG and then excited (lines 81-85 and 100-102) and substitute them for:

```julia
    psit = ITensorMPS.MPS(sites, ["Z+" for i in 1:nsite])
    initial_energy = inner(psit', H, psit)
```

alternatively, if you are feeling lazy, you can just paste the above on line 103 to overwrite the current inital state.

What do you notice about the dynamics of the quench now? Hint: think about the symmetries of the model.

3. Now try initializing the system in an anti-ferromagnetic state instead
```julia
    psit = ITensorMPS.MPS(sites, [iseven(i) ? "Z+" : "Z-" for i in 1:nsite])
```

Plot the entanglement entropy as a function of time.

```julia
julia> plot(res.times, res.entanglements, xlabel = "Time", ylabel = "Entanglement")
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

Now we are going to switch from real time to imaginary time evolution. This is incredibly easy with tensor networks, as we can just perform the substitution $dt \rightarrow - {\rm i} d \beta`. 

We will be working off the script [2-imaginary-time.jl](./2-imaginary-time.jl) which does this for you and implements the imaginary time dynamics of a random initial state under the Heisenberg Hamiltonian.


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

1. Notice how the energy is converging to that of the DMRG calculation. You can show an animation of the local $Sz$ of each spin in the chain by passing `outputlevel = 2` as a `kwarg` to `main()`. Observe how the system relaxes to a state with no local magnetisation. 

We can calculate the variance of `psit` to observe how close it is to an eigenstate of `H`. Specifically the variance is given by

```julia
julia> energy_var = inner(apply(H, psit),apply(H, psit)) - inner(psit',H,psit)^2
```

2. Calculate the variance of the energy as a function of time in your simulation and have `main` return it. Plot it. 

```julia
julia> plot(res.times, res.energy_vars, xlabel = "Imaginary Time", ylabel = "Energy Variance")
               ┌────────────────────────────────────────┐  
        4.67397│⠀⡷⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│y1
               │⠀⡇⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
               │⠀⡇⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
               │⠀⡇⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
               │⠀⡇⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
               │⠀⡇⠸⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
               │⠀⡇⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
Energy Variance│⠀⡇⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
               │⠀⡇⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
               │⠀⡇⠀⠘⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
               │⠀⡇⠀⠀⢧⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
               │⠀⡇⠀⠀⠘⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
               │⠀⡇⠀⠀⠀⢣⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
               │⠀⡇⠀⠀⠀⠀⠣⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
      -0.135919│⠤⡧⠤⠤⠤⠤⠤⠬⠶⠶⠶⠶⠦⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤│  
               └────────────────────────────────────────┘  
               ⠀-0.6⠀⠀⠀⠀⠀⠀⠀⠀Imaginary Time⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀20.6⠀  
```

The initial state we used is a random `MPS`constructed via the lines
```julia>
    rng = StableRNG(123)
    psit = random_mps(rng, sites)
```

3. Try changing the seed of the initial state. Does the result still converge to the ground state? Can you think of what initial states might prevent this happening? Hint: think about the symmetries of the model. Try to construct some. Does the variance still go to zero?


Click [here](#table-of-contents) to return to the table of contents.

</details>

<a id="tutorial-3"></a>
<details>
  <summary><h2>Tutorial 3: Finite Temperature</h2></summary>
  <hr>

We are now going to run the METTS (Minimally entangled thermal states) algorithm to extract finite temperature properties of the system whilst remaining in the pure state picture. This is done in the file [3-metts.jl](./3-metts.jl).

```julia
julia> include("3-metts.jl")

julia> res = main();
Making warmup METTS number 10
  Sampled state: ["Z-", "Z+", "Z+", "Z-", "Z+", "Z-", "Z+", "Z-", "Z+", "Z-"]
Making METTS number 10
  Energy of METTS 10 = -4.0978
  Energy of ground state from DMRG -4.2580
  Estimated Energy = -4.1171 +- 0.0318  [-4.1489,-4.0853]
  Sampled state: ["Z-", "Z-", "Z+", "Z-", "Z-", "Z+", "Z+", "Z-", "Z-", "Z+"]
Making METTS number 20
  Energy of METTS 20 = -3.5787
  Energy of ground state from DMRG -4.2580
  Estimated Energy = -4.0288 +- 0.0453  [-4.0741,-3.9836]
  Sampled state: ["Z+", "Z-", "Z-", "Z+", "Z-", "Z+", "Z+", "Z-", "Z-", "Z+"]
Making METTS number 30
  Energy of METTS 30 = -4.1557
  Energy of ground state from DMRG -4.2580
  Estimated Energy = -4.0074 +- 0.0370  [-4.0443,-3.9704]
  Sampled state: ["Z-", "Z+", "Z+", "Z+", "Z-", "Z-", "Z+", "Z-", "Z+", "Z-"]
Making METTS number 40
  Energy of METTS 40 = -3.8662
  Energy of ground state from DMRG -4.2580
  Estimated Energy = -3.9767 +- 0.0333  [-4.0100,-3.9434]
  Sampled state: ["Z-", "Z+", "Z-", "Z-", "Z-", "Z+", "Z+", "Z-", "Z+", "Z-"]
Making METTS number 50
  Energy of METTS 50 = -3.8099
  Energy of ground state from DMRG -4.2580
  Estimated Energy = -3.9556 +- 0.0301  [-3.9856,-3.9255]
  Sampled state: ["Z+", "Z+", "Z-", "Z+", "Z-", "Z-", "Z+", "Z-", "Z-", "Z+"]
Making METTS number 60
  Energy of METTS 60 = -3.8595
  Energy of ground state from DMRG -4.2580
  Estimated Energy = -3.9723 +- 0.0264  [-3.9987,-3.9459]
  Sampled state: ["Z+", "Z-", "Z+", "Z-", "Z+", "Z+", "Z-", "Z+", "Z-", "Z+"]
Making METTS number 70
  Energy of METTS 70 = -4.1383
  Energy of ground state from DMRG -4.2580
  Estimated Energy = -3.9683 +- 0.0236  [-3.9919,-3.9447]
  Sampled state: ["Z-", "Z+", "Z-", "Z+", "Z-", "Z+", "Z-", "Z-", "Z+", "Z-"]
Making METTS number 80
  Energy of METTS 80 = -4.1383
  Energy of ground state from DMRG -4.2580
  Estimated Energy = -3.9765 +- 0.0217  [-3.9982,-3.9548]
  Sampled state: ["Z-", "Z+", "Z-", "Z+", "Z-", "Z+", "Z+", "Z-", "Z+", "Z-"]
Making METTS number 90
  Energy of METTS 90 = -4.1383
  Energy of ground state from DMRG -4.2580
  Estimated Energy = -3.9760 +- 0.0203  [-3.9963,-3.9557]
  Sampled state: ["Z-", "Z+", "Z-", "Z-", "Z-", "Z+", "Z+", "Z-", "Z+", "Z-"]
Making METTS number 100
  Energy of METTS 100 = -3.8843
  Energy of ground state from DMRG -4.2580
  Estimated Energy = -3.9652 +- 0.0190  [-3.9842,-3.9462]
  Sampled state: ["Z-", "Z-", "Z+", "Z-", "Z+", "Z+", "Z-", "Z+", "Z+", "Z-"]
```

The specific heat can be approximated from the METTS algorithm via the following formula

$$C_{v}(\beta) = \frac{\beta^{2}}{\rm NMETTS}\left(\overline{\langle H^{2} \rangle} - \overline{\langle H \rangle^{2}} \right)

where $\overline$ denotes the METTS ensemble average. You can measure the square energy of a given METTS via 

```julia
julia> inner(apply(H, psi),apply(H, psi))
```

1. Modify `main()` to keep track of the square energy of each METT after it has been evolved. Average over these, and the energies (which are already kept track off) at the end of the simulation to calculate $C_{v}(\beta)$ for the given $\beta$ and have it returned by main. Check that this gives a sensible answer from `main()`. For the default parameters ($\beta = 4.0$, NMETTS $=100$) provided you should find $C_{v}(\beta = 4.0) \approx 0.19$ (the RNG for the initial state and sampling is seeded to be reproducable).

```julia
julia> res = main(; outputlevel = 0);

julia> res.specific_heat
0.1906972673732355
```

No we are going to measure the specific heat as a function of inverse temperature.

2. Construct an array of $\beta$ values spanning $0 \leq \beta \leq 8.0$, for instance

```julia

julia> betas = [0.2*i for i in 1:41];
```

and then create a vector of simulation outputs for these `betas`. E.g

```julia
julia> results = [main(; beta, betastep = 0.1, NMETTS=25) for beta in betas]
```

This might take a few minutes to run, so play around with the setting of `NMETTS`. We suggest setting NMETTS = 25 to get a coarse grained result. Plot the result.

```julia
julia> plot(betas, specific_heats, xlabel = "Beta", ylabel = "Specific Heat")
             ┌────────────────────────────────────────┐  
     0.393288│⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢰⠀⠀⠀⠀⠀⣀⠀⠀⠀⠀⡦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│y1
             │⡇⠀⠀⠀⠀⠀⠀⢰⡄⠀⠀⠀⢀⢿⠀⠀⠀⠀⢰⢹⠀⠀⠀⠀⡇⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
             │⡇⠀⠀⠀⠀⠀⠀⡜⠈⠢⠔⠢⠊⠀⡇⠀⠀⢰⠁⢸⠀⠀⠀⢀⠇⢣⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
             │⡇⠀⠀⠀⠀⠀⢀⠇⠀⠀⠀⠀⠀⠀⡇⠀⠀⡇⠀⢸⠀⠀⠀⢸⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
             │⡇⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⢇⠀⠀⡇⠀⢸⠀⠀⠀⢸⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
             │⡇⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⢸⡀⢸⠀⠀⠀⡇⠀⠀⡸⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
             │⡇⠀⠀⠀⠀⢠⠃⠀⠀⠀⠀⠀⠀⠀⠀⠈⠺⠀⠀⠀⡇⠀⠀⡇⠀⠸⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
Specific Heat│⡇⠀⠀⠀⠀⡸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⡦⡸⠀⠀⠀⡇⣦⠀⠀⠀⠀⠀⠀⠀⢀⠀⠀⠀⠀│  
             │⡇⠀⠀⠀⢀⠇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡿⠀⠀⠀⠀⠀⣿⠀⡇⠀⡠⢄⡀⠀⠀⡜⡇⠀⠀⠀│  
             │⡇⠀⠀⠀⡜⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠃⠀⢱⢰⠁⠀⢣⠀⢰⠁⢱⠀⠀⠀│  
             │⡇⠀⠀⢰⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⡸⠀⠀⠈⠒⠃⠀⠘⠲⡀⠀│  
             │⡇⠀⠀⡜⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠇⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀│  
             │⡇⠀⢠⠇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠸⠀│  
             │⡇⢀⠎⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
  -0.00416347│⣇⣎⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀│  
             └────────────────────────────────────────┘  
             ⠀-0.034⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀Beta⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀8.234⠀ 
```

The specific heat of the spin 1/2 antiferromagnetic Heisenberg model is known to display a broad peak at $T = 0.48J$ (here we have $J=  1$) with a maximum value of $~0.35J$. Do your results agree with this? 

3. The high temperature regime should display an inverse square dependence of the specific heat with temperature, i.e $C_{v} \propto \frac{1}{T^{2}}$. Use a range $0 \leq \beta \leq 0.5$ to try to confirm this. When using a finer range of betas, make sure to adjust the step size in `main` to be commensurate or you will get an error message.

```julia

julia> plot([beta*beta for beta in betas], specific_heats, xlabel = "Beta", ylabel = "Specific Heat")
             ┌────────────────────────────────────────┐  
     0.047881│⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠤⠀│y1
             │⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⡠⠒⠁⠀⠀│  
             │⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡠⠔⠁⠀⠀⠀⠀⠀│  
             │⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡠⠔⠊⠀⠀⠀⠀⠀⠀⠀⠀│  
             │⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡠⠔⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
             │⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡠⠒⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
             │⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠔⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
Specific Heat│⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠎⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
             │⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠔⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
             │⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠔⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
             │⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡠⠔⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
             │⢸⠀⠀⠀⠀⠀⠀⠀⠀⡠⠔⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
             │⢸⠀⠀⠀⠀⠀⣠⠔⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
             │⢸⠀⠀⡠⠔⠊⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
   -0.0010349│⢼⠴⠯⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤│  
             └────────────────────────────────────────┘  
             ⠀-0.004925⠀⠀⠀⠀⠀⠀⠀⠀Beta⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀0.257425⠀
```

</details>

<a id="stretch-goals"></a>
<details>
  <summary><h2>Stretch Goals</h2></summary>
  <hr>

If you completed all the tutorials and would like more of a challenge, choose from among the   following "stretch goal" activities.

1. In the low temperature regime the spin 1/2 1D Heisenberg model is known to be a gapless Luttinger Liquid which is a phase of matter characterised by a specific heat $C_{v} \propto T$. See if you can confirm this by running the METTS code in the low temperature regime (say $8.0 \leq \beta \leq 10.0$) and measuring the specific heat capacity. Note that in this low-temperature regime, fluctuations and finite-size effects will be more significant (we have been working on a small chain), so you will have to be careful about the parameters you choose and simulations could take some time. It can help to take a large enough `betastep` (say `betastep = O(0.1)`) so your simulations run in reasonable time.

Click [here](#table-of-contents) to return to the table of contents.

</details>
