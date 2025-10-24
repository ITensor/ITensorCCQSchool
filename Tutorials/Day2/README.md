# Day 2 Hands-On Tutorials

## Table of Contents

- [Tutorial 1: Real Time Evolution](#tutorial-1)
- [Tutorial 2: Imaginary Time Evolution](#tutorial-2)
- [Tutorial 3: Finite Temperature](#tutorial-3)

<a id="tutorial-1"></a>
<details>
  <summary><h2>Tutorial 1: Real Time Evolution</h2></summary>
  <hr>

In this tutorial we will run time evolution of a spin applied to the 1D spin-1/2 Heisenberg
model ground state using the time evolving block decimation (TEBD) algorithm. See
the [ITensorMPS.jl tutorial on TEBD](https://docs.itensor.org/ITensorMPS/stable/tutorials/MPSTimeEvolution.html)
for more background on the algorithm.
```julia
julia> include("1-tebd.jl")
main

julia> res = main();

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

Click [here](#table-of-contents) to return to the table of contents.

</details>

<a id="tutorial-2"></a>
<details>
  <summary><h2>Tutorial 2: Imaginary Time Evolution</h2></summary>
  <hr>

Click [here](#table-of-contents) to return to the table of contents.

</details>

<a id="tutorial-3"></a>
<details>
  <summary><h2>Tutorial 3: Finite Temperature</h2></summary>
  <hr>

Click [here](#table-of-contents) to return to the table of contents.

</details>
