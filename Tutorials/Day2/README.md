# Day 2 Hands-On Tutorials

## Table of Contents

- [Tutorial 1](#tutorial-1)
- [Tutorial 2](#tutorial-2)
- [Tutorial 3](#tutorial-3)

<a id="tutorial-1"></a>
<details>
  <summary><h2>Tutorial 1</h2></summary>
  <hr>

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
 -1.1158590813174532e-7 + 0.0im
    -11.929346421893666 + 6.609223542149849e-17im
    -11.929346417949144 + 6.4987220553615346e-15im
    -11.929346414256656 + 1.7244575134826174e-16im
                        ⋮
    -11.929346340685345 - 6.651226040999768e-15im
    -11.929346321308053 + 7.668303212388609e-16im
    -11.929346296331948 - 1.2574373847093401e-15im
    -11.929346264518774 - 7.415893341027481e-15im

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
  <summary><h2>Tutorial 2</h2></summary>
  <hr>

Click [here](#table-of-contents) to return to the table of contents.

</details>

<a id="tutorial-3"></a>
<details>
  <summary><h2>Tutorial 3</h2></summary>
  <hr>

Click [here](#table-of-contents) to return to the table of contents.

</details>
