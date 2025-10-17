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

julia> animate(; nframes = length(res.szs), fps = res.nsite) do i
           return plot(
               res.szs[i]; xlim = (1, res.nsite), ylim = (-0.5, 0.5), xlabel = "Site j",
               ylabel = "⟨Szⱼ⟩", legend = false, title = "Sweep = $(i ÷ (2 * res.nsite) + 1)"
           )
       end
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
