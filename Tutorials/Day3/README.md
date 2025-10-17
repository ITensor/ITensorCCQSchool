# Day 3 Hands-On Tutorials

## Table of Contents

- [Tutorial 1](#tutorial-1)
- [Tutorial 2](#tutorial-2)
- [Tutorial 3](#tutorial-3)

</details>

<a id="tutorial-1"></a>
<details>
  <summary><h2>Tutorial 1</h2></summary>
  <hr>

We are going to combine the `NamedGraphs.jl` and `ITensors.jl` packages to build tensor networks of varying topology. 

A simple graph `g` is just a series of vertices and edges between pairs of those vertices. There are no multiedges or self edges. The package `NamedGraphs.jl` is built around the `NamedGraph` object `g`, which can be constructed using either the pre-built graph constructors or our own via code like 

```
  julia> using NamedGraphs: NamedGraph, NamedEdge
  julia> g = NamedGraph([1,2,3])
  julia> edges = [NamedEdge(1, 2), NamedEdge(2,3)]
  julia> g = add_edges(g, edges)
```

First, lets run the  script [1-tensornetworks.jl](./1-tensornetworks.jl)

```
julia> include("1-tensornetworks.jl")
main (generic function with 1 method)
```

You will see that it builds the 3-site path graph, which can be accessed and viewed via

```
julia> res = main();
julia> @show res.g
g = NamedGraph{Int64} with 3 vertices:
3-element NamedGraphs.OrderedDictionaries.OrderedIndices{Int64}:
 1
 2
 3

and 2 edge(s):
1 => 2
2 => 3
```

1: Modify the graph construction in `main()` to create a path graph on `L` vertices, where `L` is an integer variable that can be specified as a keyword argument to main. Compare the output to the pre-written constructor `named_path_graph(L::Int)` in `NamedGraphs.jl`. Add in a `periodic` flag to your constructor to add a periodic boundary if the flag is true.
 
With this you should be able to do
```
julia> res = main(; L = 5, periodic = true);
julia> @show res.g
g = NamedGraph{Int64} with 5 vertices:
5-element NamedGraphs.OrderedDictionaries.OrderedIndices{Int64}:
 1
 2
 3
 4
 5

and 5 edge(s):
1 => 2
1 => 5
2 => 3
3 => 4
4 => 5
```

We can build a tensor network as a dictionary of tensors, one for each vertex of the `NamedGraph` `g`. The edges of the graph `g` (which are of the  type `NamedEdge`) dictate which tensors share indices to be contracted over. 

Provided in [1-tensornetworks.jl](./1-tensornetworks.jl) is a pre-built constructor for the tensor network representing the partition function of the ising model on a given `NamedGraph` g at a given inverse temperature `β`. The partition function reads 
$$Z(\beta) = \frac{1}{2}\sum_{s_{1} \in {-1, 1}}\sum_{s_{2} \in {-1, 1}} ... \sum_{s_{L}\in {-1, 1}}\exp(-\beta \sum_{ij}s_{i}.s_{j}),$$
where we have scaled by a factor of 1/2 for convenience.

You can inspect the individual tensors on each vertex of the constructed tensor network via `res.tensornetwork[v]` where `v` is the name of the vertex.
```
julia> res = main(L=3, periodic = false, beta = 0.2);

julia> @show res.tensornetwork[1]
res.tensornetwork[1] = ITensor ord=1
Dim 1: (dim=2|id=103|"e1_2")
NDTensors.Dense{Float64, Vector{Float64}}
 2-element
 1.2724995249301703
 1.2724995249301705
ITensor ord=1 (dim=2|id=103|"e1_2")
NDTensors.Dense{Float64, Vector{Float64}}
```

This tensornetwork can be contracted with the `contract_tensornetwork` function provided. Its output is pre-computed for you in `main()`

```
julia> res = main(n=3, periodic = false);
julia> @show res.z;
res.z = 2.081072371838455
```

In 1D the partition function of the Ising model is analytically computable for any system size L and both Periodic and Open Boundaries. The results are
$$Z_{L,OBC}(\beta) = 2\cosh^{L-1}(\beta)$$
for open boundaries and
$$Z_{L, PBC}(\beta) = \cosh^{L}(\beta) + \sinh^{L}(\beta)$$
for periodic boundaries.

2. Compare the output of `res.z` with these values for both periodic and open boundaries. Do they agree? If they do, then congratulations, you just solved the 1D PBC and OBC Ising model with a tensor network approach.

<a id="tutorial-2"></a>
<details>
  <summary><h2>Tutorial 2</h2></summary>
  <hr>

In the previous tutorial, the `contract_tensornetwork()` function contracted the tensor network exactly by multiplying the tensors together, vertex by vertex. This can only be done efficiently for tree-like networks (those composed of no loops, or a small number of loops) and only when taking careful care over the order of contraction.

In this tutorial we are going to use the script [2-beliefpropagation.jl](./2-beliefpropagation.jl) belief propagation to contract tensor networks in an efficient, but approximate manner.

The script now builds an $nx \times ny$ square grid tensornetwork representing the partition function of the Ising model in 2D. Inverse temperature is set via the `beta` kwarg and periodic boundaries (in both directions) can be added with the kwarg `periodic`. Returned is the number of iterations BP took to converge, and the rescaled free energy density 
$$\phi(\beta) = -\beta f(\beta) = log(Z(\beta))/(nx*ny)$$

We can do the following to get the BP computed value for $\phi$ on a 10x1 OBC square grid. This is just a path graph, like in the previous example.
```
julia> include("2-beliefpropagation.jl")
main (generic function with 1 method)

julia> res = main(; Lx = 3, Ly = 1, beta=  0.2, periodic = false);
BP Algorithm Converged after 3 iterations

julia> res.bp_phi
0.24429444141332002
```
1. Compare the result to the analytical value for 1D OBC
$$ \phi_{Lx,OBC}(\beta) = \frac{1}{Lx}ln(2\cosh^{Lx-1}(\beta))$$

They agree, even though we used BP to compute it. Why?

2. We can also get the bp approximated free energy density for a periodic ring. 
```
julia> res = main(; Lx=  3, Ly = 1, periodic = true);
BP Algorithm Converged after 8 iterations

julia> res.bp_phi
0.019868071835749606
```

2. Compare the result to the 1D scaled free energy density on PBC, 
$$\phi_{Lx,OBC}(\beta) = \frac{1}{Lx}ln(\cosh^{Lx}(\beta) + \sinh^{Lx}(\beta))$$

They don't agree. Why? Pick a finite value of $\beta$ between $0$ and $1$ and compute both the exact PBC free energy vs $Lx$ for $Lx = 3,4,...30$ and the `bp` free energy using the `main` function (set $Ly = 1$ and `periodic = true`).

Plot the error between the bp approximated free energy density and
the exact free energy density as a function of $L_{x}$ on a log scale. What's the scaling? Why?

```
julia> plot([Lx for Lx in 3:25], bp_abs_errors, yscale = :ln)
         ┌────────────────────────────────────────┐  
ℯ⁻³⸱³⁸⁴⁹⁹│⠀⠢⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│y1
         │⠀⠀⠈⠢⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
         │⠀⠀⠀⠀⠈⠑⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
         │⠀⠀⠀⠀⠀⠀⠀⠱⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠢⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠢⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠒⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⠢⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠢⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠑⠢⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠢⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠢⡀⠀⠀⠀⠀⠀⠀⠀⠀│  
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠑⠢⣀⠀⠀⠀⠀⠀│  
ℯ⁻²⁶⸱¹²⁵⁴│⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠒⠒⠤⠀│  
         └────────────────────────────────────────┘  
         ⠀2.34⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀25.66⠀ 
```

Inspect the values for `phi` returned by `bp` versus system size? Do you notice something odd? Why are they all the same value?

Now we're going to move fully into 2D. Let's compute the BP approximate free energy density on a OBC square grid with $L_{x} = L$ and $L_{y} = L$ as a function of $\beta$.

```
julia> betas =[0.05*(i-1) for i in 1:21]
julia> bp_phis = [main(; Lx=15, Ly = 15, periodic = false, beta).bp_phi for beta in betas]
```

Congratulations. You just approximately solved the 2D Ising model on a 15x15 square lattice for twenty different inverse temperatures in about 10 seconds.

3. How does the number of iterations that BP took to converge depend on the inverse temperature? Plot this. Where's the peak? Is it near the critical point of the 2D model? Or somewhere different?

Included in `[2-beliefpropagation.jl](./2-beliefpropagation.jl)` is a function for computing the exact rescaled free energy of the 2D model in the thermodynamic limit via Onsager's famous result. This is returned by `main` as `exact_phi_onsager`.

$$\phi(\beta) = \beta f(\beta) =
= -\ln 2 +

* \frac{1}{8\pi^{2}}
  \int_{0}^{2\pi}!!\int_{0}^{2\pi}
  \ln!\left[
  \cosh!\left(2\beta J_{1}\right)\cosh!\left(2\beta J_{2}\right)
  -\sinh!\left(2\beta J_{1}\right)\cos!\left(\theta_{1}\right)
  -\sinh!\left(2\beta J_{2}\right)\cos!\left(\theta_{2}\right)
  \right],
  d\theta_{1}, d\theta_{2}.$$

Lets compare our results to that.

5. Pick a small value for $\beta$ (say $\beta = 0.1$) and plot the error between `bp` and the `exact` result as a function of lattice size $L$ for $L_{x} = L$ and $L_{y} = L$. How does it scale?

Now lets move to periodic boundary conditions. 
```
julia> res = main(; Lx = 5, Ly = 5, periodic = true, beta = 0.2)
BP Algorithm Converged after 21 iterations
(bp_phi = -0.6534110369600732, exact_phi_onsager = -0.6517635488435647, niterations = 21)
```
6. What do you notice about the dependence of `bp_phi` on $L$?


As BP is letting us work directly in the thermodynamic limit with periodic boundaries, pick a small $L = 3$ and a fine-range of betas

```
julia> betas = [0.01*(i-1) for i in 1:101]
```

Plot the absolute error between BP and Onsager's result. Where does it peak? 

```
julia> plot(betas, errs)
            ┌────────────────────────────────────────┐  
   0.0181699│⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡰⡆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│y1
            │⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⢱⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
            │⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⠃⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
            │⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
            │⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡜⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
            │⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⢱⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
            │⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
            │⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡜⠀⠀⠀⠸⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
            │⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠇⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
            │⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⡸⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
            │⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⢠⠃⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
            │⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⡜⠀⠀⠀⠀⠀⠀⠘⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
            │⠀⡇⠀⠀⠀⠀⠀⠀⠀⡜⠀⠀⠀⠀⠀⠀⠀⠀⢇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
            │⠀⡇⠀⠀⠀⠀⠀⢠⠔⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠣⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  
-0.000529221│⠤⡧⠤⠤⠴⠶⠮⠥⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠬⠽⠶⠶⠦⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤│  
            └────────────────────────────────────────┘  
            ⠀-0.03⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀1.03⠀  
```
<a id="tutorial-3"></a>
<details>
  <summary><h2>Tutorial 3</h2></summary>
  <hr>

Now we are going to try to correct our BP results with a first order cluster expansion.