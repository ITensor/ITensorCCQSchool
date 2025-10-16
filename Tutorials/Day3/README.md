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
$Z = \sum_{s_{1} \in {-1, 1}}\sum_{s_{2} \in {-1, 1}} \hdots \sum_{s_{L} \in {-1, 1}}\exp(-beta \sum_{ij}s_{i}.s_{j})$

You can inspect the individual tensors on each vertex of the constructed tensor network via `res.tensornetwork[v]` where `v` is the name of the vertex.
```
julia> res = main(L=3, periodic = false, β = 0.2);

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

julia> res = main(n=3, periodic = false);
julia> @show res.z;
res.z = 4.16214474367691

In 1D the partition function of the Ising model is analytically computable for any system size L and both Periodic and Open Boundaries. The results are
  $Z_{L,OBC} = 2*(2*\cosh(\beta)^{L-1})$
for open boundaries and
  $Z_{L, PBC} = 2*\cosh(\beta)^{L} + 2*\sinh(\beta)^{L}$
for periodic boundaries.

2. Compare the output of `res.z` with these values for both periodic and open boundaries. Do they agree?