# ITensor CCQ School Hands-On Tutorials Introduction

Welcome to Day 3

01-tensornetworks.jl

We are going to use the `NamedGraphs.jl` package to build tensor networks of varying topology. A simple graph `g` is just a series of vertices and edges between pairs of 
those vertices.

`NamedGraphs.jl` designated these as a `NamedGraph` object and they can be constructed using either the pre-built graph constructors or our own via code like 

  # Create a simple path graph of 3 vertices
  g = NamedGraph([1,2,3])
  edges = [NamedEdge(1, 2), NamedEdge(2,3)]
  g = add_edges(g, edges)

Task 1: Write a function to create a path graph on `n` vertices. Compare the output to the pre-written constructor `named_path_graph(n::Int)` in `NamedGraphs.jl`. Add in a flag to your function to add a periodic boundary. 

We can then build a tensor network as a dictionary of tensors, one for each vertex of the `NamedGraph`. The edges of the graph `g` (which are of the  type `NamedEdge`) dictate which tensors share indices to be contracted over.

Provided in `01-tensornetworks.jl` is a constructor for the tensor network representing the partition function of the classical ising model on an arbitrary graph. This is a bond dimension 2 network and its contraction equals the partition function
  $Z = sum_{s_{1} \in {-1, 1}}sum_{s_{2} \in {-1, 1}}...sum_{s_{L} \in {-1, 1}}exp(-beta \sum_{ij}s_{i}.s_{j})$
For the given graph.

In 1D this partition function is given by
  $Z_{N,OBC} = 2*(2*\cosh(\beta)^{L-1})$
for open boundaries and
  $Z_{N, PBC} = 2*\cosh(\beta)^{L} + 2*\sinh(\beta)^{L}$
for periodic boundaries.

Task 2: Using the `ising_tn` constructor and the `contract_tensornetwork` function,   verify these formulas. As you scale the system size, both results should approach the formula 
  $f(\beta) = -\frac{1}{\beta}\lim_{N \rightarrow \infty}\frac{1}{N}log(Z_{N}) = -\frac{1}{\beta}\ln(2 \cosh(\beta))
for the free energy density.

02-beliefpropagation.jl

When the tensor network has loops, it can become very computationally costly to contract. We are going to use belief propagation to approximately contract the ising partition function. On loopy networks.

Lets start with a single loop, i.e .a periodic boundary chain. 

Task 3: How does the BP value for the free energy density on a PBC chain compare to the exact value? How does this difference scale as a function of the chain of the length? 