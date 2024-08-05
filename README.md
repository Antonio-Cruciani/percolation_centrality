# Fast Approximation of Percolation Centrality

Repository associated to the paper "Fast Approximation of Percolation Centrality".

# Network file format
The input file must be a edge list, one line for each edge (a edge is a tuple `(u,v)`, where `u` and `v` are two vertices). While loading the edge list you must specify if the graph is directed or undirected. The two elements of an edge can be separated by any string, which can be specified while using the `load_graph` function (see below). In the following, we assume that this string is just a tab.



# Note on multi-threading
In order to properly run the multi-threading implementations of our approaches you need to set the number of threads that you desire to use with: 
```
 julia --threads 16
```
where the number indicates the number of thread to use, for more information we refer to the official documentation: [Starting Julia with multiple threads](https://docs.julialang.org/en/v1/manual/multi-threading/)


# How to reproduce the experiments in the paper
To reproduce all the experiment you need to: 

(i) Compute the exact percolation centrality values running the command `julia --threads <nthreads> compute_exact_percolations.jl`

(ii) Compute our algorithm's approximation, running `julia --threads <nthreads> compute_progressive_cmcera.jl`

(iii) Compute the Lima et al., progressive sampling algorith's approximation, running `julia --threads <nthreads> compute_progressive_era.jl`

(iv) Compute the Lima et al., progressive sampling algorith's approximation, running `julia --threads <nthreads> compute_apx_fixed_ss.jl`

where `<nthreads>` is the number of assigned threads.
All the results will be automatically saved in the `scores` and `times` folders.


