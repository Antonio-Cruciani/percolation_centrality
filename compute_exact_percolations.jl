include("src/PERC.jl")


datasets = ["01_enron.txt"]

directed = false

separator = "\t"
normalized = true
for ds in datasets
    gf = "graphs/"*ds
    g = load_graph(gf,directed,separator)
    perc = random_percolations(nv(g))
    ds_name = string(split(ds,".txt")[1])
    save_percolation_array(ds_name,perc)
    x = parallel_percolation_centrality(g,perc,normalized)
    save_results(ds_name,"exact",x[1],x[2])
end
