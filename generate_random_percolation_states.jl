include("src/PERC.jl")



datasets = ["com-youtube.ungraph.txt","com-lj.ungraph.txt","com-orkut.ungraph.txt"]



directed = false

separator = "\t"
normalized = false
for ds in datasets
    gf = "graphs/"*ds
    g = load_graph(gf,directed,separator)
    perc = random_percolations(nv(g))
    ds_name = string(split(ds,".txt")[1])
    save_percolation_array(ds_name,perc)


end