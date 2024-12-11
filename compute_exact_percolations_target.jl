include("src/PERC.jl")

# Undirected


datasets = ["00_hiv.txt","00_ego-fb-combined-N.txt"]
directed = false

separator = "\t"
normalized = false
percolation_path = "percolation_states/"
graphs_path = "graphs/"
for ds in datasets
    gf = graphs_path*ds
    g = load_graph(gf,directed,separator)
    perc = read_percolation(percolation_path*ds)
    ds_name = string(split(ds,".txt")[1])
    @info("Computing Ground Truth percolation centrality for "*ds_name)
    flush(stderr)
    
    x = parallel_percolation_centrality_target(g,perc,normalized)
    save_results(ds_name,"exact_target",x[1],x[2])
end
#=
# Directed
#datasets = ["15_cit_hepph.txt","14_p2p_gnutella31.txt","11_soc_epinions.txt","12_soc_slashdot.txt"]
datasets = ["09_italian_twitter.txt","13_soc_pokec.txt"]

directed = true

separator = "\t"
normalized = false
for ds in datasets
    gf = graphs_path*ds
    g = load_graph(gf,directed,separator)
    perc = read_percolation(percolation_path*ds)
    ds_name = string(split(ds,".txt")[1])
    @info("Computing Ground Truth percolation centrality for "*ds_name)
    flush(stderr)
    x = parallel_percolation_centrality_target(g,perc,normalized)
    save_results(ds_name,"exact_target",x[1],x[2])
end


=#
