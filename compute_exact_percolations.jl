include("src/PERC.jl")

# Undirected
#=
#datasets = ["03_ca_astroph.txt","07_large_twitch_edges.txt"]
datasets = ["00_hiv.txt","00_ego-fb-combined-N.txt"]
directed = false
graphs_path = "graphs/"
percolation_path = "percolation_states/"
separator = "\t"
normalized = true
for ds in datasets
    gf = graphs_path*ds
    g = load_graph(gf,directed,separator)
    perc = read_percolation(percolation_path*ds)
    ds_name = string(split(ds,".txt")[1])
    @info("Computing Ground Truth percolation centrality for "*ds_name)
    flush(stderr)
    
    x = parallel_percolation_centrality(g,perc,normalized)
    save_results(ds_name,"exact",x[1],x[2])
end
=#
# Directed
#datasets = ["04_web_notredame.txt","06_web_google.txt","08_web_berkstan.txt"]
datasets = ["05_wiki_talk","13_soc_pokec.txt"]


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
    x = parallel_percolation_centrality(g,perc,normalized)
    save_results(ds_name,"exact",x[1],x[2])
end


