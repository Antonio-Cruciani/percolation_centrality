include("src/PERC.jl")

# Undirected





#epsilon_list = [0.1,0.07,0.05,0.01,0.005,0.001,0.0005]
sample_size_list = [1000,5000,10000,50000,100000]
ss_save = [1,2,3,4,5,6,7]
delta = 0.05

run = 5
graphs_path = "graphs/"
percolation_path = "percolation_states/"
output_path = ""
component_size = 50
#,"07_large_twitch_edges.txt"
datasets = [
"youtube_10000_edges.txt" 
]

directed = false
separator = " "

for ds in datasets
    gf = graphs_path*ds
    g = load_graph(gf,directed,separator)
    ds_name = string(split(ds,".txt")[1])

    x = read_percolation(percolation_path*ds_name*".txt")

    @info("Computing Ground Truth percolation centrality for "*ds_name)
    flush(stderr)
    y = parallel_percolation_centrality_new_target(g,x)
    save_results_new(ds_name,"exact_target",y[1],y[6],y[5])
    #save_results(ds_name,"exact_target_unnormalized",x[2],x[3])
end




