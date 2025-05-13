include("src/PERC.jl")

graphs_path = "graphs/"
percolation_path = "percolation_states/"
output_path = ""



graphs_path = "graphs/"
percolation_path = "percolation_states/"
output_path = ""

directed = false
separator = "\t"


#"07_large_twitch_edges.txt"
datasets = ["01_musae_facebook_edges.txt","02_email_enron.txt","03_ca_astroph.txt"]

for ds in datasets
    gf = graphs_path*ds
    g = load_graph(gf,directed,separator)

    x = random_percolations(nv(g))
    ds_name = string(split(ds,".txt")[1])*"_unif"
    save_percolation_array(ds_name,x)
    @info("Computing Ground Truth percolation centrality for "*ds_name)
    flush(stderr)
    y = parallel_percolation_centrality_new_target(g,x)
    save_results_new(ds_name,"exact_target",y[1],y[6],y[5])
end

directed = true

datasets = ["15_cit_hepph.txt" ,"14_p2p_gnutella31.txt","11_soc_epinions.txt","12_soc_slashdot.txt","04_web_notredame.txt","06_web_google.txt"]
for ds in datasets
    gf = graphs_path*ds
    g = load_graph(gf,directed,separator)
    x = random_percolations(nv(g))
    ds_name = string(split(ds,".txt")[1])*"_unif"
    save_percolation_array(ds_name,x)
    @info("Computing Ground Truth percolation centrality for "*ds_name)
    flush(stderr)
    y = parallel_percolation_centrality_new_target(g,x)
    save_results_new(ds_name,"exact_target",y[1],y[6],y[5])
end