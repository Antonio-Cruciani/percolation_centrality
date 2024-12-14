include("src/PERC.jl")

# Undirected



#=
separator = "\t"
percolation_path = "percolation_states/"
graphs_path = "graphs/"

datasets = ["01_musae_facebook_edges.txt","02_email_enron.txt", "03_ca_astroph.txt","07_large_twitch_edges.txt" ,"10_flickr.txt"]
directed = false
for ds in datasets
    gf = graphs_path*ds
    g = load_graph(gf,directed,separator)
    perc = read_percolation(percolation_path*ds)
    ds_name = string(split(ds,".txt")[1])
    @info("Computing Ground Truth percolation centrality for "*ds_name)
    flush(stderr)
    
    x = parallel_percolation_centrality_target(g,perc)
    save_results(ds_name,"exact_target",x[1],x[3])
    save_results(ds_name,"exact_target_unnormalized",x[2],x[3])

end
=#
# Directed
#datasets = ["15_cit_hepph.txt","14_p2p_gnutella31.txt","11_soc_epinions.txt","12_soc_slashdot.txt"]
datasets = [ "15_cit_hepph.txt" ,"14_p2p_gnutella31.txt","11_soc_epinions.txt" ,"12_soc_slashdot.txt","04_web_notredame.txt" ,"06_web_google.txt","08_web_berkstan.txt"  ,"09_italian_twitter","05_wiki_talk.txt","13_soc_pokec.txt"]                 


directed = true

separator = "\t"
for ds in datasets
    gf = graphs_path*ds
    g = load_graph(gf,directed,separator)
    perc = read_percolation(percolation_path*ds)
    ds_name = string(split(ds,".txt")[1])
    @info("Computing Ground Truth percolation centrality for "*ds_name)
    flush(stderr)
    x = parallel_percolation_centrality_target(g,perc)
    save_results(ds_name,"exact_target",x[1],x[3])
    save_results(ds_name,"exact_target_unnormalized",x[2],x[3])
end


