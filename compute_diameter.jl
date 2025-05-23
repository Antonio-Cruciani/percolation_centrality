include("src/PERC.jl")



graph_path = "graphs/"

# Undirected
datasets = ["01_musae_facebook_edges.txt","02_email_enron.txt","03_ca_astroph.txt","24_uselections.txt","23_twitter_pol.txt","22_obamacare.txt","21_brexit.txt","20_abortion.txt"]


graphs_path = "graphs/"
percolation_path = "percolation_states/"
separator = " "
directed = false
for ds in datasets
    gf = graphs_path*ds
    g = load_graph(gf,directed,separator)
    ds_name = string(split(ds,".txt")[1])
    @info("Computing Ground Truth Diameter for "*ds_name)
    x,y,t = parallel_random_bfs_rho(g,0)
    save_diameter(ds_name,x,y,t)
end





datasets = ["15_cit_hepph.txt" ,"14_p2p_gnutella31.txt","11_soc_epinions.txt","12_soc_slashdot.txt","04_web_notredame.txt","06_web_google.txt"]
graphs_path = "graphs/"
percolation_path = "percolation_states/"
separator = " "
directed = true
for ds in datasets
    gf = graphs_path*ds
    g = load_graph(gf,directed,separator)
    ds_name = string(split(ds,".txt")[1])
    @info("Computing Ground Truth Diameter for "*ds_name)
    x,y,t = parallel_random_bfs_rho(g,0)
    save_diameter(ds_name,x,y,t)
end

