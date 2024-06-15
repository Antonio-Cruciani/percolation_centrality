include("src/PERC.jl")

# Undirected

datasets = ["00_hiv.txt","00_ego-fb-combined-N.txt","01_musae_facebook_edges.txt","02_email_enron.txt","03_ca_astroph.txt","07_large_twitch_edges.txt"]
directed = false


separator = "\t"
normalized = false
for ds in datasets
    gf = "graphs/"*ds
    g = load_graph(gf,directed,separator)
    ds_name = string(split(ds,".txt")[1])
    @info("Computing Ground Truth distance based metrics for "*ds_name)
    flush(stderr)
    x = parallel_random_bfs_rho(g,0)
    save_diameter(ds_name,x[1],x[2],x[3])
end

# Directed
datasets = ["15_cit_hepph.txt","14_p2p_gnutella31.txt","11_soc_epinions.txt","12_soc_slashdot.txt","04_web_notredame.txt","06_web_google.txt","08_web_berkstan.txt"]
#datasets = ["06_web_google.txt","08_web_berkstan.txt","09_italian_twitter.txt","13_soc_pokec.txt","05_wiki_talk.txt"]


directed = true

separator = "\t"
normalized = false
for ds in datasets
    gf = "graphs/"*ds
    g = load_graph(gf,directed,separator)
    ds_name = string(split(ds,".txt")[1])
    @info("Computing Ground Truth distance based metrics for "*ds_name)
    flush(stderr)
    x = parallel_random_bfs_rho(g,0)
    save_diameter(ds_name,x[1],x[2],x[3])
end

