include("src/PERC.jl")

# Undirected
#=
datasets = ["00_hiv.txt","00_ego-fb-combined-N.txt","01_musae_facebook_edges.txt","02_email_enron.txt","03_ca_astroph.txt","07_large_twitch_edges.txt"]
datasets = []
directed = false

separator = "\t"
normalized = false
for ds in datasets
    gf = "graphs/"*ds
    g = load_graph(gf,directed,separator)
    perc = random_percolations(nv(g))
    ds_name = string(split(ds,".txt")[1])
    save_percolation_array(ds_name,perc)
    @info("Computing Ground Truth percolation centrality for "*ds_name)
    flush(stderr)
    
    x = parallel_percolation_centrality(g,perc,normalized)
    save_results(ds_name,"exact",x[1],x[2])
end
=#
# Directed
datasets = ["11_soc_epinions.txt","12_soc_slashdot.txt","04_web_notredame.txt","05_wiki_talk.txt","06_web_google.txt","08_web_berkstan.txt","09_italian_twitter.txt","13_soc_pokec.txt"]
directed = true

separator = "\t"
normalized = false
for ds in datasets
    gf = "graphs/"*ds
    g = load_graph(gf,directed,separator)
    perc = random_percolations(nv(g))
    ds_name = string(split(ds,".txt")[1])
    save_percolation_array(ds_name,perc)
    @info("Computing Ground Truth percolation centrality for "*ds_name)
    flush(stderr)
    x = parallel_percolation_centrality(g,perc,normalized)
    save_results(ds_name,"exact",x[1],x[2])
end