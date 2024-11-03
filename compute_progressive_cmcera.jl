include("src/PERC.jl")

#Common parameters
epsilon_list = [0.1,0.07,0.05,0.01,0.005]
ss_save = [1,2,3,4,5]
#ss_save = [6,7,8]
#epsilon_list = [0.0025,0.0001,0.00005]

delta = 0.1
run = 10

# Undirected
datasets = ["01_musae_facebook_edges.txt","02_email_enron.txt","03_ca_astroph.txt","07_large_twitch_edges.txt","10_flickr.txt","com-youtube.ungraph.txt","com-lj.ungraph.txt","com-orkut.ungraph.txt"]
#datasets = ["com-youtube.ungraph.txt","com-lj.ungraph.txt","com-orkut.ungraph.txt"]

directed = false
separator = "\t"

for ds in datasets
    gf = "graphs/"*ds
    g = load_graph(gf,directed,separator)
    perc = read_percolation("percolation_states/"*ds)
    ds_name = string(split(ds,".txt")[1])
    i = 1
    for epsilon in epsilon_list
        for _ in 1:run
            x = parallel_estimate_percolation_centrality_lock(g,perc,epsilon,delta)
            save_results_progressive_sampling(ds_name,"cmcera",x[1],x[2],x[4],ss_save[i],epsilon,x[5])
        end
        i+=1
    end
end

#=
# Directed
datasets = ["15_cit_hepph.txt","14_p2p_gnutella31.txt","11_soc_epinions.txt","12_soc_slashdot.txt","04_web_notredame.txt","06_web_google.txt","08_web_berkstan.txt","05_wiki_talk.txt","09_italian_twitter.txt","13_soc_pokec.txt"]

directed = true
separator = "\t"

for ds in datasets
    gf = "graphs/"*ds
    g = load_graph(gf,directed,separator)
    perc = read_percolation("percolation_states/"*ds)
    ds_name = string(split(ds,".txt")[1])
    i = 1
    for epsilon in epsilon_list
        for _ in 1:run
            @info("Computing Apx percolation centrality for "*ds_name)
            flush(stderr)
            x = parallel_estimate_percolation_centrality_lock(g,perc,epsilon,delta)
            save_results_progressive_sampling(ds_name,"cmcera",x[1],x[2],x[4],ss_save[i],epsilon,x[5])
        end
        i+=1
    end
end

=#