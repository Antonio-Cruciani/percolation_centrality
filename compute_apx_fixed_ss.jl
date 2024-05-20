include("src/PERC.jl")


#Common parameters
epsilon_list = [0.07,0.05,0.01,0.005]
ss_save = [2,3,4,5]
delta = 0.1
run = 10


# Undirected
#=
datasets = ["00_hiv.txt","00_ego-fb-combined-N.txt","01_musae_facebook_edges.txt","02_email_enron.txt","03_ca_astroph.txt","07_large_twitch_edges.txt"]
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
            x = parallel_estimate_percolation_centrality_fixed_sample_size(g,perc,epsilon,delta)
            save_results_progressive_sampling(ds_name,"fixed_ss",x[1],x[2],x[3],ss_save[i],epsilon)
        end
        i+=1
    end
end
=#
# Directed
#datasets = ["11_soc_epinions.txt","12_soc_slashdot.txt","04_web_notredame.txt","06_web_google.txt","08_web_berkstan.txt"]
datasets = ["08_web_berkstan.txt"]

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
            x = parallel_estimate_percolation_centrality_fixed_sample_size(g,perc,epsilon,delta)
            save_results_progressive_sampling(ds_name,"fixed_ss",x[1],x[2],x[3],ss_save[i],epsilon)
        end
        i+=1
    end
end

epsilon_list = [0.001,0.0005]
ss_save = [6,7]

run = 1



# Undirected
datasets = ["00_hiv.txt","00_ego-fb-combined-N.txt","01_musae_facebook_edges.txt","02_email_enron.txt","03_ca_astroph.txt","07_large_twitch_edges.txt"]
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
            x = parallel_estimate_percolation_centrality_fixed_sample_size(g,perc,epsilon,delta)
            save_results_progressive_sampling(ds_name,"fixed_ss",x[1],x[2],x[3],ss_save[i],epsilon)
        end
        i+=1
    end
end

# Directed
datasets = ["11_soc_epinions.txt","12_soc_slashdot.txt","04_web_notredame.txt","06_web_google.txt","08_web_berkstan.txt"]
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
            x = parallel_estimate_percolation_centrality_fixed_sample_size(g,perc,epsilon,delta)
            save_results_progressive_sampling(ds_name,"fixed_ss",x[1],x[2],x[3],ss_save[i],epsilon)
        end
        i+=1
    end
end