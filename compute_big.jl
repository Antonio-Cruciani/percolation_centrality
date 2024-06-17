include("src/PERC.jl")
epsilon_list = [0.1,0.07,0.05,0.01,0.005]
ss_save = [1,2,3,4,5]
delta = 0.1
run = 10


# Percolation states generation

# Undirected
datasets = ["10_flickr.txt"]



directed = false
separator = "\t"

for ds in datasets
    gf = "graphs/"*ds
    g = load_graph(gf,directed,separator)
    perc = random_percolations(nv(g))
    ds_name = string(split(ds,".txt")[1])
    save_percolation_array(ds_name,perc)
end

# Directed
datasets = ["05_wiki_talk.txt","09_italian_twitter.txt","13_soc_pokec.txt"]

directed = true
separator = "\t"

for ds in datasets
    gf = "graphs/"*ds
    g = load_graph(gf,directed,separator)
    perc = random_percolations(nv(g))
    ds_name = string(split(ds,".txt")[1])
    save_percolation_array(ds_name,perc)
end

# Estimation phase

# c-MCERA

# Undirected

datasets = ["10_flickr.txt"]



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
            save_results_progressive_sampling(ds_name,"cmcera",x[1],x[2],x[4],ss_save[i],epsilon)
        end
        i+=1
    end
end


# Directed
datasets = ["05_wiki_talk.txt","09_italian_twitter.txt","13_soc_pokec.txt"]

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


# ERA

# Undirected

datasets = ["10_flickr.txt"]



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
            @info("Computing Apx percolation centrality for "*ds_name)
            flush(stderr)
            x = parallel_estimate_percolation_centrality_era(g,perc,epsilon,delta,0,1.2,256,true)
            save_results_progressive_sampling(ds_name,"era",x[1],x[2][end],x[5],ss_save[i],x[4])
        end
        i+=1
    end
end





# Directed
datasets = ["05_wiki_talk.txt","09_italian_twitter.txt","13_soc_pokec.txt"]

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
            x = parallel_estimate_percolation_centrality_era(g,perc,epsilon,delta,0,1.2,256,true)
            save_results_progressive_sampling(ds_name,"era",x[1],x[2][end],x[5],ss_save[i],x[4])
        end
        i+=1
    end
end

# Fixed Sample Size
# Undirected

datasets = ["10_flickr.txt"]



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
            @info("Computing Apx percolation centrality for "*ds_name)
            flush(stderr)
            x = parallel_estimate_percolation_centrality_fixed_sample_size(g,perc,epsilon,delta)
            save_results_progressive_sampling(ds_name,"fixed_ss",x[1],x[2],x[3],ss_save[i],epsilon)
        end
        i+=1
    end
end


# Directed
datasets = ["05_wiki_talk.txt","09_italian_twitter.txt","13_soc_pokec.txt"]

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
            x = parallel_estimate_percolation_centrality_fixed_sample_size(g,perc,epsilon,delta)
            save_results_progressive_sampling(ds_name,"fixed_ss",x[1],x[2],x[3],ss_save[i],epsilon)
        end
        i+=1
    end
end
