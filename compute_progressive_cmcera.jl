include("src/PERC.jl")


datasets = ["01_email_enron.txt"]

directed = false
epsilon_list = [0.1,0.07,0.05,0.01,0.005]
ss_save = [1,2,3,4,5]
delta = 0.1
run = 10
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