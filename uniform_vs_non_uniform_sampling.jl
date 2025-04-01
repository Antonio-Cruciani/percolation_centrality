include("src/PERC.jl")


datasets = ["01_musae_facebook_edges.txt","02_email_enron.txt", "03_ca_astroph.txt"]


epsilon_list = [0.1,0.07,0.05,0.01,0.005,0.001,0.0005]

ss_save = [1,2,3,4,5,6,7]
delta = 0.05

run = 5

graphs_path = "graphs/"
percolation_path = "percolation_states/"
output_path = ""
component_size = 50
#=
directed = false
separator = "\t"

for ds in datasets
    gf = graphs_path*ds
    g = load_graph(gf,directed,separator)
    h,x = attach_gadget(g,component_size)
    ds_name = string(split(ds,".txt")[1])*"_"*string(component_size)
    save_percolation_array(ds_name,x)
    save_graph(h,ds_name,separator)
    @info("Computing Ground Truth percolation centrality for "*ds_name)
    flush(stderr)
    y = parallel_percolation_centrality_new_target(h,x)
    save_results_new(ds_name,"exact_target",y[1],y[6],y[5])
    #save_results(ds_name,"exact_target_unnormalized",x[2],x[3])
end

datasets = ["15_cit_hepph.txt" ,"14_p2p_gnutella31.txt","11_soc_epinions.txt"]

directed = true

for ds in datasets
    gf = graphs_path*ds
    g = load_graph(gf,directed,separator)
    h,x = attach_gadget(g,component_size)
    ds_name = string(split(ds,".txt")[1])*"_"*string(component_size)
    save_percolation_array(ds_name,x)
    save_graph(h,ds_name,separator)
    @info("Computing Ground Truth percolation centrality for "*ds_name)
    flush(stderr)
    y = parallel_percolation_centrality_new_target(h,x)
    save_results_new(ds_name,"exact_target",y[1],y[6],y[5])
    #save_results(ds_name,"exact_target_unnormalized",x[2],x[3])
end
=#

# Non Uniform Sampling
graphs_path = "components/"

directed = false
separator = "\t"


for ds in datasets
    ds_name = string(split(ds,".txt")[1])*"_"*string(component_size)
    gf = graphs_path*ds_name*".txt"
    g = load_graph(gf,directed,separator)
    perc = read_percolation(percolation_path*ds_name*".txt")
    i =1
    for epsilon in epsilon_list
        for _ in 1:run
            y = parallel_estimate_percolation_centrality_new_lock(g,perc,epsilon,delta)
            save_results_progressive_sampling(ds_name,"psilvan",y[1],y[2],y[4],ss_save[i],epsilon,y[5],output_path)
        end
        i+=1
    end


end


datasets = ["15_cit_hepph.txt" ,"14_p2p_gnutella31.txt","11_soc_epinions.txt"]



directed = true
separator = "\t"

for ds in datasets
    ds_name = string(split(ds,".txt")[1])*"_"*string(component_size)
    gf = graphs_path*ds_name*".txt"
    g = load_graph(gf,directed,separator)
    perc = read_percolation(percolation_path*ds_name*".txt")
    i =1
    for epsilon in epsilon_list
        for _ in 1:run
            y = parallel_estimate_percolation_centrality_new_lock(g,perc,epsilon,delta)
            save_results_progressive_sampling(ds_name,"psilvan",y[1],y[2],y[4],ss_save[i],epsilon,y[5],output_path)
        end
        i+=1
    end


end

# Uniform Sampling Rade
vc_bound = false
datasets = ["01_musae_facebook_edges.txt","02_email_enron.txt", "03_ca_astroph.txt"]

directed = false
separator = "\t"

for ds in datasets
    ds_name = string(split(ds,".txt")[1])*"_"*string(component_size)
    gf = graphs_path*ds_name*".txt"
    g = load_graph(gf,directed,separator)
    perc = read_percolation(percolation_path*ds_name*".txt")
    i =1
    for epsilon in epsilon_list
        for _ in 1:run
            y = parallel_estimate_percolation_centrality_era_new(g,perc,epsilon,delta)
            save_results_progressive_sampling(ds_name,"era",y[1],y[2][end],y[5],ss_save[i],y[4],-1.0,output_path)

        end
        i+=1
    end


end




datasets = ["15_cit_hepph.txt" ,"14_p2p_gnutella31.txt","11_soc_epinions.txt"]



directed = true
separator = "\t"


for ds in datasets
    ds_name = string(split(ds,".txt")[1])*"_"*string(component_size)
    gf = graphs_path*ds_name*".txt"
    g = load_graph(gf,directed,separator)
    perc = read_percolation(percolation_path*ds_name*".txt")
    i =1
    for epsilon in epsilon_list
        for _ in 1:run
            y = parallel_estimate_percolation_centrality_era_new(g,perc,epsilon,delta)
            save_results_progressive_sampling(ds_name,"era",y[1],y[2][end],y[5],ss_save[i],y[4],-1.0,output_path)
        end
        i+=1
    end


end

# Uniform Fixed Sample Size


datasets = ["01_musae_facebook_edges.txt","02_email_enron.txt", "03_ca_astroph.txt"]

directed = false
separator = "\t"
for ds in datasets
    ds_name = string(split(ds,".txt")[1])*"_"*string(component_size)
    gf = graphs_path*ds_name*".txt"
    g = load_graph(gf,directed,separator)
    perc = read_percolation(percolation_path*ds_name*".txt")
    i =1
    for epsilon in epsilon_list
        for _ in 1:run
            y = parallel_estimate_percolation_centrality_fixed_sample_size_new(g,perc,epsilon,delta)
            save_results_progressive_sampling(ds_name,"fixed_ss",y[1],y[2],y[3],ss_save[i],epsilon,-1.0,output_path)

        end
        i+=1
    end


end


datasets = ["15_cit_hepph.txt" ,"14_p2p_gnutella31.txt","11_soc_epinions.txt"]


directed = true
separator = "\t"


for ds in datasets
    ds_name = string(split(ds,".txt")[1])*"_"*string(component_size)
    gf = graphs_path*ds_name*".txt"
    g = load_graph(gf,directed,separator)
    perc = read_percolation(percolation_path*ds_name*".txt")
    i =1
    for epsilon in epsilon_list
        for _ in 1:run
            y = parallel_estimate_percolation_centrality_fixed_sample_size_new(g,perc,epsilon,delta)
            save_results_progressive_sampling(ds_name,"fixed_ss",y[1],y[2],y[3],ss_save[i],epsilon,-1.0,output_path)
        end
        i+=1
    end


end