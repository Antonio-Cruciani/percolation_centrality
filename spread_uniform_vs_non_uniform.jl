include("src/PERC.jl")



sample_size_list = [1000,5000,10000,50000,100000]
ss_save = [1,2,3,4,5,6,7]
delta = 0.05

run = 10
graphs_path = "graphs/"
percolation_path = "percolation_states/"
output_path = ""
#datasets = ["00_hiv.txt"]

directed = false
separator = "\t"

mass = 4.0 
max_distance = typemax(Int64)


datasets = ["01_musae_facebook_edges.txt","02_email_enron.txt", "03_ca_astroph.txt","07_large_twitch_edges.txt"]

for ds in datasets
    gf = graphs_path*ds
    g = load_graph(gf,directed,separator)
    k  = trunc(Int64,ceil(log(nv(g))))
    @info("Simulating spreading from "*string(k)*"/"*string(nv(g))*" random nodes")
    x = simulate_spreading(g,k,max_distance,mass)
    ds_name = string(split(ds,".txt")[1])*"_e_log"
    save_percolation_array(ds_name,x)
    @info("Computing Ground Truth percolation centrality for "*ds_name)
    flush(stderr)
    y = parallel_percolation_centrality_new_target(g,x)
    save_results_new(ds_name,"exact_target_e_log",y[1],y[6],y[5])
    #save_results(ds_name,"exact_target_unnormalized",x[2],x[3])
end


# Non uniform 
for ds in datasets
    ds_name = string(split(ds,".txt")[1])*"_e_log"
    gf = graphs_path*ds
    g = load_graph(gf,directed,separator)
    perc = read_percolation(percolation_path*ds_name*".txt")
    i =1
    for sample_size in sample_size_list
        for _ in 1:run
            y = parallel_estimate_percolation_centrality_non_uniform(g,perc,sample_size)
            save_results_progressive_sampling(ds_name,"non_uniform_log",y[1],y[2],y[4],ss_save[i],0.0,y[5],output_path)
        end
        i+=1
    end


end

#Uniform

for ds in datasets
    ds_name = string(split(ds,".txt")[1])*"_e_log"
    gf = graphs_path*ds
    g = load_graph(gf,directed,separator)
    perc = read_percolation(percolation_path*ds_name*".txt")
    i =1
    for sample_size in sample_size_list
        for _ in 1:run
            y = parallel_estimate_percolation_centrality_uniform(g,perc,sample_size)
            save_results_progressive_sampling(ds_name,"uniform_log",y[1],y[2],y[4],ss_save[i],0.0,y[5],output_path)
        end
        i+=1
    end


end



datasets = ["15_cit_hepph.txt" ,"14_p2p_gnutella31.txt","11_soc_epinions.txt","12_soc_slashdot.txt","04_web_notredame.txt","06_web_google.txt"]

directed = true

for ds in datasets
    gf = graphs_path*ds
    g = load_graph(gf,directed,separator)
    k  = trunc(Int64,ceil(log(nv(g))))
    @info("Simulating spreading from "*string(k)*"/"*string(nv(g))*" random nodes")
    x = simulate_spreading(g,k,max_distance,mass)
    ds_name = string(split(ds,".txt")[1])*"_e_log"
    save_percolation_array(ds_name,x)
    @info("Computing Ground Truth percolation centrality for "*ds_name)
    flush(stderr)
    y = parallel_percolation_centrality_new_target(g,x)
    save_results_new(ds_name,"exact_target_e_log",y[1],y[6],y[5])
    #save_results(ds_name,"exact_target_unnormalized",x[2],x[3])
end

# Non uniform 
for ds in datasets
    ds_name = string(split(ds,".txt")[1])*"_e_log"
    gf = graphs_path*ds
    g = load_graph(gf,directed,separator)
    perc = read_percolation(percolation_path*ds_name*".txt")
    i =1
    for sample_size in sample_size_list
        for _ in 1:run
            y = parallel_estimate_percolation_centrality_non_uniform(g,perc,sample_size)
            save_results_progressive_sampling(ds_name,"non_uniform_log",y[1],y[2],y[4],ss_save[i],0.0,y[5],output_path)
        end
        i+=1
    end


end

# Uniform
for ds in datasets
    ds_name = string(split(ds,".txt")[1])*"_e_log"
    gf = graphs_path*ds
    g = load_graph(gf,directed,separator)
    perc = read_percolation(percolation_path*ds_name*".txt")
    i =1
    for sample_size in sample_size_list
        for _ in 1:run
            y = parallel_estimate_percolation_centrality_uniform(g,perc,sample_size)
            save_results_progressive_sampling(ds_name,"uniform_log",y[1],y[2],y[4],ss_save[i],0.0,y[5],output_path)
        end
        i+=1
    end


end


# Non Uniform Sampling
#=
directed = false
separator = "\t"
datasets = ["01_musae_facebook_edges.txt","02_email_enron.txt", "03_ca_astroph.txt","07_large_twitch_edges.txt"]



for ds in datasets
    ds_name = string(split(ds,".txt")[1])*"_e_log"
    gf = graphs_path*ds
    g = load_graph(gf,directed,separator)
    perc = read_percolation(percolation_path*ds_name*".txt")
    i =1
    for sample_size in sample_size_list
        for _ in 1:run
            y = parallel_estimate_percolation_centrality_non_uniform(g,perc,sample_size)
            save_results_progressive_sampling(ds_name,"non_uniform_log",y[1],y[2],y[4],ss_save[i],0.0,y[5],output_path)
        end
        i+=1
    end


end


#datasets = ["15_cit_hepph.txt" ,"14_p2p_gnutella31.txt","11_soc_epinions.txt"]
datasets = ["15_cit_hepph.txt" ,"14_p2p_gnutella31.txt","11_soc_epinions.txt","12_soc_slashdot.txt","04_web_notredame.txt","06_web_google.txt"]
=#
#datasets = ["08_web_berkstan.txt"]

#=
directed = true
separator = "\t"


for ds in datasets
    ds_name = string(split(ds,".txt")[1])*"_e_log"
    gf = graphs_path*ds
    g = load_graph(gf,directed,separator)
    perc = read_percolation(percolation_path*ds_name*".txt")
    i =1
    for sample_size in sample_size_list
        for _ in 1:run
            y = parallel_estimate_percolation_centrality_non_uniform(g,perc,sample_size)
            save_results_progressive_sampling(ds_name,"non_uniform_log",y[1],y[2],y[4],ss_save[i],0.0,y[5],output_path)
        end
        i+=1
    end


end
=#
# Uniform
#=
directed = false
separator = "\t"
#datasets = ["01_musae_facebook_edges.txt","02_email_enron.txt", "03_ca_astroph.txt"]
datasets = ["01_musae_facebook_edges.txt","02_email_enron.txt", "03_ca_astroph.txt","07_large_twitch_edges.txt"]



for ds in datasets
    ds_name = string(split(ds,".txt")[1])*"_e_log"
    gf = graphs_path*ds
    g = load_graph(gf,directed,separator)
    perc = read_percolation(percolation_path*ds_name*".txt")
    i =1
    for sample_size in sample_size_list
        for _ in 1:run
            y = parallel_estimate_percolation_centrality_uniform(g,perc,sample_size)
            save_results_progressive_sampling(ds_name,"uniform_log",y[1],y[2],y[4],ss_save[i],0.0,y[5],output_path)
        end
        i+=1
    end


end




#datasets = ["15_cit_hepph.txt" ,"14_p2p_gnutella31.txt","11_soc_epinions.txt"]
datasets = ["15_cit_hepph.txt" ,"14_p2p_gnutella31.txt","11_soc_epinions.txt","12_soc_slashdot.txt","04_web_notredame.txt","06_web_google.txt"]

=#
#datasets = ["08_web_berkstan.txt"]

#directed = true
#separator = "\t"
#=



for ds in datasets
    ds_name = string(split(ds,".txt")[1])*"_e_log"
    gf = graphs_path*ds
    g = load_graph(gf,directed,separator)
    perc = read_percolation(percolation_path*ds_name*".txt")
    i =1
    for sample_size in sample_size_list
        for _ in 1:run
            y = parallel_estimate_percolation_centrality_uniform(g,perc,sample_size)
            save_results_progressive_sampling(ds_name,"uniform_log",y[1],y[2],y[4],ss_save[i],0.0,y[5],output_path)
        end
        i+=1
    end


end
=#