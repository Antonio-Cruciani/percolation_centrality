include("src/PERC.jl")




#epsilon_list = [0.1,0.07,0.05,0.01,0.005,0.001,0.0005]
#500000,1000000
sample_size_list = [1000,5000,10000,50000,100000]
ss_save = [1,2,3,4,5,6,7]
delta = 0.05

run = 10
graphs_path = "graphs/"
percolation_path = "percolation_states/"
output_path = ""
component_size = 50
# ,"07_large_twitch_edges.txt"
datasets = [ "10_flickr.txt"]

#datasets = ["00_hiv.txt"]

directed = false
separator = "\t"

for ds in datasets
    gf = graphs_path*ds
    g = load_graph(gf,directed,separator)
    #h =g
    h,x = plant_random_initiators(g,component_size)
    ds_name = string(split(ds,".txt")[1])*"_rnd_init_"*string(component_size)
    println(percolation_path*ds_name*".txt")
    #x = read_percolation(percolation_path*ds_name*".txt")

    save_percolation_array(ds_name,x)
    #save_graph(h,ds_name,separator)
    @info("Computing Ground Truth percolation centrality for "*ds_name)
    flush(stderr)
    y = parallel_percolation_centrality_new_target(h,x)
    save_results_new(ds_name,"exact_target",y[1],y[6],y[5])
    #save_results(ds_name,"exact_target_unnormalized",x[2],x[3])
end


datasets = [ "08_web_berkstan.txt","13_soc_pokec.txt"]

#datasets = ["00_hiv.txt"]

directed = true
separator = "\t"

for ds in datasets
    gf = graphs_path*ds
    g = load_graph(gf,directed,separator)
    #h =g
    h,x = plant_random_initiators(g,component_size)
    ds_name = string(split(ds,".txt")[1])*"_rnd_init_"*string(component_size)
    println(percolation_path*ds_name*".txt")
    #x = read_percolation(percolation_path*ds_name*".txt")

    save_percolation_array(ds_name,x)
    #save_graph(h,ds_name,separator)
    @info("Computing Ground Truth percolation centrality for "*ds_name)
    flush(stderr)
    y = parallel_percolation_centrality_new_target(h,x)
    save_results_new(ds_name,"exact_target",y[1],y[6],y[5])
    #save_results(ds_name,"exact_target_unnormalized",x[2],x[3])
end


#graphs_path = "components/"

# Non uniform
#=
for ds in datasets
    ds_name = string(split(ds,".txt")[1])*"_rnd_init_"*string(component_size)
    gf = graphs_path*ds
    g = load_graph(gf,directed,separator)
    perc = read_percolation(percolation_path*ds_name*".txt")
    i =1
    for sample_size in sample_size_list
        for _ in 1:run
            y = parallel_estimate_percolation_centrality_non_uniform(g,perc,sample_size)
            save_results_progressive_sampling(ds_name,"non_uniform",y[1],y[2],y[4],ss_save[i],0.0,y[5],output_path)
            #Uniform
            h = parallel_estimate_percolation_centrality_uniform(g,perc,sample_size)
            save_results_progressive_sampling(ds_name,"uniform",h[1],h[2],h[4],ss_save[i],0.0,h[5],output_path)
        end
        i+=1
    end


end

# Uniform
for ds in datasets
    ds_name = string(split(ds,".txt")[1])*"_rnd_init_"*string(component_size)
    gf = graphs_path*ds
    g = load_graph(gf,directed,separator)
    perc = read_percolation(percolation_path*ds_name*".txt")
    i =1
    for sample_size in sample_size_list
        for _ in 1:run
            y = parallel_estimate_percolation_centrality_uniform(g,perc,sample_size)
            save_results_progressive_sampling(ds_name,"uniform",y[1],y[2],y[4],ss_save[i],0.0,y[5],output_path)
        end
        i+=1
    end


end

#graphs_path = "graphs/"




datasets = ["15_cit_hepph.txt" ,"14_p2p_gnutella31.txt","11_soc_epinions.txt","12_soc_slashdot.txt","04_web_notredame.txt","06_web_google.txt"]
#datasets = ["15_cit_hepph.txt" ]
separator = "\t"

directed = true

for ds in datasets
    gf = graphs_path*ds
    g = load_graph(gf,directed,separator)
    h,x = plant_random_initiators(g,component_size)
    ds_name = string(split(ds,".txt")[1])*"_rnd_init_"*string(component_size)
    save_percolation_array(ds_name,x)
    #save_graph(h,ds_name,separator)
    #nods = [nv(g)+i for i in 1:nv(h)-nv(g)]
    @info("Computing Ground Truth percolation centrality for "*ds_name)
    y = parallel_percolation_centrality_new_target(h,x)
    save_results_new(ds_name,"exact_target",y[1],y[6],y[5])
    #save_results(ds_name,"exact_target_unnormalized",x[2],x[3])
end

#graphs_path = "components/"
# Non uniform
for ds in datasets
    ds_name = string(split(ds,".txt")[1])*"_rnd_init_"*string(component_size)
    gf = graphs_path*ds
    g = load_graph(gf,directed,separator)
    perc = read_percolation(percolation_path*ds_name*".txt")
    i =1
    for sample_size in sample_size_list
        for _ in 1:run
            y = parallel_estimate_percolation_centrality_non_uniform(g,perc,sample_size)
            save_results_progressive_sampling(ds_name,"non_uniform",y[1],y[2],y[4],ss_save[i],0.0,y[5],output_path)
            #Uniform
            h = parallel_estimate_percolation_centrality_uniform(g,perc,sample_size)
            save_results_progressive_sampling(ds_name,"uniform",h[1],h[2],h[4],ss_save[i],0.0,h[5],output_path)
        end
        i+=1
    end


end
=#

#=
# Uniform
for ds in datasets
    ds_name = string(split(ds,".txt")[1])*"_rnd_init_"*string(component_size)
    gf = graphs_path*ds
    g = load_graph(gf,directed,separator)
    perc = read_percolation(percolation_path*ds_name*".txt")
    i =1
    for sample_size in sample_size_list
        for _ in 1:run
            y = parallel_estimate_percolation_centrality_uniform(g,perc,sample_size)
            save_results_progressive_sampling(ds_name,"uniform",y[1],y[2],y[4],ss_save[i],0.0,y[5],output_path)
        end
        i+=1
    end


end
=#
#=
# Non Uniform Sampling
graphs_path = "components/"

directed = false
separator = "\t"
#datasets = ["01_musae_facebook_edges.txt","02_email_enron.txt", "03_ca_astroph.txt"]
#datasets = ["07_large_twitch_edges.txt"]
datasets = ["01_musae_facebook_edges.txt","02_email_enron.txt", "03_ca_astroph.txt","07_large_twitch_edges.txt"]


for ds in datasets
    ds_name = string(split(ds,".txt")[1])*"_lcc_"*string(component_size)
    gf = graphs_path*ds_name*".txt"
    g = load_graph(gf,directed,separator)
    perc = read_percolation(percolation_path*ds_name*".txt")
    i =1
    for sample_size in sample_size_list
        for _ in 1:run
            y = parallel_estimate_percolation_centrality_non_uniform(g,perc,sample_size)
            save_results_progressive_sampling(ds_name,"non_uniform",y[1],y[2],y[4],ss_save[i],0.0,y[5],output_path)
        end
        i+=1
    end


end

=#
#datasets = ["15_cit_hepph.txt" ,"14_p2p_gnutella31.txt","11_soc_epinions.txt"]
#datasets = ["12_soc_slashdot.txt","04_web_notredame.txt","06_web_google.txt","08_web_berkstan.txt"]

#datasets = ["15_cit_hepph.txt" ,"14_p2p_gnutella31.txt","11_soc_epinions.txt","12_soc_slashdot.txt","04_web_notredame.txt","06_web_google.txt"]
datasets = ["15_cit_hepph.txt" ]
graphs_path = "components/"
#=
directed = true
separator = "\t"

for ds in datasets
    ds_name = string(split(ds,".txt")[1])*"_lcc_"*string(component_size)
    gf = graphs_path*ds_name*".txt"
    g = load_graph(gf,directed,separator)
    perc = read_percolation(percolation_path*ds_name*".txt")
    i =1
    for sample_size in sample_size_list
        for _ in 1:run
            y = parallel_estimate_percolation_centrality_non_uniform(g,perc,sample_size)
            save_results_progressive_sampling(ds_name,"non_uniform",y[1],y[2],y[4],ss_save[i],0.0,y[5],output_path)
        end
        i+=1
    end


end
=#
#=
# Uniform

directed = false
separator = "\t"
#datasets = ["01_musae_facebook_edges.txt","02_email_enron.txt", "03_ca_astroph.txt"]
datasets = ["01_musae_facebook_edges.txt","02_email_enron.txt", "03_ca_astroph.txt","07_large_twitch_edges.txt"]



for ds in datasets
    ds_name = string(split(ds,".txt")[1])*"_lcc_"*string(component_size)
    gf = graphs_path*ds_name*".txt"
    g = load_graph(gf,directed,separator)
    perc = read_percolation(percolation_path*ds_name*".txt")
    i =1
    for sample_size in sample_size_list
        for _ in 1:run
            y = parallel_estimate_percolation_centrality_uniform(g,perc,sample_size)
            save_results_progressive_sampling(ds_name,"uniform",y[1],y[2],y[4],ss_save[i],0.0,y[5],output_path)
        end
        i+=1
    end


end
=#

#datasets = ["15_cit_hepph.txt" ,"14_p2p_gnutella31.txt","11_soc_epinions.txt"]
#,"08_web_berkstan.txt"
#datasets = ["15_cit_hepph.txt" ,"14_p2p_gnutella31.txt","11_soc_epinions.txt","12_soc_slashdot.txt","04_web_notredame.txt","06_web_google.txt"]
datasets = ["15_cit_hepph.txt" ]
#=
directed = true
separator = "\t"

for ds in datasets
    ds_name = string(split(ds,".txt")[1])*"_lcc_"*string(component_size)
    gf = graphs_path*ds_name*".txt"
    g = load_graph(gf,directed,separator)
    perc = read_percolation(percolation_path*ds_name*".txt")
    i =1
    for sample_size in sample_size_list
        for _ in 1:run
            y = parallel_estimate_percolation_centrality_uniform(g,perc,sample_size)
            save_results_progressive_sampling(ds_name,"uniform",y[1],y[2],y[4],ss_save[i],0.0,y[5],output_path)
        end
        i+=1
    end


end
=#
