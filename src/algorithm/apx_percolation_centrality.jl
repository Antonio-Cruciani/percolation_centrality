
function estimate_percolation_centrality(g,percolation_states::Array{Float64},epsilon::Float64,delta::Float64, alpha_sampling::Float64 = 1.0, mc_trials::Int64 = 25, empirical_peeling_a::Float64 = 2.0,sample_size_diam::Int64 = 256 )
    n::Int64 = nv(g)
    m::Int64 = ne(g)
    directed::Bool = is_directed(g)
    sg::static_graph = static_graph(adjacency_list(g),incidency_list(g))
    q::Array{Int64} = zeros(Int64,n)
    ball::Array{Int16} = zeros(Int16,n)
    n_paths::Array{Int64} = zeros(Int64,n)
    dist::Array{Int64} = zeros(Int64,n)
    pred::Array{Array{Int64}} = [Array{Int64}([]) for _ in 1:n]
    num_paths::Array{Int64} = zeros(Int64,n)
    percolation_centrality::Array{Float64} = zeros(Float64,n)
    wimpy_variance::Array{Float64} = zeros(Float64,n)
    shortest_path_length::Array{Int64} = zeros(Int64,n+1)
    percolated_path_length::Array{Float64} = zeros(Float64,n+1)
    mcrade::Array{Float64} = zeros(Float64,(n+1)*mc_trials)
    tmp_perc_states::Array{Float64} = copy(percolation_states)
    percolation_data::Tuple{Float64,Array{Float64}} = percolation_differences(sort(tmp_perc_states),n)
    betweenness::Array{Float64} = zeros(Float64,n)
    start_time::Float64 = time()
    emp_wimpy_node::Float64 = 0.0
    min_inv_wimpy_node::Float64 = 0.0
    max_perc::Float64 = 0.0
    max_wv::Float64 = 0.0
    # Partitions 
    number_of_non_empty_partitions::Int64 = 0
    non_empty_partitions::Dict{Int64,Int64} = Dict{Int64,Int64}()
    partitions_ids_map::Dict{Int64,Int64} = Dict{Int64,Int64}()
    partition_index::Array{Int64} = zeros(n)
    part_idx::Int64 = 1
    # Diameter approximation using Propagate/RandomBFS
    # Using Random BFS
    println("Approximating diameter using Random BFS algorithm")
    flush(stdout)
    diam,time_diam = random_bfs(g,sample_size_diam)
    finish_diam::String = string(round(time_diam; digits=4))
    println("Estimated diameter "*string(diam)*" in "*finish_diam*" seconds")
    flush(stdout)
    # Sample size related
    omega::Float64 = 1000
    max_num_samples::Float64 = 0.0
    tau::Int64 = trunc(Int64,max(1. / epsilon * (log(1. / delta)) , 100.))
    tau =  trunc(Int64,max(tau,2*(diam -1) * (log(1. / delta))) )
    println("Bootstrap phase "*string(tau)*" iterations")
    flush(stdout)
    for _ in 1:tau
        _random_path!(sg,n,q,ball,n_paths,dist,pred,num_paths,percolation_centrality,wimpy_variance,percolation_states,percolation_data,shortest_path_length,percolated_path_length,mcrade,mc_trials,alpha_sampling,betweenness,true)
    end
    println("Empirical peeling phase:")
    flush(stdout)
    for i in 1:n
        max_perc = max(max_perc,percolation_centrality[i])
        max_wv = max(max_wv,wimpy_variance[i])
        emp_w_node = wimpy_variance[i] * 1. /tau
        min_inv_w_node = min(1. /emp_w_node,tau)
        node_partition_idx = trunc(Int,log(min_inv_w_node)/log(empirical_peeling_a)+1)
        partition_index[i] = node_partition_idx
        if haskey(non_empty_partitions,node_partition_idx)
            non_empty_partitions[node_partition_idx] += 1
        else
            non_empty_partitions[node_partition_idx] = 1
        end
    end
    number_of_non_empty_partitions = 0
    for key in keys(non_empty_partitions)
        partitions_ids_map[key] = part_idx
        part_idx+=1
        number_of_non_empty_partitions+=1
    end
    # To do the rest
    #WIP

    #= Test
    test_couples = []
    for i in 1:n
        for j in i:n
            if i!=j
                push!(test_couples,(i,j))
            end
        end
    end
    println(" number of couples ",length(test_couples))
    #println(test_couples)
    for tc in test_couples
        _random_path!(tc[1],tc[2],betweenness,sg,n,q,ball,n_paths,dist,pred,num_paths,percolation_centrality,wimpy_variance,percolation_states,percolation_data,shortest_path_length,percolated_path_length,mcrade,mc_trials,alpha_sampling)
    end
    #println(betweenness)
    
    return (betweenness)
    =#

    println("---------------")
    println(percolation_centrality)
    println(wimpy_variance)
    println("-------------------------")
    println(shortest_path_length)
    println(percolated_path_length)

end