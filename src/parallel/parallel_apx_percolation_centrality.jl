
function parallel_estimate_percolation_centrality(g,percolation_states::Array{Float64},epsilon::Float64,delta::Float64, alpha_sampling::Float64 = 0.1, mc_trials::Int64 = 25, empirical_peeling_a::Float64 = 2.0,sample_size_diam::Int64 = 256 )
    n::Int64 = nv(g)
    m::Int64 = ne(g)
    directed::Bool = is_directed(g)
    @info("----------------------------------------| Stats |--------------------------------------------------")
    @info("Analyzing graph")
    @info("Number of nodes "*string(n))
    @info("Number of edges "*string(m))
    @info("Directed ? "*string(directed))
    @info("Maximum Percolated state "*string(maximum(percolation_states)))
    @info("Minimum Percolated state "*string(minimum(percolation_states)))
    @info("Average Percolated state "*string(mean(percolation_states))*" std "*string(std(percolation_states)))
    @info("Using "*string(nthreads())* " Threads")
    @info("---------------------------------------------------------------------------------------------------")
    flush(stdout)
    ntasks = nthreads()
    sg::static_graph = static_graph(adjacency_list(g),incidency_list(g))
   
    percolation_centrality::Array{Array{Float64}} = [zeros(Float64,n) for _ in 1:ntasks]
    final_percolation_centrality::Array{Float64} = zeros(Float64,n)
    wimpy_variance::Array{Array{Float64}} = [zeros(Float64,n) for _ in 1:ntasks]
    final_wimpy_variance::Array{Float64} = zeros(Float64,n)
    shortest_path_length::Array{Array{Int64}} = [zeros(Int64,n+1) for _ in 1:ntasks]
    final_shortest_path_length::Array{Int64} = zeros(Int64,n+1)
    percolated_path_length::Array{Array{Float64}} = [zeros(Float64,n+1) for _ in 1:ntasks]
    final_percolated_path_length::Array{Float64} = zeros(Float64,n+1)
    mcrade::Array{Array{Float64}} = [zeros(Float64,(n+1)*mc_trials) for _ in 1:ntasks]
    final_mcrade::Array{Float64} = zeros(Float64,(n+1)*mc_trials)
    tmp_perc_states::Array{Float64} = copy(percolation_states)
    percolation_data::Tuple{Float64,Array{Float64}} = percolation_differences(sort(tmp_perc_states),n)
    #betweenness::Array{Float64} = zeros(Float64,n)
    start_time::Float64 = time()
    emp_wimpy_node::Float64 = 0.0
    min_inv_wimpy_node::Float64 = 0.0
    max_perc::Float64 = 0.0
    max_bc::Float64 = 0.0
    max_wv::Float64 = 0.0
    # Partitions 
    number_of_non_empty_partitions::Int64 = 0
    non_empty_partitions::Dict{Int64,Int64} = Dict{Int64,Int64}()
    partitions_ids_map::Dict{Int64,Int64} = Dict{Int64,Int64}()
    partition_index::Array{Int64} = zeros(n)
    part_idx::Int64 = 1
    # Diameter approximation using Propagate/RandomBFS
    # Using Random BFS
    @info("Approximating diameter using Random BFS algorithm")
    flush(stdout)
    diam,time_diam = random_bfs(g,sample_size_diam)
    finish_diam::String = string(round(time_diam; digits=4))
    @info("Estimated diameter "*string(diam)*" in "*finish_diam*" seconds")
    flush(stdout)
    # Sample size related
    omega::Float64 = 1000
    max_num_samples::Float64 = 0.0
    tau::Int64 = trunc(Int64,max(1. / epsilon * (log(1. / delta)) , 100.))
    tau =  trunc(Int64,max(tau,2*(diam -1) * (log(1. / delta))) )
    start_time_bootstrap::Float64 = time()
    @info("Bootstrap phase "*string(tau)*" iterations")
    flush(stdout)
    new_diam_estimate::Array{Array{Int64}} = [[diam] for _ in 1:ntasks]
    final_new_diam_estimate::Int64 = diam
    task_size = cld(tau, ntasks)
    vs_active = [i for i in 1:tau]
    @sync for (t, task_range) in enumerate(Iterators.partition(1:tau, task_size))
        Threads.@spawn for _ in @view(vs_active[task_range])
            _parallel_random_path!(sg,n,percolation_centrality[t],wimpy_variance[t],percolation_states,percolation_data,shortest_path_length[t],percolated_path_length[t],mcrade[t],mc_trials,alpha_sampling,new_diam_estimate[t],true)
        end
    end
  
    # Reducing
    changed::Bool = false
    old_diam::Int64 = final_new_diam_estimate[1]
    for u in 1:n
        for t in 1:ntasks
            final_mcrade[u] += mcrade[t][u]
            final_percolation_centrality[u] += percolation_centrality[t][u]
            final_wimpy_variance[u] += wimpy_variance[t][u]
            final_shortest_path_length[u] += shortest_path_length[t][u]
            final_percolated_path_length[u] += percolated_path_length[t][u]
            if final_new_diam_estimate[1] < new_diam_estimate[t][1]
                changed = true
                final_new_diam_estimate[1] =  new_diam_estimate[t][1]
            end
        end
        #final_mcrade[u] = final_mcrade[u] /tau
        final_percolation_centrality[u] = final_percolation_centrality[u] /tau
        #final_wimpy_variance[u] = final_wimpy_variance[u] /tau
        #final_shortest_path_length[u] = final_shortest_path_length[u] /tau
        #final_percolated_path_length[u] = final_percolated_path_length[u] /tau
    end
    

    #betweenness = betweenness .* [1/tau]
    @info("Empirical peeling phase:")
    flush(stdout)
    if changed
        @info("Diameter estimation refined from "*string(old_diam)*" to "*string(final_new_diam_estimate[1]))
        diam = final_new_diam_estimate[1]
        flush(stdout)
    end
    for i in 1:n
        max_perc = max(max_perc,final_percolation_centrality[i])
        #max_bc = max(max_bc,betweenness[i])

        max_wv = max(max_wv,final_wimpy_variance[i])
        emp_w_node = final_wimpy_variance[i] * 1. /tau
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
    avg_diam_ub::Float64 = upper_bound_average_diameter(delta/8,diam,final_shortest_path_length,tau,true)
    top1bc_upper_bound::Float64 = upper_bound_top_1_bc(max_perc,delta/8,tau)
    wimpy_var_upper_bound::Float64 = upper_bound_top_1_bc(max_wv/tau,delta/8,tau)
    @info("Average shortest path length (upper bound) "*string(avg_diam_ub))
    flush(stdout)
    max_num_samples = upper_bound_samples(top1bc_upper_bound,wimpy_var_upper_bound,avg_diam_ub,epsilon,delta/2 ,false)
    omega = 0.5/epsilon/epsilon * (log2(diam-1)+1+log(2/delta))
    @info("Maximum number of samples "*string(trunc(Int64,floor(max_num_samples)))*" VC Bound "*string(trunc(Int64,floor(omega))))
    @info("Sup perc. centr. estimation "*string(max_perc))
    #println("Sup bc. estimation "*string(max_perc))

    @info("Sup empirical wimpy variance "*string(max_wv/tau))
    flush(stdout)
    omega = 0.5/epsilon/epsilon * (log2(diam-1)+1+log(2/delta))
    if max_num_samples > 0
        omega = max_num_samples
    end
    max_num_samples = omega
    first_stopping_samples::Float64 = 0.0
    eps_guess::Float64 = 1.0
    first_sample_lower::Float64 = 1/epsilon *log(2/delta)
    first_sample_upper::Float64 = omega
    sup_emp_wimpy_var_norm::Float64  = max_wv/tau +1/tau
    finish_bootstrap = string(round(time() - start_time_bootstrap; digits=4))
    @info("Bootstrap completed in "*finish_bootstrap)
    @info("Inferring initial sample size for the geometric sampler")
    flush(stdout)
    num_samples_bs::Float64 = 0.0
    while first_sample_upper - first_sample_lower> 10
        num_samples_bs = (first_sample_upper+first_sample_lower)÷2
        eps_guess = sqrt(2*sup_emp_wimpy_var_norm*log(2/delta) /num_samples_bs) + log(2/delta)/num_samples_bs/3
        if eps_guess > epsilon
            first_sample_lower = num_samples_bs
        else
            first_sample_upper = num_samples_bs
        end
    end
    first_stopping_samples = trunc(Int64,floor(num_samples_bs))
    last_stopping_samples = trunc(Int64,floor(omega))
    @info("Maximum number of iterations "*string(last_stopping_samples))
    @info("Initial sample size "*string(first_stopping_samples))
    if first_stopping_samples >= last_stopping_samples/4
        first_stopping_samples = trunc(Int64,floor(last_stopping_samples/4))
        @info("Initial sample size dropped to "*string(first_stopping_samples))
    end
    flush(stdout)
    # Freeing all the arrays
    percolation_centrality = Array{Array{Float64}}([zeros(Float64,n) for _ in 1:ntasks])
    final_percolation_centrality = Array{Float64}(zeros(Float64,n))
    wimpy_variance = Array{Array{Float64}}([zeros(Float64,n) for _ in 1:ntasks])
    final_wimpy_variance = Array{Float64}(zeros(Float64,n))
    shortest_path_length = Array{Array{Int64}}([zeros(Int64,n+1) for _ in 1:ntasks])
    final_shortest_path_length = Array{Int64}(zeros(Int64,n+1))
    percolated_path_length = Array{Array{Float64}}([zeros(Float64,n+1) for _ in 1:ntasks])
    final_percolated_path_length = Array{Float64}(zeros(Float64,n+1))
    mcrade = Array{Array{Float64}}([zeros(Float64,(n+1)*mc_trials) for _ in 1:ntasks])
    final_mcrade = Array{Float64}(zeros(Float64,(n+1)*mc_trials))

    next_stopping_samples::Float64 = first_stopping_samples
    prev_stopping_samples::Float64 = 0.0
    has_to_stop::Bool = false
    sample_i::Int64 = 0
    num_samples::Int64 = 0
    iteration_index::Int64 =1 
    while !has_to_stop
        sample_i = trunc(Int,next_stopping_samples-prev_stopping_samples)
        task_size = cld(sample_i, ntasks)
        vs_active = [i for i in 1:sample_i]
        @sync for (t, task_range) in enumerate(Iterators.partition(1:sample_i, task_size))
            Threads.@spawn for _ in @view(vs_active[task_range])
                _parallel_random_path!(sg,n,percolation_centrality[t],wimpy_variance[t],percolation_states,percolation_data,shortest_path_length[t],percolated_path_length[t],mcrade[t],mc_trials,alpha_sampling,new_diam_estimate[t],true)
                if (Sys.free_memory() / Sys.total_memory() < 0.1)
                    clean_gc()
                    sleep(0.01)
                end
            end
        end

        num_samples += sample_i
        if num_samples >= omega
            @info("----------------------------| Iteration : "*string(iteration_index)*" |----------------------------")
            has_to_stop = true
            finish_partial = string(round(time() - start_time; digits=4))
            @info("Completed, sampled "*string(num_samples)*"/"*string(omega)* " couples in "*finish_partial)
            @info("---------------------------------------------------------------------------------------------------")
            flush(stdout)
        end
        if !has_to_stop & (num_samples < last_stopping_samples)&(num_samples >= next_stopping_samples)
            changed = false
            old_diam = diam
            for u in 1:n
                for t in 1:ntasks
                    final_mcrade[u] += mcrade[t][u]
                    final_percolation_centrality[u] += percolation_centrality[t][u]
                    final_wimpy_variance[u] += wimpy_variance[t][u]
                    final_shortest_path_length[u] += shortest_path_length[t][u]
                    final_percolated_path_length[u] += percolated_path_length[t][u]
                    if final_new_diam_estimate[1] < new_diam_estimate[t][1]
                        changed = true
                        final_new_diam_estimate[1] =  new_diam_estimate[t][1]
                    end
                end
                #final_mcrade[u] = final_mcrade[u] /tau
                final_percolation_centrality[u] = final_percolation_centrality[u] /tau
                #final_wimpy_variance[u] = final_wimpy_variance[u] /tau
                #final_shortest_path_length[u] = final_shortest_path_length[u] /tau
                #final_percolated_path_length[u] = final_percolated_path_length[u] /tau
            end
            @info("----------------------------| Iteration : "*string(iteration_index)*" |----------------------------")
            if changed
                @info("Diameter estimation refined from "*string(old_diam)*" to "*string(final_new_diam_estimate[1]))
                diam = final_new_diam_estimate[1]
                flush(stdout)
            end
            tmp_omega = Array{Float64}([omega])
            tmp_has_to_stop = Array{Bool}([false])

            _check_stopping_condition!(final_percolation_centrality,final_wimpy_variance,last_stopping_samples,num_samples,epsilon,delta,iteration_index,true,diam,final_shortest_path_length,num_samples,mc_trials,partition_index,partitions_ids_map,final_mcrade,number_of_non_empty_partitions,tmp_omega,tmp_has_to_stop)
            omega = tmp_omega[1]
            has_to_stop = tmp_has_to_stop[1]
            if has_to_stop
                @info("Progressive sampler converged!")
                flush(stdout)
            else
                prev_stopping_samples = next_stopping_samples
                next_stopping_samples,iteration_index = get_next_stopping_sample(next_stopping_samples,iteration_index )
                @info("Increasing sample size to "*string(next_stopping_samples))
                flush(stdout)
            end
            @info("---------------------------------------------------------------------------------------------------")
            flush(stdout)
        end
    end
    final_percolation_centrality = reduce(+,percolation_centrality)

    @info("Estimation completed "*string(round(time() - start_time; digits=4)))
    flush(stdout)
    return final_percolation_centrality .*[1/num_samples],num_samples,max_num_samples,time()-start_time
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
    #println(betweenness)
    #println("---------------")
    #println(percolation_centrality)
    #println(wimpy_variance)
    #println("-------------------------")
    # println(shortest_path_length)
    #println(percolated_path_length)

end