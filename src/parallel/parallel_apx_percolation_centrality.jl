
function parallel_estimate_percolation_centrality(g,percolation_states::Array{Float64},epsilon::Float64,delta::Float64, alpha_sampling::Float64 = 0.1, mc_trials::Int64 = 25, empirical_peeling_a::Float64 = 2.0,sample_size_diam::Int64 = 256 )
    @assert nv(g) == lastindex(percolation_states) "Length of the percolation state array must be the same as the number of nodes"
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
    run_perc::Bool = true
    if std(percolation_states) == 0
        run_perc = false
        @info("Percolation states are all the same, percolation centrality is equal to the betweenness centrality")
    end
    percolation_centrality::Array{Array{Float64}} = [zeros(Float64,n) for _ in 1:ntasks]
    final_percolation_centrality::Array{Float64} = zeros(Float64,n)
    wimpy_variance::Array{Array{Float64}} = [zeros(Float64,n) for _ in 1:ntasks]
    final_wimpy_variance::Array{Float64} = zeros(Float64,n)
    shortest_path_length::Array{Array{Int64}} = [zeros(Int64,n+1) for _ in 1:ntasks]
    final_shortest_path_length::Array{Int64} = zeros(Int64,n+1)
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
    diam,time_diam = parallel_random_bfs(g,sample_size_diam)
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
            _parallel_random_path!(sg,n,percolation_centrality[t],wimpy_variance[t],percolation_states,percolation_data,shortest_path_length[t],mcrade[t],mc_trials,alpha_sampling,new_diam_estimate[t],run_perc,true)
        end
    end
  
    # Reducing
    changed::Bool = false
    old_diam::Int64 = final_new_diam_estimate[1]
    for t in 1:ntasks
        if final_new_diam_estimate[1] < new_diam_estimate[t][1]
            changed = true
            final_new_diam_estimate[1] =  new_diam_estimate[t][1]
        end
    end

    # reducing
    task_size = cld(n, ntasks)
    vs_active = [i for i in 1:n]
    @sync for (t, task_range) in enumerate(Iterators.partition(1:n, task_size))
        Threads.@spawn for u in @view(vs_active[task_range])
            _reduce_data!(u,ntasks,percolation_centrality,wimpy_variance,shortest_path_length,final_percolation_centrality,final_wimpy_variance,final_shortest_path_length)
        end
    end
    #final_percolation_centrality = reduce(+,percolation_centrality) .* [1/tau]
    #final_wimpy_variance = reduce(+,wimpy_variance)
    #final_shortest_path_length = reduce(+,shortest_path_length)
    #=
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
    =#

    #betweenness = betweenness .* [1/tau]
    @info("Empirical peeling phase:")
    flush(stdout)
    if changed
        @info("Diameter estimation refined from "*string(old_diam)*" to "*string(final_new_diam_estimate[1]))
        diam = final_new_diam_estimate[1]
        flush(stdout)
    end
    for i in 1:n
        max_perc = max(max_perc,final_percolation_centrality[i]/tau)
        #max_bc = max(max_bc,betweenness[i])

        max_wv = max(max_wv,final_wimpy_variance[i]/tau)
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
                _parallel_random_path!(sg,n,percolation_centrality[t],wimpy_variance[t],percolation_states,percolation_data,shortest_path_length[t],mcrade[t],mc_trials,alpha_sampling,new_diam_estimate[t],run_perc,false)
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
            for t in 1:ntasks
                if final_new_diam_estimate[1] < new_diam_estimate[t][1]
                    changed = true
                    final_new_diam_estimate[1] =  new_diam_estimate[t][1]
                end
            end
            # reducing
            #r_time = time()
            task_size = cld(n, ntasks)
            vs_active = [i for i in 1:n]
            @sync for (t, task_range) in enumerate(Iterators.partition(1:n, task_size))
                Threads.@spawn for u in @view(vs_active[task_range])
                    _reduce_data!(u,ntasks,percolation_centrality,wimpy_variance,shortest_path_length,final_percolation_centrality,final_wimpy_variance,final_shortest_path_length)
                end
            end
            #println("Iteration "*string(iteration_index)*" reduce time "*string(time()-r_time))
            #final_percolation_centrality = reduce(+,percolation_centrality)
            #final_wimpy_variance = reduce(+,wimpy_variance)
            #final_shortest_path_length = reduce(+,shortest_path_length)
            #=
            for u in 1:n
                for t in 1:ntasks
                    final_mcrade[u] += mcrade[t][u]
                    final_percolation_centrality[u] += percolation_centrality[t][u]
                    final_wimpy_variance[u] += wimpy_variance[t][u]
                    final_shortest_path_length[u] += shortest_path_length[t][u]
                    #final_percolated_path_length[u] += percolated_path_length[t][u]
                    if final_new_diam_estimate[1] < new_diam_estimate[t][1]
                        changed = true
                        final_new_diam_estimate[1] =  new_diam_estimate[t][1]
                    end
                end
                =#
                #final_mcrade[u] = final_mcrade[u] /tau
                #final_percolation_centrality[u] = final_percolation_centrality[u] /num_samples
                #final_wimpy_variance[u] = final_wimpy_variance[u] /tau
                #final_shortest_path_length[u] = final_shortest_path_length[u] /tau
                #final_percolated_path_length[u] = final_percolated_path_length[u] /tau
            #end
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

end



function parallel_estimate_percolation_centrality_lock(g,percolation_states::Array{Float64},epsilon::Float64,delta::Float64, alpha_sampling::Float64 = 0.1, mc_trials::Int64 = 25, empirical_peeling_a::Float64 = 2.0,sample_size_diam::Int64 = 256 )
    @assert nv(g) == lastindex(percolation_states) "Length of the percolation state array must be the same as the number of nodes"
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
    run_perc::Bool = true
    if std(percolation_states) == 0
        run_perc = false
        @info("Percolation states are all the same, percolation centrality is equal to the betweenness centrality")
    end
    #percolation_centrality::Array{Array{Float64}} = [zeros(Float64,n) for _ in 1:ntasks]
    final_percolation_centrality::Array{Float64} = zeros(Float64,n)
    #wimpy_variance::Array{Array{Float64}} = [zeros(Float64,n) for _ in 1:ntasks]
    final_wimpy_variance::Array{Float64} = zeros(Float64,n)
    #shortest_path_length::Array{Array{Int64}} = [zeros(Int64,n+1) for _ in 1:ntasks]
    final_shortest_path_length::Array{Int64} = zeros(Int64,n+1)
    #mcrade::Array{Array{Float64}} = [zeros(Float64,(n+1)*mc_trials) for _ in 1:ntasks]
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
    diam,time_diam = parallel_random_bfs(g,sample_size_diam)
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
    #new_diam_estimate::Array{Array{Int64}} = [[diam] for _ in 1:ntasks]
    final_new_diam_estimate::Array{Int64} = [diam]
    lk::ReentrantLock = ReentrantLock()
    task_size = cld(tau, ntasks)
    vs_active = [i for i in 1:tau]
    @sync for (t, task_range) in enumerate(Iterators.partition(1:tau, task_size))
        Threads.@spawn for _ in @view(vs_active[task_range])
            _parallel_random_path_lk!(sg,n,final_percolation_centrality,final_wimpy_variance,percolation_states,percolation_data,final_shortest_path_length,final_mcrade,mc_trials,alpha_sampling,final_new_diam_estimate,lk,run_perc,true)
        end
    end
  
    # Reducing
    changed::Bool = false
    old_diam::Int64 = diam
    if diam < final_new_diam_estimate[1]
        changed = true
    end
    #=
    for t in 1:ntasks
        if final_new_diam_estimate[1] < new_diam_estimate[t][1]
            changed = true
            final_new_diam_estimate[1] =  new_diam_estimate[t][1]
        end
    end
   

    # reducing
    task_size = cld(n, ntasks)
    vs_active = [i for i in 1:n]
    @sync for (t, task_range) in enumerate(Iterators.partition(1:n, task_size))
        Threads.@spawn for u in @view(vs_active[task_range])
            _reduce_data!(u,ntasks,percolation_centrality,wimpy_variance,shortest_path_length,final_percolation_centrality,final_wimpy_variance,final_shortest_path_length)
        end
    end
     =#
    #final_percolation_centrality = reduce(+,percolation_centrality) .* [1/tau]
    #final_wimpy_variance = reduce(+,wimpy_variance)
    #final_shortest_path_length = reduce(+,shortest_path_length)
    #=
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
    =#

    #betweenness = betweenness .* [1/tau]
    @info("Empirical peeling phase:")
    flush(stdout)
    if changed
        @info("Diameter estimation refined from "*string(old_diam)*" to "*string(final_new_diam_estimate[1]))
        diam = final_new_diam_estimate[1]
        flush(stdout)
    end
    for i in 1:n
        max_perc = max(max_perc,final_percolation_centrality[i]/tau)
        #max_bc = max(max_bc,betweenness[i])

        max_wv = max(max_wv,final_wimpy_variance[i]/tau)
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
    #percolation_centrality = Array{Array{Float64}}([zeros(Float64,n) for _ in 1:ntasks])
    final_percolation_centrality = Array{Float64}(zeros(Float64,n))
    #wimpy_variance = Array{Array{Float64}}([zeros(Float64,n) for _ in 1:ntasks])
    final_wimpy_variance = Array{Float64}(zeros(Float64,n))
    #shortest_path_length = Array{Array{Int64}}([zeros(Int64,n+1) for _ in 1:ntasks])
    final_shortest_path_length = Array{Int64}(zeros(Int64,n+1))
    #mcrade = Array{Array{Float64}}([zeros(Float64,(n+1)*mc_trials) for _ in 1:ntasks])
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
                _parallel_random_path_lk!(sg,n,final_percolation_centrality,final_wimpy_variance,percolation_states,percolation_data,final_shortest_path_length,final_mcrade,mc_trials,alpha_sampling,final_new_diam_estimate,lk,run_perc,false)
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
            if old_diam < final_new_diam_estimate[1]
                changed = true
            end
            #=
            for t in 1:ntasks
                if final_new_diam_estimate[1] < new_diam_estimate[t][1]
                    changed = true
                    final_new_diam_estimate[1] =  new_diam_estimate[t][1]
                end
            end
            
            # reducing
            #r_time = time()
            task_size = cld(n, ntasks)
            vs_active = [i for i in 1:n]
            @sync for (t, task_range) in enumerate(Iterators.partition(1:n, task_size))
                Threads.@spawn for u in @view(vs_active[task_range])
                    _reduce_data!(u,ntasks,percolation_centrality,wimpy_variance,shortest_path_length,final_percolation_centrality,final_wimpy_variance,final_shortest_path_length)
                end
            end
            =#
            #println("Iteration "*string(iteration_index)*" reduce time "*string(time()-r_time))
            #final_percolation_centrality = reduce(+,percolation_centrality)
            #final_wimpy_variance = reduce(+,wimpy_variance)
            #final_shortest_path_length = reduce(+,shortest_path_length)
            #=
            for u in 1:n
                for t in 1:ntasks
                    final_mcrade[u] += mcrade[t][u]
                    final_percolation_centrality[u] += percolation_centrality[t][u]
                    final_wimpy_variance[u] += wimpy_variance[t][u]
                    final_shortest_path_length[u] += shortest_path_length[t][u]
                    #final_percolated_path_length[u] += percolated_path_length[t][u]
                    if final_new_diam_estimate[1] < new_diam_estimate[t][1]
                        changed = true
                        final_new_diam_estimate[1] =  new_diam_estimate[t][1]
                    end
                end
                =#
                #final_mcrade[u] = final_mcrade[u] /tau
                #final_percolation_centrality[u] = final_percolation_centrality[u] /num_samples
                #final_wimpy_variance[u] = final_wimpy_variance[u] /tau
                #final_shortest_path_length[u] = final_shortest_path_length[u] /tau
                #final_percolated_path_length[u] = final_percolated_path_length[u] /tau
            #end
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
    #final_percolation_centrality = reduce(+,percolation_centrality)

    @info("Estimation completed "*string(round(time() - start_time; digits=4)))
    flush(stdout)
    return final_percolation_centrality .*[1/num_samples],num_samples,max_num_samples,time()-start_time

end


function parallel_estimate_percolation_centrality_era(g,percolation_states::Array{Float64},epsilon::Float64,delta::Float64,initial_sample::Int64 = 0,geo::Float64 = 1.2,sample_size_diam::Int64 = 256 )
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
    if initial_sample == 0
        @info("Inferring the size of the first sample in the schedule")
        initial_sample = trunc(Int64,floor((1+8*epsilon + sqrt(1+16*epsilon)*log(6/delta))/(4*epsilon*epsilon)))
        @info("The size of the first sample in the schedule is "*string(initial_sample))
        flush(stdout)
    end

    local_B::Array{Dict{Float64,Float64}} =  [Dict{Float64,Float64}() for i in 1:ntasks]
    local_B_1::Array{Array{Float64}} =  [zeros(n) for i in 1:ntasks]
    final_B_1::Array{Float64} = zeros(Float64,n)
    local_B_2::Array{Array{Float64}} =  [zeros(n) for i in 1:ntasks]
    B_vectorized::Array{Float64} =Array{Float64}([])
    final_B_2::Array{Float64} = zeros(Float64,n)
    start_time::Float64 = time()
    k::Int64 = 0
    j::Int64 = 2
    keep_sampling::Bool = true
    sampled_so_far::Int64 = 0
    new_sample::Int64 = 0
    sample_size_schedule::Array{Int64} = [0,initial_sample]
    xi::Float64 = 0
    tmp_perc_states::Array{Float64} = copy(percolation_states)
    percolation_data::Tuple{Float64,Array{Float64}} = percolation_differences(sort(tmp_perc_states),n)
    @info("Approximating diameter using Random BFS algorithm")
    flush(stdout)
    diam,time_diam = parallel_random_bfs(g,sample_size_diam)
    finish_diam::String = string(round(time_diam; digits=4))
    @info("Estimated diameter "*string(diam)*" in "*finish_diam*" seconds")
    flush(stdout)
    max_sample::Int64 = trunc(Int64,0.5/epsilon/epsilon * (log2(diam-1)+1+log(2/delta)))
    if diam == 0
        max_sample = Inf
    end
    while keep_sampling
        k+=1
        if (k >= 2)
            new_sample = trunc(Int,geo^k*sample_size_schedule[2])
            push!(sample_size_schedule,new_sample)
        end

        sample_i = sample_size_schedule[j]-sample_size_schedule[j-1]
        task_size = cld(sample_i, ntasks)
        vs_active = [i for i in 1:sample_i]
        @sync for (t, task_range) in enumerate(Iterators.partition(1:sample_i, task_size))
            Threads.@spawn for _ in @view(vs_active[task_range])
                _parallel_sz_bfs!(g,percolation_states,percolation_data,local_B[t],local_B_1[t],local_B_2[t])
            end
        end
        sampled_so_far += sample_i
        B_vectorized = reduce_dictionary(local_B)
        if length(B_vectorized)>0
            omega = compute_xi(B_vectorized,sample_size_schedule[j])
            delta_i = delta/2^k
            xi = 2*omega + (log(3/delta_i)+sqrt((log(3/delta_i)+4*sample_size_schedule[j]*omega)*log(3/delta_i)))/sample_size_schedule[j]  + sqrt(log(3/delta_i)/(2*sample_size_schedule[j]))
        else
            xi = Inf
        end
        @info("ERA Upperbound "*string(xi)*" Target SD "*string(epsilon)*" #Sampled pairs "*string(sample_size_schedule[j]))
        if xi <= epsilon || sample_size_schedule[j] >= max_sample
            keep_sampling = false
        else
            j+=1
        end

    end
    final_B_1 = reduce(+,local_B_1) .*[1/sample_size_schedule[end]]
    finish_time::Float64 = time()-start_time
    @info("Converged! Sampled "*string(sample_size_schedule[end])*"/"*string(max_sample)*" couples in "*string(round(finish_time;digits = 4))*" seconds ")
    return final_B_1,sample_size_schedule,max_sample,xi,finish_time
end


function _parallel_sz_bfs!(g,percolation_states::Array{Float64},percolation_data::Tuple{Float64,Array{Float64}},B::Dict{Float64,Float64},B_1::Array{Float64},B_2::Array{Float64})
    n::Int64 = nv(g)
    q::Queue{Int64} = Queue{Int64}()
    ball::Array{Int16} = zeros(Int16,n)
    n_paths::Array{Int64} = zeros(Int64,n)
    dist::Array{Int64} = zeros(Int64,n)
    pred::Array{Array{Int64}} = [Array{Int64}([]) for _ in 1:n]

    s::Int64 = sample(1:n)
    z::Int64 = sample(1:n)
    w::Int64 = 0
    d_z_min::Float64 = Inf
    q_backtrack::Queue{Int64} = Queue{Int64}()
    while (s == z)
        z = sample(1:n)
    end
    enqueue!(q,s)
    dist[s] = 0
    n_paths[s] = 1
    ball[s] = 1
    while length(q) != 0
        w = dequeue!(q)
        if dist[w] < d_z_min
            for v in outneighbors(g,w)
                if (ball[v] == 0)
                    dist[v] = dist[w] +1
                    ball[v] = 1
                    n_paths[v] = n_paths[w]
                    if (v == z)
                        if dist[v] < d_z_min
                            d_z_min = dist[v]
                        end
                        if length(q_backtrack) == 0
                            enqueue!(q_backtrack,z)
                        end
                    end
                    push!(pred[v],w)
                    enqueue!(q,v)
                elseif (dist[v] == dist[w] + 1)
                    n_paths[v] += n_paths[w]
                    push!(pred[v],w)
                end

            end
        end
    end
    #backtrack
   while length(q_backtrack)!=0
        w = dequeue!(q_backtrack)
        if w != s && w != z
            summand = (n_paths[w]/n_paths[z]) *(ramp(percolation_states[s],percolation_states[z])/percolation_data[2][w])  
            # Updating phase
            b = B_2[w]
            b_1 = b + summand^2
            if !haskey(B,b_1) 
                B[b_1] = 1
            else
                B[b_1] += 1
            end
            if b > 0 && B[b] >= 1
                B[b] -= 1
            end
            if b > 0 && B[b] == 0
                delete!(B, b)
            end
            B_1[w] += summand
            B_2[w] += summand^2
        end
        for p in pred[w]
            enqueue!(q_backtrack,p)
        end
   end
   

   return nothing
end


function reduce_dictionary(B::Array{Dict{Float64,Float64}})::Vector{Float64}
    B_new::Dict{Float64,Float64} = Dict{Float64,Float64}()
    for cur_B in B
        for b_1 in keys(cur_B)
            if !haskey(B_new,b_1)
                B_new[b_1] = 1
            else
                B_new[b_1] +=1
            end
        end
    end
    return collect(keys(B_new))
end


# Fixed sample size
function parallel_estimate_percolation_centrality_fixed_sample_size(g,percolation_states::Array{Float64},epsilon::Float64,delta::Float64,sample_size_diam::Int64 = 256 )
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
    start_time::Float64 = time()
    
    

    final_percolation_centrality::Array{Float64} = zeros(Float64,n)
    percolation_centrality::Array{Array{Float64}} =  [zeros(Float64,n) for _ in 1:ntasks]
    
    #sg::static_graph = static_graph(adjacency_list(g),incidency_list(g))
    tmp_perc_states::Array{Float64} = copy(percolation_states)
    percolation_data::Tuple{Float64,Array{Float64}} = percolation_differences(sort(tmp_perc_states),n)
    @info("Approximating diameter using Random BFS algorithm")
    flush(stdout)
    diam,time_diam = parallel_random_bfs(g,sample_size_diam)
    finish_diam::String = string(round(time_diam; digits=4))
    @info("Estimated diameter "*string(diam)*" in "*finish_diam*" seconds")
    flush(stdout)
    max_sample::Int64 = trunc(Int64,0.5/epsilon/epsilon * (log2(diam-1)+1+log(2/delta)))
    @info("Maximum sample size "*string(max_sample))
    flush(stdout)
    task_size = cld(max_sample, ntasks)
    vs_active = [i for i in 1:max_sample]
    @sync for (t, task_range) in enumerate(Iterators.partition(1:max_sample, task_size))
        Threads.@spawn for _ in @view(vs_active[task_range])
            _parallel_sz_bfs_fss!(g,percolation_states,percolation_data,percolation_centrality[t])
        end
    end

     # reducing
     task_size = cld(n, ntasks)
     vs_active = [i for i in 1:n]
     @sync for (t, task_range) in enumerate(Iterators.partition(1:n, task_size))
         Threads.@spawn for u in @view(vs_active[task_range])
             _reduce_data_a!(u,ntasks,percolation_centrality,final_percolation_centrality)
         end
     end


    finish_time::Float64 = time()-start_time
    @info("Completed! Sampled "*string(max_sample)*" couples in "*string(round(finish_time;digits = 4))*" seconds ")
    flush(stdout)
    return final_percolation_centrality.*[1/max_sample],max_sample,finish_time
end


function _parallel_sz_bfs_fss!(g,percolation_states::Array{Float64},percolation_data::Tuple{Float64,Array{Float64}},percolation_centrality::Array{Float64})
    n::Int64 = nv(g)
    s::Int64 = sample(1:n)
    z::Int64 = sample(1:n)
    w::Int64 = 0
    q::Queue{Int64} = Queue{Int64}()
    ball::Array{Int16} = zeros(Int16,n)
    n_paths::Array{Int64} = zeros(Int64,n)
    dist::Array{Int64} = zeros(Int64,n)
    pred::Array{Array{Int64}} = [Array{Int64}([]) for _ in 1:n]
    d_z_min::Float64 = Inf
    tot_weight::Int64 = 0
    random_edge::Int64 = 0
    cur_edge::Int64 = 0
    path::Array{Int64} = Array{Int64}([])
    q_backtrack::Queue{Int64} = Queue{Int64}()
    while (s == z)
        z = sample(1:n)
    end
   
    enqueue!(q,s)
    dist[s] = 0
    n_paths[s] = 1
    ball[s] = 1
    while length(q) != 0
        w = dequeue!(q)
        if dist[w] < d_z_min
            for v in outneighbors(g,w)
                if (ball[v] == 0)
                    dist[v] = dist[w] +1
                    ball[v] = 1
                    n_paths[v] = n_paths[w]
                    if (v == z)
                        if dist[v] < d_z_min
                            d_z_min = dist[v]
                        end
                        if length(q_backtrack) == 0
                            enqueue!(q_backtrack,z)
                        end
                    end
                    push!(pred[v],w)
                    enqueue!(q,v)
                elseif (dist[v] == dist[w] + 1)
                    n_paths[v] += n_paths[w]
                    push!(pred[v],w)
                end

            end
        end
    end
    #backtrack
    if length(q_backtrack) != 0
        w = dequeue!(q_backtrack)
        tot_weight = n_paths[z]
        random_edge = rand(0:tot_weight-1)
        cur_edge = 0
        for p in pred[z]
            cur_edge += n_paths[p]
            if cur_edge > random_edge
                _backtrack_path!(s,z,p,path,n_paths,pred)
                break;
            end
        end
        for w in path
            percolation_centrality[w] += (ramp(percolation_states[s],percolation_states[z])/percolation_data[2][w])  
        end
    end  

   return nothing
end