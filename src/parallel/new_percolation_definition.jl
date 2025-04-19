
function parallel_estimate_percolation_centrality_new_lock(g,percolation_states::Array{Float64},epsilon::Float64,delta::Float64, alpha_sampling::Float64 = 0.1, mc_trials::Int64 = 25, empirical_peeling_a::Float64 = 2.0,sample_size_diam::Int64 = 256 )
    @assert nv(g) == lastindex(percolation_states) "Length of the percolation state array must be the same as the number of nodes"
    n::Int64 = nv(g)
    m::Int64 = ne(g)
    directed::Bool = is_directed(g)
    dv,_,_,_=compute_d_max(nv(g) ,percolation_states)
    max_d_v::Float64 = maximum(dv)
    @info("----------------------------------------| Stats |--------------------------------------------------")
    @info("Analyzing graph")
    @info("Number of nodes "*string(n))
    @info("Number of edges "*string(m))
    @info("Directed ? "*string(directed))
    @info("Maximum Percolated state "*string(maximum(percolation_states)))
    @info("Minimum Percolated state "*string(minimum(percolation_states)))
    @info("Average Percolated state "*string(mean(percolation_states))*" std "*string(std(percolation_states)))
    @info("Epsilon "*string(epsilon)*" Delta "*string(delta)*" Alpha "*string(alpha_sampling)*" MC Trials "*string(mc_trials)*" Empirical peeling param. "*string(empirical_peeling_a))
    @info("Maximum d_v "*string(max_d_v))
    @info("Using "*string(nthreads())* " Threads")
    @info("---------------------------------------------------------------------------------------------------")
    flush(stderr)
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
    tmp_perc_states::Dict{Int64,Float64} = Dict(v => percolation_states[v] for v in 1:n)
    sorted_dict = OrderedDict(sort(collect(tmp_perc_states), by = kv -> kv[2]))
    percolation_data::Tuple{Float64,Dict{Int64,Float64}} = percolation_differences(sorted_dict,n)
    #tmp_perc_states::Array{Float64} = copy(percolation_states)
    #percolation_data::Tuple{Float64,Array{Float64}} = percolation_differences(sort(tmp_perc_states),n)
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
    flush(stderr)
    diam,time_diam = parallel_random_bfs(g,sample_size_diam)
    finish_diam::String = string(round(time_diam; digits=4))
    @info("Estimated diameter "*string(diam)*" in "*finish_diam*" seconds")
    flush(stderr)
    # Sample size related
    omega::Float64 = 1000
    max_num_samples::Float64 = 0.0
    tau::Int64 = trunc(Int64,max(1. / epsilon * (log(1. / delta)) , 100.))
    tau =  trunc(Int64,max(tau,2*(diam -1) * (log(1. / delta))) )
    start_time_bootstrap::Float64 = time()
    @info("Bootstrap phase "*string(tau)*" iterations")
    flush(stderr)
    #new_diam_estimate::Array{Array{Int64}} = [[diam] for _ in 1:ntasks]
    final_new_diam_estimate::Array{Int64} = [diam]
    lk::ReentrantLock = ReentrantLock()
    task_size = cld(tau, ntasks)
    vs_active = [i for i in 1:tau]
    @sync for (t, task_range) in enumerate(Iterators.partition(1:tau, task_size))
        Threads.@spawn for _ in @view(vs_active[task_range])
            _parallel_random_path_weighted_lk!(sg,n,final_percolation_centrality,final_wimpy_variance,percolation_states,percolation_data,final_shortest_path_length,final_mcrade,mc_trials,alpha_sampling,final_new_diam_estimate,lk,run_perc,true)
        end
    end
  
    # Reducing
    changed::Bool = false
    old_diam::Int64 = diam
    if diam < final_new_diam_estimate[1]
        changed = true
    end
   

    @info("Empirical peeling phase:")
    flush(stderr)
    if changed
        @info("Diameter estimation refined from "*string(old_diam)*" to "*string(final_new_diam_estimate[1]))
        diam = final_new_diam_estimate[1]
        flush(stderr)
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
    flush(stderr)
    max_num_samples = upper_bound_samples(top1bc_upper_bound,wimpy_var_upper_bound,avg_diam_ub,epsilon,delta/2 ,false)
    omega = 0.5/epsilon/epsilon * (log2(diam-1)+1+log(2/delta))
    @info("Maximum number of samples "*string(trunc(Int64,floor(max_num_samples)))*" VC Bound "*string(trunc(Int64,floor(omega))))
    @info("Sup perc. centr. estimation "*string(max_perc))
    #println("Sup bc. estimation "*string(max_perc))

    @info("Sup empirical wimpy variance "*string(max_wv/tau))
    flush(stderr)
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
    flush(stderr)
    time_estimation::Float64 = time()
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
    if first_stopping_samples == 0
        first_stopping_samples = 10
    end
    last_stopping_samples = trunc(Int64,floor(omega))
    @info("Maximum number of iterations "*string(last_stopping_samples))
    @info("Initial sample size "*string(first_stopping_samples))
    if first_stopping_samples >= last_stopping_samples/4
        first_stopping_samples = trunc(Int64,floor(last_stopping_samples/4))
        @info("Initial sample size dropped to "*string(first_stopping_samples))
    end
    flush(stderr)
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
                _parallel_random_path_weighted_lk!(sg,n,final_percolation_centrality,final_wimpy_variance,percolation_states,percolation_data,final_shortest_path_length,final_mcrade,mc_trials,alpha_sampling,final_new_diam_estimate,lk,run_perc,false)
                #=
                if (Sys.free_memory() / Sys.total_memory() < 0.1)
                    clean_gc()
                    sleep(0.01)
                end
                =#
            end
        end

        num_samples += sample_i
        if num_samples >= omega
            @info("----------------------------| Iteration : "*string(iteration_index)*" |----------------------------")
            has_to_stop = true
            finish_partial = string(round(time() - start_time; digits=4))
            @info("Completed, sampled "*string(num_samples)*"/"*string(omega)* " couples in "*finish_partial)
            @info("---------------------------------------------------------------------------------------------------")
            flush(stderr)
        end
        if !has_to_stop & (num_samples < last_stopping_samples)&(num_samples >= next_stopping_samples)
            changed = false
            old_diam = diam
            if old_diam < final_new_diam_estimate[1]
                changed = true
            end
           
            @info("----------------------------| Iteration : "*string(iteration_index)*" |----------------------------")
            if changed
                @info("Diameter estimation refined from "*string(old_diam)*" to "*string(final_new_diam_estimate[1]))
                diam = final_new_diam_estimate[1]
                flush(stderr)
            end
            tmp_omega = Array{Float64}([omega])
            tmp_has_to_stop = Array{Bool}([false])

            _check_stopping_condition!(final_percolation_centrality,final_wimpy_variance,last_stopping_samples,num_samples,epsilon,delta,iteration_index,true,diam,final_shortest_path_length,num_samples,mc_trials,partition_index,partitions_ids_map,final_mcrade,number_of_non_empty_partitions,tmp_omega,tmp_has_to_stop)
            omega = tmp_omega[1]
            has_to_stop = tmp_has_to_stop[1]
            if has_to_stop
                @info("Progressive sampler converged!")
                flush(stderr)
            else
                prev_stopping_samples = next_stopping_samples
                next_stopping_samples,iteration_index = get_next_stopping_sample(next_stopping_samples,iteration_index )
                @info("Increasing sample size to "*string(next_stopping_samples))
                flush(stderr)
            end
            @info("---------------------------------------------------------------------------------------------------")
            flush(stderr)
        end
    end
    #final_percolation_centrality = reduce(+,percolation_centrality)

    @info("Estimation completed "*string(round(time() - start_time; digits=4)))
    flush(stderr)
    return final_percolation_centrality .*[1/num_samples],num_samples,max_num_samples,time()-start_time,time()-time_estimation

end


function _random_path_non_unif_silv!(sg::static_graph,n, percolation_centrality::Array{Float64},percolation_states::Array{Float64},percolation_data::Tuple{Float64,Dict{Int64,Float64}},alpha_sampling::Float64,lk::ReentrantLock,uniform_sampling::Bool)
    ball_indicator::Array{Int16} = zeros(Int16,n)
    n_paths::Array{UInt128} = zeros(UInt128,n)
    dist::Array{Int64} = zeros(Int64,n)
    q::Array{Int64} = zeros(Int64,n)
    pred::Array{Array{Int64}} = [Array{Int64}([]) for _ in 1:n]
    u::Int64 = sample(1:n)
    v::Int64 = u
    if !uniform_sampling
        u,v = weighted_sample_kappa(percolation_states)
    else
        while (u == v)
            v = sample(1:n)
        end
    end
    if percolation_states[u] > percolation_states[v] 
        end_q = 2
        #q[1] = u
        #q[2] = v
        path_map::Dict{Int64,UInt128} = Dict{Int64,UInt128}()
        num_path_to_sample::UInt128 = 0
        ball_indicator[u] = @visited_s
        ball_indicator[v] = @visited_z
        n_paths[u] = 1
        n_paths[v] = 1
        dist[u] = 0
        dist[v] = 0

        sp_edges = Vector{Tuple{Int64, Int64}}()
        have_to_stop = false
        start_u, end_u = 1, 2
        start_v, end_v = 2, 3
        vis_edges = 0

        sum_degs_u = 0
        sum_degs_v = 0
        degrees = sg.degrees_adj
        q_s::Queue{Int64} = Queue{Int64}()
        q_z::Queue{Int64} = Queue{Int64}()
        enqueue!(q_s,u)
        enqueue!(q_z,v)
        to_expand = @adjacency

        while !have_to_stop
            if sum_degs_u <= sum_degs_v
                start_cur, end_cur = start_u, end_u
                start_cur = 0
                end_cur = length(q_s)
                start_u = end_q
                new_end_cur = Ref(end_u)
                end_u = end_q
                sum_degs_u = 0
                sum_degs_cur = Ref(sum_degs_u)
                degrees = sg.degrees_adj
                adj = sg.adjacency
                to_expand = @adjacency

            else
                start_cur, end_cur = start_v, end_v
                start_cur = 0
                end_cur = length(q_z)
                start_v = end_q
                new_end_cur = Ref(end_v)
                end_v = end_q
                sum_degs_v = 0
                sum_degs_cur = Ref(sum_degs_v)
                degrees = sg.degrees_idj
                adj = sg.incidency
                to_expand = @incidency

            end

            while start_cur < end_cur
                #x = q[start_cur]
                if to_expand == @adjacency
                    x = dequeue!(q_s)
                else
                    x = dequeue!(q_z)
                end
                start_cur += 1
                neigh_num = degrees[x]

                for j in 1:neigh_num
                    vis_edges += 1
                    y = adj[x][j]
                    if ball_indicator[y] == 0
                        sum_degs_cur[] += degrees[y]
                        n_paths[y] = n_paths[x]
                        ball_indicator[y] = ball_indicator[x]
                        if (to_expand == @adjacency)
                            enqueue!(q_s,y)
                        else
                            enqueue!(q_z,y)
                        end
                        #q[end_q] = y
                        end_q += 1
                        new_end_cur[] += 1
                        push!(pred[y],x)
                        dist[y] = dist[x] + 1
                    elseif ball_indicator[y] != ball_indicator[x]
                        have_to_stop = true
                        push!(sp_edges, (x, y))
                    elseif dist[y] == dist[x] + 1
                        n_paths[y] += n_paths[x]
                        push!(pred[y],x)
                    end
                end
            end

            if sum_degs_cur[] == 0
                have_to_stop = true
            end
        end

        if length(sp_edges) != 0
            #println("IT IS NOT 0")
            tot_weight = sum(n_paths[x] * n_paths[y] for (x, y) in sp_edges)
            num_path_to_sample = 1
            if (alpha_sampling > 0 && tot_weight > 1)
                num_path_to_sample = trunc(Int64,floor(alpha_sampling * tot_weight))
            end
            for j in 1:num_path_to_sample
                random_edge = rand(0:tot_weight-1)
                cur_edge = 0
                path::Array{Int64} =Array{Int64}([])

                for (x, y) in sp_edges
                    cur_edge += n_paths[x] * n_paths[y]
                    if cur_edge > random_edge
                        ##println("CUR EDG ",cur_edge, " rnd edg ",random_edge, " x ",x," y ",y)
                        #println("S ",u, " Z ",v)
                        _backtrack_path!(u, v, x, path, n_paths, pred)
                        _backtrack_path!(u, v, y, path, n_paths, pred)
                        break
                    end
                end
                for w in path
                    if haskey(path_map,w)
                        path_map[w] +=1
                    else
                        path_map[w] = 1
                    end
                end
                #for i in 1:end_q
                #    ball_indicator[q[i]] = @unvisited
                #end
                #println(path)
                
            end
            begin
                lock(lk)
                try
                    for w in keys(path_map)
                        percolation_value = path_map[w]/num_path_to_sample
                        
                        if uniform_sampling
                            percolation_value =  (ramp(percolation_states[u],percolation_states[v])/percolation_data[1])  *  path_map[w]/num_path_to_sample
                        end
                        percolation_centrality[w] += percolation_value
                    end
                finally
                    unlock(lk)
                end
            end
        end
    end
    return nothing
end

#=

function _random_path_non_unif!(sg::static_graph,n, percolation_centrality::Array{Float64},percolation_states::Array{Float64},percolation_data::Tuple{Float64,Dict{Int64,Float64}},alpha_sampling::Float64,lk::ReentrantLock,uniform_sampling::Bool)
    
    u::Int64 = sample(1:n)
    v::Int64 = u
    if !uniform_sampling
        u,v = weighted_sample_kappa(percolation_states)
    else
        while (u == v)
            v = sample(1:n)
        end
    end
    if percolation_states[u] > percolation_states[v] 
        # Variables for the visit
        ball_indicator::Array{Int16} = zeros(Int16,n)
        n_paths::Array{UInt128} = zeros(UInt128,n)
        dist::Array{Int64} = zeros(Int64,n)
        pred::Array{Array{Int64}} = [Array{Int64}([]) for _ in 1:n]
        end_q = 2
        #q[1] = u
        #q[2] = v
        tot_weight::UInt128 = 0
        random_edge::UInt128 = 0
        cur_edge::UInt128 = 0
        ball_indicator[u] = @visited_s
        ball_indicator[v] = @visited_z
        n_paths[u] = 1
        n_paths[v] = 1
        dist[u] = 0
        dist[v] = 0

        sp_edges::Vector{Tuple{Int64, Int64}} = Vector{Tuple{Int64, Int64}}()
        have_to_stop::Bool = false
        start_u, end_u = 1, 2
        start_v, end_v = 2, 3
        vis_edges::Int64 = 0

        sum_degs_u::Int64 = 0
        sum_degs_v::Int64 = 0
        degrees = sg.degrees_adj
        q_s::Queue{Int64} = Queue{Int64}()
        q_z::Queue{Int64} = Queue{Int64}()
        enqueue!(q_s,u)
        enqueue!(q_z,v)
        to_expand::Int32 = @adjacency

        while !have_to_stop
            if sum_degs_u <= sum_degs_v
                start_cur, end_cur = start_u, end_u
                start_cur = 0
                end_cur = length(q_s)
                start_u = end_q
                new_end_cur = Ref(end_u)
                end_u = end_q
                sum_degs_u = 0
                sum_degs_cur = Ref(sum_degs_u)
                degrees = sg.degrees_adj
                adj = sg.adjacency
                to_expand = @adjacency

            else
                start_cur, end_cur = start_v, end_v
                start_cur = 0
                end_cur = length(q_z)
                start_v = end_q
                new_end_cur = Ref(end_v)
                end_v = end_q
                sum_degs_v = 0
                sum_degs_cur = Ref(sum_degs_v)
                degrees = sg.degrees_idj
                adj = sg.incidency
                to_expand = @incidency

            end

            while start_cur < end_cur
                #x = q[start_cur]
                if to_expand == @adjacency
                    x = dequeue!(q_s)
                else
                    x = dequeue!(q_z)
                end
                start_cur += 1
                neigh_num = degrees[x]

                for j in 1:neigh_num
                    vis_edges += 1
                    y = adj[x][j]
                    if ball_indicator[y] == 0
                        sum_degs_cur[] += degrees[y]
                        n_paths[y] = n_paths[x]
                        ball_indicator[y] = ball_indicator[x]
                        if (to_expand == @adjacency)
                            enqueue!(q_s,y)
                        else
                            enqueue!(q_z,y)
                        end
                        #q[end_q] = y
                        end_q += 1
                        new_end_cur[] += 1
                        push!(pred[y],x)
                        dist[y] = dist[x] + 1
                    elseif ball_indicator[y] != ball_indicator[x]
                        have_to_stop = true
                        push!(sp_edges, (x, y))
                    elseif dist[y] == dist[x] + 1
                        n_paths[y] += n_paths[x]
                        push!(pred[y],x)
                    end
                end
            end

            if sum_degs_cur[] == 0
                have_to_stop = true
            end
        end

        if length(sp_edges) != 0
            #println("IT IS NOT 0")
            tot_weight = sum(n_paths[x] * n_paths[y] for (x, y) in sp_edges)
            random_edge = rand(0:tot_weight-1)
            cur_edge = 0
            path::Array{Int64} =Array{Int64}([])

            for (x, y) in sp_edges
                cur_edge += n_paths[x] * n_paths[y]
                if cur_edge > random_edge
                    ##println("CUR EDG ",cur_edge, " rnd edg ",random_edge, " x ",x," y ",y)
                    #println("S ",u, " Z ",v)
                    _backtrack_path!(u, v, x, path, n_paths, pred)
                    _backtrack_path!(u, v, y, path, n_paths, pred)
                    break
                end
            end

            #for i in 1:end_q
            #    ball_indicator[q[i]] = @unvisited
            #end
            #println(path)
            #println("Number of path to backtrack is ",trunc(UInt128,floor(alpha_sampling * tot_weight)), " for thread ",(threadid()))
            begin
                lock(lk)
                try
                    for w in path
                        percolation_value = 1
                        
                        if uniform_sampling
                            percolation_value =  (ramp(percolation_states[u],percolation_states[v])/percolation_data[1])  
                        end
                        percolation_centrality[w] += percolation_value
                    end
                finally
                    unlock(lk)
                end
            end
        end
    end
    return nothing
end

=#

function _parallel_random_path_weighted_lk!(sg::static_graph,n::Int64,percolation_centrality::Array{Float64},wimpy_variance::Array{Float64},percolation_states::Array{Float64},percolation_data::Tuple{Float64,Dict{Int64,Float64}} ,shortest_path_length::Array{Int64},mcrade::Array{Float64},mc_trials::Int64,alpha_sampling::Float64,new_diam_estimate::Array{Int64},lk::ReentrantLock,run_perc::Bool = true,boostrap_phase::Bool = false,uniform_sampling::Bool = false)

    #q::Array{Int64} = zeros(Int64,n)
    ball::Array{Int16} = zeros(Int16,n)
    n_paths::Array{UInt128} = zeros(UInt128,n)
    dist::Array{Int64} = zeros(Int64,n)
    pred::Array{Array{Int64}} = [Array{Int64}([]) for _ in 1:n]
    #num_paths::Array{Int64} = zeros(Int64,n)
    end_q::UInt32 = 1
    tot_weight::UInt128 = 0
    cur_edge::UInt128 = 0
    random_edge::UInt128 = 0
    s::Int64 = sample(1:n)
    z::Int64 = s
    if !uniform_sampling
        s,z = weighted_sample_kappa(percolation_states)
    else
        while (s == z)
            z = sample(1:n)
        end
    end
    #println("s = ",s," percolation s ",percolation_states[s]," z = ",z," percolation z ",percolation_states[z])
    x::Int64  = 0
    y::Int64 = 0
    longest_dist::Int64 = 0
    have_to_stop::Bool = false;
    start_s::Int64 = 1
    start_z::Int64 = 2
    end_s::UInt64 = 2
    end_z::UInt64 = 3
    start_cur::UInt64 = 0
    end_cur::UInt64 = 0
    new_end_cur::UInt64 = 0
    sum_degs_s::UInt64 = 0
    sum_degs_z::UInt64 = 0
    sum_degs_cur::UInt64 = 0
    neigh_num::UInt64 = 0
    to_expand::Int16 = 0
    vis_edges::Int64 = 0
    num_path_to_sample::UInt128 = 1
    num_paths::UInt128 = 0
    sp_edges::Array{Tuple{Int64,Int64}} = Array{Tuple{Int64,Int64}}([])
    path::Array{Int64} = Array{Int64}([])
    path_map::Dict{Int64,UInt128} = Dict{Int64,UInt128}()
    pm_value::UInt128 = 0
    percolation_value::Float64 = 0.0
    lambdas::Array{Int64} = zeros(mc_trials)
    maxval_lambdas::Int64 = 100000000
    if percolation_states[s] > percolation_states[z]
        maxval_half::Int64 = maxval_lambdas/2
        if !boostrap_phase
            for j in 1:mc_trials
                lambdas[j] = (sample(0:maxval_lambdas-1) >= maxval_half)*2-1
            end
        end
        q_s::Queue{Int64} = Queue{Int64}()
        q_z::Queue{Int64} = Queue{Int64}()
        enqueue!(q_s,s)
        enqueue!(q_z,z)

        #end_q = 2
        #q[1] = s
        #q[2] = z
        #println("q ",q)
        ball[s] = @visited_s
        ball[z] = @visited_z

        n_paths[s] = 1
        n_paths[z] = 1
        dist[s] = 0
        dist[z] = 0
        while (!have_to_stop)
            
            # decide what ball we have to expand
            if (sum_degs_s <= sum_degs_z)
                #start_cur = start_s
                start_cur = 0
                end_cur = length(q_s)
                #end_cur = end_s
                start_s = end_q
                new_end_cur = end_s
                end_s = end_q
                sum_degs_s = 0
                sum_degs_cur = sum_degs_s 
                to_expand = @adjacency
                #println("SUM DEG S ",sum_degs_s)
            else
                #start_cur = start_z
                start_cur = 0
                end_cur = length(q_z)
                #end_cur = end_z
                start_z = end_q
                new_end_cur = end_z
                end_z = end_q
                sum_degs_z = 0
                sum_degs_cur = sum_degs_z
                to_expand = @incidency
                #println("SUM DEG Z ",sum_degs_z)
            end
            
            while (start_cur < end_cur)
            # println("START CUR ",start_cur," END CUR ",end_cur)
                # to test
                #x = q[start_cur]
                if to_expand == @adjacency
                    x = dequeue!(q_s)
                else
                    x = dequeue!(q_z)
                end
                start_cur += 1
                if (to_expand == @adjacency)
                    neigh_num = lastindex(sg.adjacency[x])
                    #println("Expanding forward from ",x)
                else
                    neigh_num = lastindex(sg.incidency[x])
                    #println("Expanding backward from ",x)
                end
                for j in 1:neigh_num
                    #println("now ",q)
                    vis_edges +=1
                    if (to_expand == @adjacency)
                        y = sg.adjacency[x][j]
                    else
                        y = sg.incidency[x][j]
                    end
                    #println(" y ",y)
                    if (ball[y] == @unvisited)
                        if (to_expand == @adjacency)
                            sum_degs_cur += lastindex(sg.adjacency[y])
                            enqueue!(q_s,y)
                        else
                            sum_degs_cur += lastindex(sg.incidency[y])
                            enqueue!(q_z,y)
                        end
                        n_paths[y] = n_paths[x]
                        ball[y] = ball[x]
                        #end_q += 1
                        #q[end_q] = y
                        new_end_cur += 1
                        push!(pred[y],x)
                        dist[y] = dist[x] + 1
                    elseif (ball[y] != ball[x])
                        have_to_stop = true
                        push!(sp_edges,(x,y))
                        #println("The bfs met at ",(x,y))

                    elseif (dist[y] == dist[x] +1)
                        n_paths[y] += n_paths[x]
                        push!(pred[y],x)
                    end
                    #println("q after visiting one neig ",q)
                end

            end
            if (sum_degs_cur == 0)
                have_to_stop = true
            end
            if (to_expand == @adjacency)
                sum_degs_s = sum_degs_cur
                start_s = start_cur
                end_s = new_end_cur
            else
                sum_degs_z = sum_degs_cur
                start_z = start_cur
                end_z = new_end_cur
            end
            #end_cur = new_end_cur
        end
        if (length(sp_edges) != 0)
            for p in sp_edges
                tot_weight += n_paths[p[1]] * n_paths[p[2]]
            end
            if (alpha_sampling > 0 && tot_weight > 1)
                num_path_to_sample = trunc(UInt128,floor(alpha_sampling * tot_weight))
            end
            num_paths = num_path_to_sample
            #Debug
            #println("Total Weight ", tot_weight," | PATHS TO SAMPLE: ",num_paths)

            for j in 1:num_path_to_sample
                path = Array{Int64}([])
                random_edge = rand(0:tot_weight-1)
                cur_edge = 0
                #println("PREDS ",pred)
                for p in sp_edges
                    cur_edge += n_paths[p[1]] * n_paths[p[2]]
                    if (cur_edge > random_edge)
                        println("TOT WEIGHT ",tot_weight," cur edge ",cur_edge," random edge ",random_edge)
                        #println("randm edge ",p)
                        #println("s ",s," z ",z," p ",p)
                        _backtrack_path!(s,z,p[1],path,n_paths,pred)
                        _backtrack_path!(s,z,p[2],path,n_paths,pred)
                        break
                    end
                end
                #println("PATH ",path)
                if (j == 1)
                    path_length = length(path)
                    begin
                        lock(lk)
                        try
                            shortest_path_length[path_length+1] +=1
                            if path_length+2 > new_diam_estimate[1]
                                new_diam_estimate[1] = path_length
                            end
                        finally 
                            unlock(lk)
                        end
                    end
                end
                for u in path
                    if haskey(path_map,u)
                        path_map[u] +=1
                    else
                        path_map[u] = 1
                    end
                    #=
                        if ball[u] == @visited_s
                            if run_perc
                                percolated_path_length[dist[u]+1] += ramp(percolation_states[s],percolation_states[z])/percolation_data[2][u]
                            end
                        else
                            if run_perc
                                longest_dist = dist[sp_edges[1][1]] + dist[sp_edges[1][2]]
                                percolated_path_length[longest_dist-dist[u]+1] += ramp(percolation_states[s],percolation_states[z])/percolation_data[2][u]
                            end
                        end
                    =#
                end
            end
            
        end
        # Updating Percolation centrality, wimpy variance and c-MCERA
        begin
            lock(lk)
            try
                for u in keys(path_map)
                    pm_value = path_map[u]
            
                  
                    percolation_value = pm_value/num_path_to_sample

                    if uniform_sampling
                        percolation_value =  percolation_value * (ramp(percolation_states[s],percolation_states[z])/percolation_data[1])  
                    end
                    percolation_centrality[u] += percolation_value
                    #println(" percolation_centrality of  ",u," increasing by ",percolation_value)
                    #betweenness[u] += pm_value/num_path_to_sample
                    # check whether the denominator must be squared too
                    wimpy_variance[u] += percolation_value*percolation_value 
                    #wimpy_variance[u] += pm_value*pm_value/num_path_to_sample/num_path_to_sample

                    # Updating c-Monte Carlo Empirical Rademacher (Only during the estimation phase)
                    if !boostrap_phase
                        for j in 1:mc_trials
                            mcrade[(u*mc_trials) + j] += lambdas[j] * percolation_value
                            #mcrade[(u*mc_trials) + j] += lambdas[j] * pm_value/num_path_to_sample

                        end
                    end
                    # Normalizing the percolated_path_length
                    #=
                        if run_perc
                            percolated_path_length[u] = percolated_path_length[u]/num_path_to_sample
                        end
                    =#
                end
            finally
                unlock(lk)
            end
        end
    end
    return nothing
end



#=

# Exact algorithm for our definition of percolation centrality

function parallel_percolation_centrality_new_target(g,percolation_states::Array{Float64},subgraph = nothing)::Tuple{Array{Float64},Array{Float64},Float64,Array{Float64},Float64,Float64}
    @assert nv(g) == lastindex(percolation_states) "Vertex set and percolation array must have the same dimensions"
    dv,_,_,_=compute_d_max(nv(g) ,percolation_states)
    max_d_v::Float64 = maximum(dv)
    @info("----------------------------------------| Stats |--------------------------------------------------")
    @info("Analyzing graph")
    @info("Number of nodes "*string(nv(g)))
    @info("Number of edges "*string(ne(g)))
    @info("Directed ? "*string(is_directed(g)))
    @info("Maximum Percolated state "*string(maximum(percolation_states)))
    @info("Minimum Percolated state "*string(minimum(percolation_states)))
    @info("Average Percolated state "*string(mean(percolation_states))*" std "*string(std(percolation_states)))
    @info("Maximum d_v "*string(max_d_v))
    @info("Using "*string(nthreads())* " Threads")
    @info("---------------------------------------------------------------------------------------------------")
    @info("Computing percolation centrality considering source and targets")
    flush(stderr)
    H = g
    induced_sg::Bool = false
    if !isnothing(subgraph)
        H = induced_subgraph(g,subgraph)
        induced_sg = true
    end
    start_time::Float64 = time()
    n::Int64 = nv(g)
    n_old::Int64 = n
    #percent_status::Array{Tuple{Float64,Int64,String}} = [(0.25,trunc(Int64,n*0.25),"Analyzed 25% of the graph"),(0.50,trunc(Int64,n*0.50),"Analyzed 50% of the graph"),(0.75,trunc(Int64,n*0.75),"Analyzed 75% of the graph"),(0.9,trunc(Int64,n*0.9),"Analyzed 90% of the graph")]
    #percent_index::Int16 = 1
    vs_active = [i for i in 1:n]
    d, r = divrem(n, nthreads())
    ntasks = d == 0 ? r : nthreads()
    task_size = cld(n, ntasks)
    percolation::Array{Array{Float64}} = [zeros(Float64,n) for _ in 1:ntasks]
    final_percolation::Array{Float64} = zeros(Float64,n)
    #tmp_perc_states::Array{Float64} = copy(percolation_states)
    #percolation_data::Tuple{Float64,Array{Float64}} = percolation_differences(sort(tmp_perc_states),n)    
    tmp_perc_states::Dict{Int64,Float64} = Dict(v => percolation_states[v] for v in 1:n)
    sorted_dict = OrderedDict(sort(collect(tmp_perc_states), by = kv -> kv[2]))
    percolation_data::Tuple{Float64,Dict{Int64,Float64}} = percolation_differences(sorted_dict,n)
    verbose_step::Int64 = trunc(Int64,floor(n*0.25))
    processed_so_far::Int64 = 0
    if induced_sg
        g = H[1] 
        n = nv(g)
        percolation_states = percolation_states[subgraph]
        percolation = [zeros(Float64,n) for _ in 1:ntasks]
        final_percolation = zeros(Float64,n)
        task_size = cld(n, ntasks)
    end
    @sync for (t, task_range) in enumerate(Iterators.partition(1:n, task_size))
        Threads.@spawn for s in @view(vs_active[task_range])
            if length(outneighbors(g,s)) > 0 && percolation_states[s] > 0.0
                _parallel_sz_bfs_exact!(g,s,percolation_states,percolation[t])
            end
            #=
            if (Sys.free_memory() / Sys.total_memory() < 0.1)
                clean_gc()
                sleep(0.01)DimensionMismatch
            end
            =#
            processed_so_far +=1
            if (verbose_step > 0 && processed_so_far % verbose_step == 0)
                finish_partial::String = string(round(time() - start_time; digits=4))
                time_to_finish::String = string(round((n*(time() - start_time) / processed_so_far )-(time() - start_time) ; digits=4))
                @info("Percolation Centrality. Processed " * string(processed_so_far) * "/" * string(n) * " nodes in " * finish_partial * " seconds | Est. remaining time : "*time_to_finish)
                flush(stderr)
            end
        end
    end
    final_percolation = reduce(+,percolation)
    # Normalizing using percolation states 
    unnormalized_scores::Array{Float64} = copy(final_percolation)
    for v in 1:n
        final_percolation[v] = final_percolation[v]/percolation_data[1]

    end
    finish_time::Float64 = time()-start_time
    @info("Percolation centrality s to z computed in "*string(finish_time))
    flush(stderr)
    if induced_sg
        tmp_final_percolation::Array{Float64} = zeros(Float64,n_old)
        println(H[2])
        println(final_percolation)
        println(percolation_states)
        for i in 1:lastindex(H[2])

            tmp_final_percolation[H[2][i]] = final_percolation[i]
        end
        final_percolation = tmp_final_percolation
    end
    return final_percolation,unnormalized_scores,percolation_data[1],dv,max_d_v,finish_time
end

=#


# ERA Approximation algorithm



function parallel_estimate_percolation_centrality_era_new(g,percolation_states::Array{Float64},epsilon::Float64,delta::Float64,initial_sample::Int64 = 0,geo::Float64 = 1.2,sample_size_diam::Int64 = 256,vc_upperbound::Bool = false )
    n::Int64 = nv(g)
    m::Int64 = ne(g)
    directed::Bool = is_directed(g)
    dv,_,_,_=compute_d_max(nv(g) ,percolation_states)
    max_d_v::Float64 = maximum(dv)
    @info("----------------------------------------| Stats |--------------------------------------------------")
    @info("Analyzing graph")
    @info("Number of nodes "*string(n))
    @info("Number of edges "*string(m))
    @info("Directed ? "*string(directed))
    @info("Maximum Percolated state "*string(maximum(percolation_states)))
    @info("Minimum Percolated state "*string(minimum(percolation_states)))
    @info("Average Percolated state "*string(mean(percolation_states))*" std "*string(std(percolation_states)))
    @info("Maximum d_v "*string(max_d_v))
    @info("Epsilon "*string(epsilon)*" Delta "*string(delta))
    @info("Using "*string(nthreads())* " Threads")
    @info("---------------------------------------------------------------------------------------------------")
    flush(stderr)
    ntasks = nthreads()
    if initial_sample == 0
        @info("Inferring the size of the first sample in the schedule")
        initial_sample = trunc(Int64,floor((1+8*epsilon + sqrt(1+16*epsilon)*log(6/delta))/(4*epsilon*epsilon)))
        @info("The size of the first sample in the schedule is "*string(initial_sample))
        flush(stderr)
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
    #tmp_perc_states::Array{Float64} = copy(percolation_states)
    diam::Int64 = 0
    tmp_perc_states::Dict{Int64,Float64} = Dict(v => percolation_states[v] for v in 1:n)
    sorted_dict = OrderedDict(sort(collect(tmp_perc_states), by = kv -> kv[2]))
    percolation_data::Tuple{Float64,Dict{Int64,Float64}} = percolation_differences(sorted_dict,n)
    #percolation_data::Tuple{Float64,Array{Float64}} = percolation_differences(sort(tmp_perc_states),n)
    if vc_upperbound
        @info("Approximating diameter using Random BFS algorithm")
        flush(stderr)
        diam,time_diam = parallel_random_bfs(g,sample_size_diam)
        finish_diam::String = string(round(time_diam; digits=4))
        @info("Estimated diameter "*string(diam)*" in "*finish_diam*" seconds")
        flush(stderr)
        max_sample::Int64 = trunc(Int64,0.5/epsilon/epsilon * (log2(diam-1)+1+log(2/delta)))
        @info("VC-Upperbound on the maximum sample size "*string(max_sample))
        flush(stderr)
    end
    if diam == 0
        max_sample = typemax(Int64)
        @info("Deterministic ERA Upperbound on the maximum sample size ")
        flush(stderr)
           
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
                _parallel_sz_bfs_new!(g,percolation_states,percolation_data,local_B[t],local_B_1[t],local_B_2[t])
            end
        end
        sampled_so_far += sample_i
        B_vectorized = reduce_dictionary(local_B)
        if length(B_vectorized)>0
            omega = compute_xi(B_vectorized .*[n*(n-1)],sample_size_schedule[j])
            delta_i = delta/2^k
            xi = 2*omega + (log(3/delta_i)+sqrt((log(3/delta_i)+4*sample_size_schedule[j]*omega)*log(3/delta_i)))/sample_size_schedule[j]  + sqrt(log(3/delta_i)/(2*sample_size_schedule[j]))
        else
            xi = Inf
        end
        @info("ERA Upperbound "*string(xi)*" Target SD "*string(epsilon)*" #Sampled pairs "*string(sample_size_schedule[j]) *" in "*string(round(time()-start_time;digits = 4))*" Next Sample size "*string(trunc(Int,geo^(k+1)*sample_size_schedule[2])))
        flush(stderr)
        if xi <= epsilon || sample_size_schedule[j] >= max_sample
            keep_sampling = false
        else
            j+=1
        end

    end
    final_B_1 = reduce(+,local_B_1) .*[(n*(n-1))/sample_size_schedule[end]]
    finish_time::Float64 = time()-start_time
    @info("Converged! Sampled "*string(sample_size_schedule[end])*"/"*string(max_sample)*" couples in "*string(round(finish_time;digits = 4))*" seconds ")
    flush(stderr)
    return final_B_1,sample_size_schedule,max_sample,xi,finish_time
end



function _parallel_sz_bfs_new!(g,percolation_states::Array{Float64},percolation_data::Tuple{Float64,Dict{Int64,Float64}},B::Dict{Float64,Float64},B_1::Array{Float64},B_2::Array{Float64})
    n::Int64 = nv(g)
    q::Queue{Int64} = Queue{Int64}()
    ball::Array{Int16} = zeros(Int16,n)
    n_paths::Array{UInt128} = zeros(UInt128,n)
    dist::Array{Int64} = zeros(Int64,n)
    pred::Array{Array{Int64}} = [Array{Int64}([]) for _ in 1:n]

    s::Int64 = sample(1:n)
    z::Int64 = sample(1:n)
    w::Int64 = 0
    d_z_min::Float64 = Inf
    #q_backtrack::Stack{Int64} = Stack{Int64}()
    while (s == z)
        z = sample(1:n)
    end
    if percolation_states[s] > percolation_states[z]
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
                            #=
                            if length(q_backtrack) == 0
                                enqueue!(q_backtrack,z)
                            end
                            =#
                        end
                        push!(pred[v],w)
                        enqueue!(q,v)
                        #push!(q_backtrack,v)
                    elseif (dist[v] == dist[w] + 1)
                        n_paths[v] += n_paths[w]
                        push!(pred[v],w)
                    end

                end
            end
        end
        #backtrack
        if ball[z] == 1
            q_backtrack::Queue{Int64} = Queue{Int64}()
            queued::Array{Int16} = zeros(Int16,n)
            #number_of_sps_to_target::Dict{Int64,Int128} = Dict{Int64,Int128}()
            number_of_sps_to_target::Array{UInt128} = zeros(UInt128,n)
            number_of_paths_through_curr::UInt128 = 0
            for p in pred[z]
                enqueue!(q_backtrack,p)
                number_of_sps_to_target[p] = 1
                queued[p] = 1
            end
            w = dequeue!(q_backtrack)
            while w != s
                number_of_paths_through_curr = n_paths[w] * number_of_sps_to_target[w]
                summand = (number_of_paths_through_curr/n_paths[z]) *(ramp(percolation_states[s],percolation_states[z])/percolation_data[1])  
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
                
                

                for p in pred[w]
                    if queued[p] == 0
                        enqueue!(q_backtrack,p)
                        queued[p] = 1
                    end
                    number_of_sps_to_target[p] += number_of_sps_to_target[w]
                end
                w  = dequeue!(q_backtrack)
            end
        end
    end

   return nothing
end



# Fixed sample size
function parallel_estimate_percolation_centrality_fixed_sample_size_new(g,percolation_states::Array{Float64},epsilon::Float64,delta::Float64,sample_size_diam::Int64 = 256 )
    n::Int64 = nv(g)
    m::Int64 = ne(g)
    directed::Bool = is_directed(g)
    dv,_,_,_=compute_d_max(nv(g) ,percolation_states)
    max_d_v::Float64 = maximum(dv)
    @info("----------------------------------------| Stats |--------------------------------------------------")
    @info("Analyzing graph")
    @info("Number of nodes "*string(n))
    @info("Number of edges "*string(m))
    @info("Directed ? "*string(directed))
    @info("Maximum Percolated state "*string(maximum(percolation_states)))
    @info("Minimum Percolated state "*string(minimum(percolation_states)))
    @info("Average Percolated state "*string(mean(percolation_states))*" std "*string(std(percolation_states)))
    @info("Epsilon "*string(epsilon)*" Delta "*string(delta))
    @info("Maximum d_v "*string(max_d_v))
    @info("Using "*string(nthreads())* " Threads")
    @info("---------------------------------------------------------------------------------------------------")
    flush(stderr)
    ntasks = nthreads()
    start_time::Float64 = time()
    
    

    final_percolation_centrality::Array{Float64} = zeros(Float64,n)
    percolation_centrality::Array{Array{Float64}} =  [zeros(Float64,n) for _ in 1:ntasks]
    
    #sg::static_graph = static_graph(adjacency_list(g),incidency_list(g))
    tmp_perc_states::Dict{Int64,Float64} = Dict(v => percolation_states[v] for v in 1:n)
    sorted_dict = OrderedDict(sort(collect(tmp_perc_states), by = kv -> kv[2]))
    percolation_data::Tuple{Float64,Dict{Int64,Float64}} = percolation_differences(sorted_dict,n)
    #tmp_perc_states::Array{Float64} = copy(percolation_states)
    #percolation_data::Tuple{Float64,Array{Float64}} = percolation_differences(sort(tmp_perc_states),n)
    @info("Approximating diameter using Random BFS algorithm")
    flush(stderr)
    diam,time_diam = parallel_random_bfs(g,sample_size_diam)
    finish_diam::String = string(round(time_diam; digits=4))
    @info("Estimated diameter "*string(diam)*" in "*finish_diam*" seconds")
    flush(stderr)
    max_sample::Int64 = trunc(Int64,0.5/epsilon/epsilon * (log2(diam-1)+1+log(2/delta)))
    @info("Maximum sample size "*string(max_sample))
    flush(stderr)
    task_size = cld(max_sample, ntasks)
    vs_active = [i for i in 1:max_sample]
    @sync for (t, task_range) in enumerate(Iterators.partition(1:max_sample, task_size))
        Threads.@spawn for _ in @view(vs_active[task_range])
             Random.seed!(1234 + t)
             _parallel_sz_bfs_fss_new!(g,percolation_states,percolation_data,percolation_centrality[t])
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
    flush(stderr)
    return final_percolation_centrality.*[(n*(n-1))/max_sample],max_sample,finish_time
end


function _parallel_sz_bfs_fss_new!(g,percolation_states::Array{Float64},percolation_data::Tuple{Float64,Dict{Int64,Float64}},percolation_centrality::Array{Float64},lk::ReentrantLock,uniform_sampling::Bool)
    n::Int64 = nv(g)
    s::Int64 = sample(1:n)
    z::Int64 = sample(1:n)
    w::Int64 = 0
    q::Queue{Int64} = Queue{Int64}()
    ball::Array{Int16} = zeros(Int16,n)
    n_paths::Array{UInt128} = zeros(UInt128,n)
    dist::Array{Int64} = zeros(Int64,n)
    pred::Array{Array{Int64}} = [Array{Int64}([]) for _ in 1:n]
    d_z_min::Float64 = Inf
    tot_weight::Int64 = 0
    random_edge::Int64 = 0
    cur_edge::Int64 = 0
    path::Array{Int64} = Array{Int64}([])
    q_backtrack::Queue{Int64} = Queue{Int64}()
    if !uniform_sampling
        s,z = weighted_sample_kappa(percolation_states)
    else
        while (s == z)
            z = sample(1:n)
        end
    end
    if percolation_states[s] > percolation_states[z]
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
                    #println("TOT WEIGHT ",tot_weight," cur edge ",cur_edge," random edge ",random_edge)

                    _backtrack_path!(s,z,p,path,n_paths,pred)
                    break
                end
            end
            begin
                lock(lk)
                try
                    for w in path
                        if !uniform_sampling
                            percolation_centrality[w] += 1
                        else
                            percolation_centrality[w] += (ramp(percolation_states[s],percolation_states[z])/percolation_data[1])  
                        end
                    end
                finally 
                    unlock(lk)
                end
            end    
            
        end  
    end

   return nothing
end



#Uniform and non Uniform Samplers Fixed Sample Size

# Non-Uniform 


function parallel_estimate_percolation_centrality_non_uniform(g,percolation_states::Array{Float64},sample_size::Int64,alpha_sampling::Float64 = 0.1)
    @assert nv(g) == lastindex(percolation_states) "Length of the percolation state array must be the same as the number of nodes"
    n::Int64 = nv(g)
    m::Int64 = ne(g)
    directed::Bool = is_directed(g)
    dv,_,_,_= compute_d_max(nv(g) ,percolation_states)
    max_d_v::Float64 = maximum(dv)
    @info("----------------------------------------| Stats |--------------------------------------------------")
    @info("Analyzing graph")
    @info("Number of nodes "*string(n))
    @info("Number of edges "*string(m))
    @info("Directed ? "*string(directed))
    @info("Maximum Percolated state "*string(maximum(percolation_states)))
    @info("Minimum Percolated state "*string(minimum(percolation_states)))
    @info("Average Percolated state "*string(mean(percolation_states))*" std "*string(std(percolation_states)))
    @info("Sample Size "*string(sample_size))
    @info("Sampler : Non Uniform")
    @info("Maximum d_v "*string(max_d_v))
    @info("Using "*string(nthreads())* " Threads")
    @info("---------------------------------------------------------------------------------------------------")
    flush(stderr)
    ntasks = nthreads()
    sg::static_graph = static_graph(adjacency_list(g),incidency_list(g))
    run_perc::Bool = true
    mc_trials = 1
    final_percolation_centrality::Array{Float64} = zeros(Float64,n)
    #wimpy_variance::Array{Array{Float64}} = [zeros(Float64,n) for _ in 1:ntasks]
    final_wimpy_variance::Array{Float64} = zeros(Float64,n)
    #shortest_path_length::Array{Array{Int64}} = [zeros(Int64,n+1) for _ in 1:ntasks]
    final_shortest_path_length::Array{Int64} = zeros(Int64,n+1)
    #mcrade::Array{Array{Float64}} = [zeros(Float64,(n+1)*mc_trials) for _ in 1:ntasks]
    final_mcrade::Array{Float64} = zeros(Float64,(n+1)*mc_trials)
    tmp_perc_states::Dict{Int64,Float64} = Dict(v => percolation_states[v] for v in 1:n)
    sorted_dict = OrderedDict(sort(collect(tmp_perc_states), by = kv -> kv[2]))
    percolation_data::Tuple{Float64,Dict{Int64,Float64}} = percolation_differences(sorted_dict,n)
    #tmp_perc_states::Array{Float64} = copy(percolation_states)
    #percolation_data::Tuple{Float64,Array{Float64}} = percolation_differences(sort(tmp_perc_states),n)
    #betweenness::Array{Float64} = zeros(Float64,n)
    start_time::Float64 = time()
    lk::ReentrantLock = ReentrantLock()
    task_size = cld(sample_size, ntasks)
    vs_active = [i for i in 1:sample_size]
    @sync for (t, task_range) in enumerate(Iterators.partition(1:sample_size, task_size))
        Threads.@spawn for _ in @view(vs_active[task_range])
            #_parallel_random_path_weighted_lk!(sg,n,final_percolation_centrality,final_wimpy_variance,percolation_states,percolation_data,final_shortest_path_length,final_mcrade,1,alpha_sampling,final_new_diam_estimate,lk,run_perc,false)
            #_parallel_sz_bfs_replace!(g,percolation_states,percolation_data,final_percolation_centrality,lk,false)
            #_parallel_sz_bfs_fss_new!(g,percolation_states,percolation_data,final_percolation_centrality,lk,false)
            _random_path_non_unif!(sg,n, final_percolation_centrality,percolation_states,percolation_data,alpha_sampling,lk,false)
        end
    end
    finish_time::Float64 = time()-start_time
    @info("Completed in "*string(round(finish_time,digits = 4))) 
    flush(stderr)
    return final_percolation_centrality .*[1/sample_size],sample_size,sample_size,finish_time,finish_time

end



function parallel_estimate_percolation_centrality_uniform(g,percolation_states::Array{Float64},sample_size::Int64,alpha_sampling::Float64 = 0.1)
    @assert nv(g) == lastindex(percolation_states) "Length of the percolation state array must be the same as the number of nodes"
    n::Int64 = nv(g)
    m::Int64 = ne(g)
    directed::Bool = is_directed(g)
    dv,_,_,_= compute_d_max(nv(g) ,percolation_states)
    max_d_v::Float64 = maximum(dv)
    @info("----------------------------------------| Stats |--------------------------------------------------")
    @info("Analyzing graph")
    @info("Number of nodes "*string(n))
    @info("Number of edges "*string(m))
    @info("Directed ? "*string(directed))
    @info("Maximum Percolated state "*string(maximum(percolation_states)))
    @info("Minimum Percolated state "*string(minimum(percolation_states)))
    @info("Average Percolated state "*string(mean(percolation_states))*" std "*string(std(percolation_states)))
    @info("Sample Size "*string(sample_size))
    @info("Sampler : Uniform")
    @info("Maximum d_v "*string(max_d_v))
    @info("Using "*string(nthreads())* " Threads")
    @info("---------------------------------------------------------------------------------------------------")
    flush(stderr)
    ntasks = nthreads()
    sg::static_graph = static_graph(adjacency_list(g),incidency_list(g))
    run_perc::Bool = true
    mc_trials = 1
    final_percolation_centrality::Array{Float64} = zeros(Float64,n)
    #wimpy_variance::Array{Array{Float64}} = [zeros(Float64,n) for _ in 1:ntasks]
    final_wimpy_variance::Array{Float64} = zeros(Float64,n)
    #shortest_path_length::Array{Array{Int64}} = [zeros(Int64,n+1) for _ in 1:ntasks]
    final_shortest_path_length::Array{Int64} = zeros(Int64,n+1)
    #mcrade::Array{Array{Float64}} = [zeros(Float64,(n+1)*mc_trials) for _ in 1:ntasks]
    final_mcrade::Array{Float64} = zeros(Float64,(n+1)*mc_trials)
    tmp_perc_states::Dict{Int64,Float64} = Dict(v => percolation_states[v] for v in 1:n)
    sorted_dict = OrderedDict(sort(collect(tmp_perc_states), by = kv -> kv[2]))
    percolation_data::Tuple{Float64,Dict{Int64,Float64}} = percolation_differences(sorted_dict,n)
    #tmp_perc_states::Array{Float64} = copy(percolation_states)
    #percolation_data::Tuple{Float64,Array{Float64}} = percolation_differences(sort(tmp_perc_states),n)
    #betweenness::Array{Float64} = zeros(Float64,n)
    start_time::Float64 = time()
    lk::ReentrantLock = ReentrantLock()
    task_size = cld(sample_size, ntasks)
    vs_active = [i for i in 1:sample_size]
    @sync for (t, task_range) in enumerate(Iterators.partition(1:sample_size, task_size))
        Threads.@spawn for _ in @view(vs_active[task_range])
            #_parallel_random_path_weighted_lk!(sg,n,final_percolation_centrality,final_wimpy_variance,percolation_states,percolation_data,final_shortest_path_length,final_mcrade,1,alpha_sampling,final_new_diam_estimate,lk,run_perc,true,true)
            #_parallel_sz_bfs_replace!(g,percolation_states,percolation_data,final_percolation_centrality,lk,true)
            #_parallel_sz_bfs_fss_new!(g,percolation_states,percolation_data,final_percolation_centrality,lk,true)
            _random_path_non_unif!(sg,n, final_percolation_centrality,percolation_states,percolation_data,alpha_sampling,lk,true)

        end
    end
    finish_time::Float64 = time()-start_time
    @info("Completed in "*string(round(finish_time,digits = 4))) 
    flush(stderr)
    return final_percolation_centrality .*[(n*(n-1))/sample_size],sample_size,sample_size,finish_time,finish_time

end



function _parallel_sz_bfs_replace!(g,percolation_states::Array{Float64},percolation_data::Tuple{Float64,Dict{Int64,Float64}},percolation_centrality::Array{Float64},lk::ReentrantLock,uniform_sampling::Bool)
    n::Int64 = nv(g)
    q::Queue{Int64} = Queue{Int64}()
    ball::Array{Int16} = zeros(Int16,n)
    n_paths::Array{UInt128} = zeros(UInt128,n)
    dist::Array{Int64} = zeros(Int64,n)
    pred::Array{Array{Int64}} = [Array{Int64}([]) for _ in 1:n]

    s::Int64 = sample(1:n)
    z::Int64 = sample(1:n)
    w::Int64 = 0
    d_z_min::Float64 = Inf
    #q_backtrack::Stack{Int64} = Stack{Int64}()
    if !uniform_sampling
        s,z = weighted_sample_kappa(percolation_states)
    else
        while (s == z)
            z = sample(1:n)
        end
    end
    if percolation_states[s] > percolation_states[z]
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
                        
                        end
                        push!(pred[v],w)
                        enqueue!(q,v)
                        #push!(q_backtrack,v)
                    elseif (dist[v] == dist[w] + 1)
                        n_paths[v] += n_paths[w]
                        push!(pred[v],w)
                    end

                end
            end
        end
        #backtrack

    if ball[z] == 1
            q_backtrack::Queue{Int64} = Queue{Int64}()
            queued::Array{Int16} = zeros(Int16,n)
            #number_of_sps_to_target::Dict{Int64,Int128} = Dict{Int64,Int128}()
            number_of_sps_to_target::Array{UInt128} = zeros(UInt128,n)
            number_of_paths_through_curr::UInt128 = 0
            for p in pred[z]
                enqueue!(q_backtrack,p)
                number_of_sps_to_target[p] = 1
                queued[p] = 1
            end
            w = dequeue!(q_backtrack)
            while w != s
                number_of_paths_through_curr = n_paths[w] * number_of_sps_to_target[w]
                begin
                    lock(lk)
                    try
                        if !uniform_sampling
                            percolation_centrality[w] += (number_of_paths_through_curr/n_paths[z])
                        else
                            percolation_centrality[w] += (number_of_paths_through_curr/n_paths[z]) *  ramp(percolation_states[s],percolation_states[z])/percolation_data[1]
                        end
                    finally 
                        unlock(lk)
                    end
                end    
                for p in pred[w]
                    if queued[p] == 0
                        enqueue!(q_backtrack,p)
                        queued[p] = 1
                    end
                    number_of_sps_to_target[p] += number_of_sps_to_target[w]
                end
                w  = dequeue!(q_backtrack)
            end
        end
    end


   return nothing
end



