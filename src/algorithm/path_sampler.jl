macro unvisited()
    return 0
end
macro visited_s()
    return 1
end
macro visited_z()
    return 2
end

macro adjacency()
    return 0
end

macro incidency()
    return 1
end


function _random_path!(sg::static_graph,n::Int64,q::Array{Int64},ball::Array{Int16},n_paths::Array{UInt128},dist::Array{Int64},pred::Array{Array{Int64}},num_paths::Array{Int64},percolation_centrality::Array{Float64},wimpy_variance::Array{Float64},percolation_states::Array{Float64},percolation_data::Tuple{Float64,Array{Float64}},shortest_path_length::Array{Int64},mcrade::Array{Float64},mc_trials::Int64,alpha_sampling::Float64,new_diam_estimate::Array{Int64},run_perc::Bool =true,boostrap_phase::Bool = false)

    end_q::UInt32 = 1
    tot_weight::Int64 = 0
    cur_edge::UInt64 = 0
    random_edge::UInt32 = 0
    s::Int64 = sample(1:n)
    z::Int64 = sample(1:n)
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
    maxval_half::Int64 = maxval_lambdas/2
    if !boostrap_phase
        for j in 1:mc_trials
            lambdas[j] = (sample(0:maxval_lambdas-1) >= maxval_half)*2-1
        end
    end
    while (s == z)
        z = sample(1:n)
    end
    q_s::Queue{Int64} = Queue{Int64}()
    q_z::Queue{Int64} = Queue{Int64}()
    enqueue!(q_s,s)
    enqueue!(q_z,z)

    end_q = 2
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
            start_cur = start_s
            end_cur = end_s
            start_s = end_q
            new_end_cur = end_s
            end_s = end_q
            sum_degs_s = 0
            sum_degs_cur = sum_degs_s 
            to_expand = @adjacency
            #println("SUM DEG S ",sum_degs_s)
        else
            start_cur = start_z
            end_cur = end_z
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
                    end_q += 1
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
        end_cur = new_end_cur
    end
    if (length(sp_edges) == 0)
        for u in 1:end_q
            ball[q[u]] = @unvisited
            dist[q[u]] = 0
            pred[q[u]] = Array{Int64}([])
            # Check this npaths and q, they should be set back to 0 once finished
            n_paths[q[u]] = 0
            q[u] = 0
        end
    else
        for p in sp_edges
            tot_weight += n_paths[p[1]] * n_paths[p[2]]
        end
        if (alpha_sampling > 0 && tot_weight > 1)
            num_path_to_sample = trunc(Int64,floor(alpha_sampling * tot_weight))
        end
        num_paths = num_path_to_sample

        for j in 1:num_path_to_sample
            path = Array{Int64}([])
            random_edge = rand(0:tot_weight-1)
            cur_edge = 0
            #println("PREDS ",pred)
            for p in sp_edges
                cur_edge += n_paths[p[1]] * n_paths[p[2]]
                if (cur_edge > random_edge)
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
                shortest_path_length[path_length+1] +=1
                if path_length+2 > new_diam_estimate[1]
                    new_diam_estimate[1] = path_length
                end
            end
            for u in path
                if haskey(path_map,u)
                    path_map[u] +=1
                else
                    path_map[u] = 1
                end
            end
        end

        for u in 1:end_q
            ball[q[u]] = @unvisited
            dist[q[u]] = 0
            pred[q[u]] = Array{Int64}([])
            # Check this npaths and q, they should be set back to 0 once finished
            n_paths[q[u]] = 0
            q[u] = 0
        end
    end
    # Updating Percolation centrality, wimpy variance and c-MCERA
    #println("PATH MAP ",path_map," PERC STATES ",percolation_states)
    for u in keys(path_map)
        pm_value = path_map[u]
        #println("PM VAL ",pm_value, " RAMP ",ramp(percolation_states[s],percolation_states[z])," DENOM ",percolation_data[2][u], " num path sample ",num_path_to_sample)
        if run_perc
            percolation_value = pm_value*(ramp(percolation_states[s],percolation_states[z])/percolation_data[2][u]) /num_path_to_sample
        else
            percolation_value = pm_value/num_path_to_sample
        end
        percolation_centrality[u] += percolation_value
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

    return nothing
end


function _parallel_random_path!(sg::static_graph,n::Int64,percolation_centrality::Dict{Int64,Float64},wimpy_variance::Array{Float64},percolation_states::Array{Float64},percolation_data::Tuple{Float64,Array{Float64}},shortest_path_length::Array{Int64},mcrade::Array{Float64},mc_trials::Int64,alpha_sampling::Float64,new_diam_estimate::Array{Int64},run_perc::Bool = true,boostrap_phase::Bool = false)

    #q::Array{Int64} = zeros(Int64,n)
    ball::Array{Int16} = zeros(Int16,n)
    n_paths::Array{UInt128} = zeros(UInt128,n)
    dist::Array{Int64} = zeros(Int64,n)
    pred::Array{Array{Int64}} = [Array{Int64}([]) for _ in 1:n]
    #num_paths::Array{Int64} = zeros(Int64,n)
    end_q::UInt32 = 1
    tot_weight::Int64 = 0
    cur_edge::UInt64 = 0
    random_edge::UInt32 = 0;
    s::Int64 = sample(1:n)
    z::Int64 = sample(1:n)
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
    maxval_half::Int64 = maxval_lambdas/2
    if !boostrap_phase
        for j in 1:mc_trials
            lambdas[j] = (sample(0:maxval_lambdas-1) >= maxval_half)*2-1
        end
    end
    while (s == z)
        z = sample(1:n)
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
            start_cur = start_s
            end_cur = end_s
            start_s = end_q
            new_end_cur = end_s
            end_s = end_q
            sum_degs_s = 0
            sum_degs_cur = sum_degs_s 
            to_expand = @adjacency
            #println("SUM DEG S ",sum_degs_s)
        else
            start_cur = start_z
            end_cur = end_z
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
        end_cur = new_end_cur
    end
    if (length(sp_edges) != 0)
        for p in sp_edges
            tot_weight += n_paths[p[1]] * n_paths[p[2]]
        end
        if (alpha_sampling > 0 && tot_weight > 1)
            num_path_to_sample = trunc(Int64,floor(alpha_sampling * tot_weight))
        end
        num_paths = num_path_to_sample

        for j in 1:num_path_to_sample
            path = Array{Int64}([])
            random_edge = rand(0:tot_weight-1)
            cur_edge = 0
            #println("PREDS ",pred)
            for p in sp_edges
                cur_edge += n_paths[p[1]] * n_paths[p[2]]
                if (cur_edge > random_edge)
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
                shortest_path_length[path_length+1] +=1
                if path_length+2 > new_diam_estimate[1]
                    new_diam_estimate[1] = path_length
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
        #=x - y
           println("a ",a," x ",x)
           @reduce(s += a, t += b)
  
        for u in 1:end_q
            ball[q[u]] = @unvisited
            dist[q[u]] = 0
            pred[q[u]] = Array{Int64}([])
            # Check this npaths and q, they should be set back to 0 once finished
            n_paths[q[u]] = 0
            q[u] = 0
        end
        =#
    end
    # Updating Percolation centrality, wimpy variance and c-MCERA
    #println("PATH MAP ",path_map," PERC STATES ",percolation_states)
    for u in keys(path_map)
        pm_value = path_map[u]
        #println("PM VAL ",pm_value, " RAMP ",ramp(percolation_states[s],percolation_states[z])," DENOM ",percolation_data[2][u], " num path sample ",num_path_to_sample)
        if run_perc
            # κ(s,z,v)
            percolation_value = pm_value* (ramp(percolation_states[s],percolation_states[z])/percolation_data[2][u]) /num_path_to_sample

            #percolation_value = pm_value*(ramp(percolation_states[s],percolation_states[z])/percolation_data[2][u]) /num_path_to_sample
        else
            percolation_value = pm_value/num_path_to_sample
        end
        percolation_centrality[u] += percolation_value
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

    return nothing
end

function _parallel_random_path_lk!(sg::static_graph,n::Int64,percolation_centrality::Array{Float64},wimpy_variance::Array{Float64},percolation_states::Array{Float64},percolation_data::Tuple{Float64,Dict{Int64,Float64}} ,shortest_path_length::Array{Int64},mcrade::Array{Float64},mc_trials::Int64,alpha_sampling::Float64,new_diam_estimate::Array{Int64},lk::ReentrantLock,run_perc::Bool = true,boostrap_phase::Bool = false)

    q::Array{Int64} = zeros(Int64,n)
    ball::Array{Int16} = zeros(Int16,n)
    n_paths::Array{UInt128} = zeros(UInt128,n)
    dist::Array{Int64} = zeros(Int64,n)
    pred::Array{Array{Int64}} = [Array{Int64}([]) for _ in 1:n]
    #num_paths::Array{Int64} = zeros(Int64,n)
    end_q::UInt32 = 1
    tot_weight::Int64 = 0
    cur_edge::UInt64 = 0
    random_edge::UInt32 = 0;
    s::Int64 = sample(1:n)
    z::Int64 = sample(1:n)
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
    maxval_half::Int64 = maxval_lambdas/2
    if !boostrap_phase
        for j in 1:mc_trials
            lambdas[j] = (sample(0:maxval_lambdas-1) >= maxval_half)*2-1
        end
    end
    while (s == z)
        z = sample(1:n)
    end

    end_q = 2
    q[1] = s
    q[2] = z
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
            start_cur = start_s
            end_cur = end_s
            start_s = end_q
            new_end_cur = end_s
            end_s = end_q
            sum_degs_s = 0
            sum_degs_cur = sum_degs_s 
            to_expand = @adjacency
            #println("SUM DEG S ",sum_degs_s)
        else
            start_cur = start_z
            end_cur = end_z
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
            x = q[start_cur]
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
                    else
                        sum_degs_cur += lastindex(sg.incidency[y])
                    end
                    n_paths[y] = n_paths[x]
                    ball[y] = ball[x]
                    end_q += 1
                    q[end_q] = y
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
        end_cur = new_end_cur
    end
    if (length(sp_edges) != 0)
        for p in sp_edges
            tot_weight += n_paths[p[1]] * n_paths[p[2]]
        end
        if (alpha_sampling > 0 && tot_weight > 1)
            num_path_to_sample = trunc(Int64,floor(alpha_sampling * tot_weight))
        end
        num_paths = num_path_to_sample

        for j in 1:num_path_to_sample
            path = Array{Int64}([])
            random_edge = rand(0:tot_weight-1)
            cur_edge = 0
            #println("PREDS ",pred)
            for p in sp_edges
                cur_edge += n_paths[p[1]] * n_paths[p[2]]
                if (cur_edge > random_edge)
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
                shortest_path_length[path_length+1] +=1
                if path_length+2 > new_diam_estimate[1]
                    new_diam_estimate[1] = path_length
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
        #=x - y
           println("a ",a," x ",x)
           @reduce(s += a, t += b)
  
        for u in 1:end_q
            ball[q[u]] = @unvisited
            dist[q[u]] = 0
            pred[q[u]] = Array{Int64}([])
            # Check this npaths and q, they should be set back to 0 once finished
            n_paths[q[u]] = 0
            q[u] = 0
        end
        =#
    end
    # Updating Percolation centrality, wimpy variance and c-MCERA
    #println("PATH MAP ",path_map," PERC STATES ",percolation_states)
    begin
        lock(lk)
        try
            for u in keys(path_map)
                pm_value = path_map[u]
                #println("PM VAL ",pm_value, " RAMP ",ramp(percolation_states[s],percolation_states[z])," DENOM ",percolation_data[2][u], " num path sample ",num_path_to_sample)
                if run_perc
                    percolation_value = pm_value*(ramp(percolation_states[s],percolation_states[z])/percolation_data[2][u]) /num_path_to_sample
                else
                    percolation_value = pm_value/num_path_to_sample
                end
                percolation_centrality[u] += percolation_value
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

    return nothing
end

function _backtrack_path!(s::Int64,z::Int64,w::Int64,path::Array{Int64},n_paths::Array{UInt128},pred::Array{Array{Int64}})
    tot_weight::UInt128 = n_paths[w]
    random_pred::UInt128 = 0
    cur_pred::UInt128 = 0
    v::Int64 = 0
    if (w == s || w == z)
        return nothing
    end
    #=
    tot_weight_test::UInt128 = 0
    for p in pred[w]
        tot_weight_test += n_paths[p]
    end
    if tot_weight_test != tot_weight
        println("IL TEST NON PASSA!!!!")
        exit(1)
    end
    =#
    push!(path,w)
    random_pred = rand(0:tot_weight-1)
    for p in pred[w]
        v = p
        cur_pred += n_paths[v]
        if (cur_pred > random_pred)
            break
        end
    end
    if (s != v && z != v)
        _backtrack_path!(s,z,v,path,n_paths,pred)
    end
end


