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


function _random_path!(sg::static_graph,n::Int64,q::Array{Int64},ball::Array{Int16},n_paths::Array{Int64},dist::Array{Int64},pred::Array{Int64,Array{Int64}},q::Array{Int64},num_paths::Array{Int64},percolation_cepercolation_centralityntrality::Array{Float64},wimpy_varaince::Array{Float64},percolation_states::Array{Float64},percolation_data::Tuple{Float64,Array{Float64}})

    end_q::UInt32 = 1
    tot_weight::UInt64 = 0
    cur_edge::UInt64 = 0
    random_edge::UInt32 = 0;
    s::UInt64 = sample(1:n)
    z::UInt64 = sample(1:n)
    #ball::Array{Int16} = zeros(Int32,n)
    x::UInt32  = 0
    y::UInt32 = 0
    have_to_stop::Boolean = false;
    start_s::UInt64 = 1
    start_z::UInt64 = 2
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
    num_path_to_sample::Int64 = 1
    num_paths::Int64 = 0
    sp_edges::Array{Tuple{Int64,Int64}} = Array{Tuple{Int64,Int64}}([])
    path::Array{Int64} = Array{Int64}([])
    path_map::Dict{Int64,Int64} = Dict{Int64,Int64}()
    pm_value::Int64 = 0
    while (s == z)
        z = sample(1:n)
    end

    end_q = 3
    q[1] = s
    q[2] = z

    ball[s] = @visited_s
    bass[z] = @visited_z

    n_paths[s] = 1
    n_paths[z] = 1
    dist[s] = 0
    dist[0] = 0

    while (!have_to_stop)
        
        # decide what ball we have to expand
        if (sum_degs_s <= sum_degs_z)
            start_cur = start_s
            end_cur = end_s
            start_s = end_q
            new_end_cur = end_s
            end_s = end_q
            sum_degs_s = 0
            to_expand = @adjacency
        else
            start_cur = start_z
            end_cur = end_z
            start_z = end_q
            new_end_cur = end_z
            end_z = end_q
            sum_degs_z = 0
            to_expand = @incidency
        end

        while (start_cur < end_cur)
            # wip
            x = q[start_cur+=1]
            if (to_expand == @adjacency)
                neigh_num = lastindex(sg.adjacency[x])
            else
                neigh_num = lastindex(sg.incidency[x])
            end
            for j in 1:neigh_num
                vis_edges +=1
                if (to_expand == @adjacency)
                    y = sg.adjacency[x][j]
                else
                    y = sg.incidency[x][j]
                end
                if (ball[y] == @unvisited)
                    if (to_expand == @adjacency)
                        sum_degs_cur += lastindex(sg.adjacency[y])
                    else
                        sum_degs_cur += lastindex(sg.incidency[y])
                    end
                    n_paths[y] = n_paths[x]
                    ball[y] = ball[x]
                    q[end_q+=1] = y
                    new_end_cur += 1
                    push!(pred[y],x)
                    dist[y] = dist[x] + 1
                elseif (ball[y] != ball[x])
                    have_to_stop = true
                    push!(sp_edges,(x,y))

                elseif (dist[y] == dist[x] +1)
                    n_paths[y] += n_paths[x]
                    push!(pred[y],x)
                end
               
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
    end
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

        for p in sp_edges
            cur_edge += n_paths[p[1]] * n_paths[p[2]]
            if (cur_edge > random_edge)
                _backtrack_path!(s,z,p[1],path)
                _backtrack_path!(s,z,p[2],path)
                break
            end
        end
        if (j == 1)
            path_length = length(path)
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
    # Updating Percolation and wimpy variance
    for u in keys(path_map)
        pm_value =  path_map[u]
        percolation_centrality[u] += pm_value*(ramp(percolation_states[s],percolation_states[z])/percolation_data[2][u]) /num_path_to_sample
        # check whether the denominator must be squared too
        wimpy_varaince[u] += (pm_value*(ramp(percolation_states[s],percolation_states[z])/percolation_data[2][u]))^2 /num_path_to_sample
    end
    # Updating c-Monte Carlo Rademacher (todo)

    # Updating average percolated path length (todo)
    return nothing
end
# We need to add all the betweenness update part in the bidirectional bfs

function _backtrack_path!(s::Int64,z::Int64,w::Int64,path::Array{Int64})
    #wip
    tot_weight::UInt64 = n_paths[w]
    random_pred::Int64 = 0
    cur_pred::Int64 = 0
    v::Int64 = 0
    if (w == s || w == z)
        return nothing
    end
    random_pred = rand(0:tot_weight-1)
    for p in pred[w]
        v = p
        cur_pred += n_paths[z]
        if (cur_pred > random_pred)
            push!(path, v)
            break
        end
    end
    if (s != w && z != w)
        _backtrack_path!(s,z,v,path)
    end

end