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


function _random_path!(sg::static_graph,n::Int64,q::Array{Int64},ball::Array{Int16},n_paths::Array{Int64},dist::Array{Int64},pred::Array{Int64,Array{Int64}},q::Array{Int64})

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
    sp_edges::Array{Tuple{Int64,Int64}} = Array{Tuple{Int64,Int64}}([])
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
    end

    for u in 1:end_q
        ball[u] = @unvisited
        dist[u] = 0
        pred[u] = Array{Int64}([])
        # Check this npaths and q, they should be set back to 0 once finished
        n_paths[u] = 0
        q[u] = 0
    end
    return nothing
end