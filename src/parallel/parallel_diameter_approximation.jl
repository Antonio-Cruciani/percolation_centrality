function parallel_random_bfs(g,sample_size::Int64)::Tuple{Int64,Float64}
    start_time::Float64 = time()
    n::Int64 = nv(g)
    ntasks = nthreads()
    diam_lb::Array{Array{Int64}} = [[0] for _ in 1:ntasks]
    reachable_pairs::Array{Array{Int64}} = [[0] for _ in 1:ntasks]
    lb::Int64 = 0
    exact::Bool = false
    @info("Running (parallel) Random BFS algorithm")
    # If sample size set to 0, compute the exact diameter
    if sample_size == 0
        sample_size = n
        exact = true
        @info("Computing Exact Diameter")
        flush(stderr)
    else
        @info("Computing Approximated Diameter")
        flush(stderr)
    end
    task_size = cld(sample_size, ntasks)
    vs_active = [i for i in 1:sample_size]
    s::Int64 = -1
    @sync for (t, task_range) in enumerate(Iterators.partition(1:sample_size, task_size))
        Threads.@spawn for u in @view(vs_active[task_range])
            if exact
                s = u
            else
                s = sample(1:n)
            end
            _bfs!(s,g,diam_lb[t],reachable_pairs[t])
        end
    end
    rp::Int64 = 0
    for t in 1:ntasks
        if lb < diam_lb[t][1]
            lb = diam_lb[t][1]
        end
        rp += reachable_pairs[t][1]
    end
    alpha::Float64 = n*rp/(sample_size*n*(n-1))
    @info("Connectivity rate "*string(round(alpha;digits = 5)))
    flush(stderr)
    return lb,time()-start_time
end


function _bfs!(s::Int64,g,diam::Array{Int64},reachable_pairs::Array{Int64})
    n::Int64 = nv(g)
    q::Queue{Int64} = Queue{Int64}()
    dist::Array{Int64} = zeros(Int64,n)
    # 0 not discovered, 1 discovered

    ball::Array{Int16} = zeros(Int16,n)
    w::Int64 = -1
    enqueue!(q,s)
    dist[s] = 0
    ball[s] = 1
    while length(q)!= 0
        w = dequeue!(q)
        for v in outneighbors(g,w)
            if ball[v] == 0
                dist[v] = dist[w] +1
                reachable_pairs[1] += 1
                if dist[v] > diam[1]
                    diam[1] = dist[v]
                end
                ball[v] = 1
                enqueue!(q,v)
            end
        end
    end

    return nothing
end


function parallel_random_bfs_rho(g,sample_size::Int64)::Tuple{Int64,Float64,Float64}
    start_time::Float64 = time()
    n::Int64 = nv(g)
    ntasks = nthreads()
    diam_lb::Array{Array{Int64}} = [[0] for _ in 1:ntasks]
    reachable_pairs::Array{Array{Int64}} = [[0] for _ in 1:ntasks]
    lb::Int64 = 0
    rho::Array{Array{Int64}} = [zeros(Int64,n+1) for _ in 1:ntasks]
    final_rho::Array{Float64} = zeros(Float64,n+1)
    exact::Bool = false
    @info("Running (parallel) Random BFS algorithm")
    # If sample size set to 0, compute the exact diameter
    if sample_size == 0
        sample_size = n
        exact = true
        @info("Computing Exact Diameter")
        flush(stderr)
    else
        @info("Computing Approximated Diameter")
        flush(stderr)
    end
    task_size = cld(sample_size, ntasks)
    vs_active = [i for i in 1:sample_size]
    s::Int64 = -1
    @sync for (t, task_range) in enumerate(Iterators.partition(1:sample_size, task_size))
        Threads.@spawn for u in @view(vs_active[task_range])
            if exact
                s = u
            else
                s = sample(1:n)
            end
            _bfs_rho!(s,g,diam_lb[t],rho[t],reachable_pairs[t])
        end
    end
    rp::Int64 = 0
    for t in 1:ntasks
        if lb < diam_lb[t][1]
            lb = diam_lb[t][1]
        end
        rp += reachable_pairs[t][1]
    end
    final_rho = reduce(+,rho)
    avg_dist::Float64 = 0.0
    temporal_hop_table = zeros(lb+1)
    accum = 0
    for h in 1:(lb+1)
        accum += final_rho[h]
        temporal_hop_table[h] = n * accum/sample_size
    end
    #=
    for d in 1:lb
        avg_dist += (d-1)*final_rho[d]
    end
    =#
    avg_dist = average_distance(temporal_hop_table)
    #avg_dist = avg_dist /(n*rp/sample_size)
    alpha::Float64 = n*rp/(sample_size*n*(n-1))
    @info("Connectivity rate "*string(round(alpha;digits = 5)))
    flush(stderr)
    return lb,avg_dist,time()-start_time
end



function _bfs_rho!(s::Int64,g,diam::Array{Int64},rho::Array{Int64},reachable_pairs::Array{Int64})
    n::Int64 = nv(g)
    q::Queue{Int64} = Queue{Int64}()
    dist::Array{Int64} = zeros(Int64,n)
    # 0 not discovered, 1 discovered

    ball::Array{Int16} = zeros(Int16,n)
    w::Int64 = -1
    enqueue!(q,s)
    dist[s] = 0
    ball[s] = 1
    rho[1] +=1
    while length(q)!= 0
        w = dequeue!(q)
        for v in outneighbors(g,w)
            if ball[v] == 0
                dist[v] = dist[w] +1
                reachable_pairs[1] += 1
                rho[dist[v]] += 1
                if dist[v] > diam[1]
                    diam[1] = dist[v]
                end
                ball[v] = 1
                enqueue!(q,v)
            end
        end
    end

    return nothing
end

function average_distance(ht::Array{Float64})::Float64
    if length(ht) == 0
        return 0.0
    end
    distance::Array{Float64} = distance_function(ht)
    m::Float64 = 0.0
    for i in 1:lastindex(distance)
        m += (distance[i] * (i-1))
    end
    return m/ht[length(ht)]
end


function distance_function(ht::Array{Float64})::Array{Float64}
    table::Array{Float64} = copy(ht)
    for i in lastindex(table):-1:2
        table[i] -= table[i-1]
    end
    return table
end