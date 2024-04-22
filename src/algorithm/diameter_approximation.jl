function random_bfs(g,sample_size::Int64)::Tuple{Int64,Float64}
    start_time::Float64 = time()
    n::Int64 = nv(g)
    diam_lb::Int64 = 0
    q::Queue{Int64} = Queue{Int64}()
    dist::Array{Int64} = zeros(Int64,n)
    # 0 not discovered, 1 discovered
    ball::Array{Int16} = zeros(Int16,n)
    s::Int64 = -1
    w::Int64 = -1
    reachable_pairs::Int64 = 0
    # If sample size set to 0, compute the exact diameter
    if sample_size == 0
        sample_size = n
        @info("Computing Exact Diameter")
        flush(stdout)
    end
    for i in 1:sample_size
        if sample_size == n
            s = i
        else
            s = sample(1:n)
        end
        enqueue!(q,s)
        dist[s] = 0
        ball[s] = 1
        while length(q)!= 0
            w = dequeue!(q)
            for v in outneighbors(g,w)
                if ball[v] == 0
                    dist[v] = dist[w] +1
                    reachable_pairs += 1
                    if dist[v] > diam_lb
                        diam_lb = dist[v]
                    end
                    ball[v] = 1
                    enqueue!(q,v)
                end
            end
        end
        for v in 1:n
            dist[v] = 0
            ball[v] = 0
        end
    end
    alpha::Float64 = n*reachable_pairs/(sample_size*n*(n-1))
    @info("Connectivity rate "*string(round(alpha;digits = 5)))
    flush(stdout)
    return diam_lb,time()-start_time
end