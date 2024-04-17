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
    # If sample size set to 0, compute the exact diameter
    if sample_size == 0
        sample_size = n
    end
    for _ in 1:sample_size
        s = sample(1:n)
        enqueue!(q,s)
        dist[s] = 0
        ball[s] = 1
        while length(q)!= 0
            w = dequeue!(q)
            for v in outneighbors(g,w)
                if ball[v] == 0
                    dist[v] = dist[w] +1
                    if dist[v] > diam_lb
                        diam_lb = dist[v]
                    end
                end
            end
        end
        for v in 1:n
            dist[v] = 0
            ball[v] = 0
        end
    end
    return diam_lb,time()-start_time
end