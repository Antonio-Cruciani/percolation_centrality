function parallel_random_bfs(g,sample_size::Int64)::Tuple{Int64,Float64}
    start_time::Float64 = time()
    n::Int64 = nv(g)
    ntasks = nthreads()
    diam_lb::Array{Array{Int64}} = [[0] for _ in 1:ntasks]
    lb::Int64 = 0
    exact::Bool = false
    @info("Running (parallel) Random BFS algorithm")
    # If sample size set to 0, compute the exact diameter
    if sample_size == 0
        sample_size = n
        exact = true
        @info("Computing Exact Diameter")
        flush(stdout)
    else
        @info("Computing Approximated Diameter")
        flush(stdout)
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
            _bfs!(s,g,diam_lb[t])
        end
    end

    for t in 1:ntasks
        if lb < diam_lb[t][1]
            lb = diam_lb[t][1]
        end
    end

    
    return lb,time()-start_time
end


function _bfs!(s::Int64,g,diam::Array{Int64})
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