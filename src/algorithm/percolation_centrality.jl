function percolation_centrality(g,percolation_states::Array{Float64},normalized::Bool = true)::Tuple{Array{Float64},Float64}
    @assert nv(g) == lastindex(percolation_states) "Vertex set and percolation array must have the same dimensions"
    println("Computing percolation centrality")
    flush(stdout)
    start_time::Float64 = time()
    n::Int64 = nv(g)
    delta::Array{Float64} = zeros(Float64,n)
    dist::Array{Int64} = zeros(Int64,n)
    ball::Array{Int16} = zeros(Int16,n)
    n_paths::Array{Int64} = zeros(Int64,n)
    coeff::Float64 = 0.0
    pred::Array{Array{Int64}} = Array{Array{Int64}}([[] for _ in 1:n])
    percolation::Array{Float64} = zeros(Float64,n)
    sum_percolations::Float64 = sum(percolation_states)
    q::Queue{Int64} = Queue{Int64}()
    s::Stack{Int64} = Stack{Int64}()
    for u in 1:n
        enqueue!(q,u)
        dist[u] = 0
        ball[u] = 1
        n_paths[u] = 1
        while length(q)!= 0
            w = dequeue!(q)
            for v in outneighbors(g,w)
                if ball[v] == 0
                    dist[v] = dist[w]+1
                    ball[v] = 1
                    enqueue!(q,v)
                    push!(s,v)
                end
                if (ball[v] == 1 && dist[v] == dist[w]+1)
                    n_paths[v] += n_paths[w]
                    push!(pred[v],w)
                end
            end
        end
        while length(s)!= 0
            w = pop!(s)
            coeff = (1 + delta[w]) / n_paths[w]
            for p in pred[w]
                delta[p] += n_paths[p] * coeff
            end
            if w!= u
                percolation[w] += delta[w] * percolation_states[u]/(sum_percolations-percolation_states[w])
            end
        end
        for v in 1:n
            dist[v] = 0
            ball[v] = 0
            n_paths[v] = 0
            delta[v] = 0
            pred[v] = Array{Int64}([])
        end
        # This to be 100% sure
        q = Queue{Int64}()
        s = Stack{Int64}()
        #println("Processed "*string(u)*"/"*string(n)*" in time "*string(time()-start_time))
    end
    if normalized
        percolation = percolation .* [1/(n*(n-1))]
    end
    finish_time::Float64 = time()-start_time
    println("Percolation centrality computed in "*string(finish_time))
    flush(stdout)
    return percolation,finish_time
end