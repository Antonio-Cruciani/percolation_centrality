function parallel_percolation_centrality(g,percolation_states::Array{Float64},normalized::Bool = true)::Tuple{Array{Float64},Float64}
    @assert nv(g) == lastindex(percolation_states) "Vertex set and percolation array must have the same dimensions"
    @info("----------------------------------------| Stats |--------------------------------------------------")
    @info("Analyzing graph")
    @info("Number of nodes "*string(nv(g)))
    @info("Number of edges "*string(ne(g)))
    @info("Directed ? "*string(is_directed(g)))
    @info("Maximum Percolated state "*string(maximum(percolation_states)))
    @info("Minimum Percolated state "*string(minimum(percolation_states)))
    @info("Average Percolated state "*string(mean(percolation_states))*" std "*string(std(percolation_states)))
    @info("Using "*string(nthreads())* " Threads")
    @info("---------------------------------------------------------------------------------------------------")
    @info("Computing percolation centrality")
    flush(stderr)
    start_time::Float64 = time()
    n::Int64 = nv(g)
    #percent_status::Array{Tuple{Float64,Int64,String}} = [(0.25,trunc(Int64,n*0.25),"Analyzed 25% of the graph"),(0.50,trunc(Int64,n*0.50),"Analyzed 50% of the graph"),(0.75,trunc(Int64,n*0.75),"Analyzed 75% of the graph"),(0.9,trunc(Int64,n*0.9),"Analyzed 90% of the graph")]
    #percent_index::Int16 = 1
    vs_active = [i for i in 1:n]
    d, r = divrem(n, nthreads())
    ntasks = d == 0 ? r : nthreads()
    task_size = cld(n, ntasks)
    percolation::Array{Array{Float64}} = [zeros(Float64,n) for _ in 1:ntasks]
    final_percolation::Array{Float64} = zeros(Float64,n)
    sum_percolations::Float64 = sum(percolation_states)
    verbose_step::Int64 = trunc(Int64,floor(n*0.25))
    processed_so_far::Int64 = 0
    @sync for (t, task_range) in enumerate(Iterators.partition(1:n, task_size))
        Threads.@spawn for s in @view(vs_active[task_range])
            _per_cent_accumulate!(s,g,percolation_states,sum_percolations,percolation[t])
            if (Sys.free_memory() / Sys.total_memory() < 0.1)
                clean_gc()
                sleep(0.01)
            end
            processed_so_far +=1
            if (verbose_step > 0 && processed_so_far % verbose_step == 0)
                finish_partial::String = string(round(time() - start_time; digits=4))
                time_to_finish::String = string(round((n*(time() - start_time) / processed_so_far )-(time() - start_time) ; digits=4))
                @info("Percolation Centrality. Processed " * string(processed_so_far) * "/" * string(n) * " nodes in " * finish_partial * " seconds | Est. remaining time : "*time_to_finish)
                flush(stderr)
            end
            #=
            if percent_index <= length(percent_status)
                if trunc(Int64,s*percent_status[percent_index][1]) == percent_status[percent_index][2]
                    finish_partial = string(round(time()-start_time; digits=4))
                    est =string(round( (time()-start_time )* n/s;digits = 4))
                    @info(percent_status[percent_index][3]*" in "*finish_partial*" seconds | Estimated remaining time "*est* " seconds")
                    percent_index+=1
                end
            end
            =#
        end
    end
    final_percolation = reduce(+,percolation)
    if normalized
        final_percolation = final_percolation .* [1/(n*(n-1))]
    end
    finish_time::Float64 = time()-start_time
    @info("Percolation centrality computed in "*string(finish_time))
    flush(stderr)
    return final_percolation,finish_time
end


function _per_cent_accumulate!(u::Int64,g,percolation_states::Array{Float64},sum_percolations::Float64,percolation::Array{Float64})
    coeff::Float64 = 0.0
    n::Int64 = nv(g)
    w::Int64 = -1
    delta::Array{Float64} = zeros(Float64,n)
    dist::Array{Int64} = zeros(Int64,n)
    ball::Array{Int16} = zeros(Int16,n)
    n_paths::Array{UInt128} = zeros(UInt128,n)
    pred::Array{Array{Int64}} = Array{Array{Int64}}([[] for _ in 1:n])
    q::Queue{Int64} = Queue{Int64}()
    s::Stack{Int64} = Stack{Int64}()
    
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
        
        #println("Processed "*string(u)*"/"*string(n)*" in time "*string(time()-start_time))
    
    
    return nothing

end

