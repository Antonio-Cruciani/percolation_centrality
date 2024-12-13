@inline function compute_d_max(n::Int64,X::Array{Float64})::Tuple{Array{Float64},Float64}
    start_time::Float64 = time()
    @info("Computing d_v = max_{sv} κ/̃κ")
    tmp_perc_states::Dict{Int64,Float64} = Dict(v => X[v] for v in 1:n)
    sorted_dict = OrderedDict(sort(collect(tmp_perc_states), by = kv -> kv[2]))
    percolation_data::Tuple{Float64,Dict{Int64,Float64}} = percolation_differences(sorted_dict,n)



    d::Array{Float64} = ones(Float64,n)
    

    # Extract keys and values from the sorted dictionary
    keys_sorted = collect(keys(X))
    values_sorted = collect(values(X))

    Y = Dict{Int, Float64}()

    # Prefix sum of values
    prefix_sum = zeros(Float64, n)
    prefix_sum[1] = values_sorted[1]
    for i in 2:n
        prefix_sum[i] = prefix_sum[i - 1] + values_sorted[i]
    end

    total_sum = prefix_sum[n]

    for i in 1:n
        # Sum of differences for values_sorted[i]
        diff_sum = values_sorted[i] * (i - 1) - (i > 1 ? prefix_sum[i - 1] : 0)
        diff_sum += (total_sum - prefix_sum[i]) - values_sorted[i] * (n - i)
        Y[keys_sorted[i]] = diff_sum
    end
 
    for v in 1:n
        d[v] += Y[v]/percolation_data[2][v]
    end
    



    end_time::Float64 = time()-start_time
    @info("minimum d_v = "*string(minimum(d))*" maximum d_v = "*string(maximum(d))*" computed in "*string(end_time)*" seconds")
    return d,end_time
end










#=


@inline function approx_d_max(n::Int64,X::Array{Float64},sample_size::Int64)::Tuple{Array{Float64},Float64}
    start_time::Float64 = time()
    @info("Approximatin computing d_v = max_{sv} κ/̃κ")
    @info("Sample size = "*string(sample_size))
    tmp_perc_states::Dict{Int64,Float64} = Dict(v => X[v] for v in 1:n)
    sorted_dict = OrderedDict(sort(collect(tmp_perc_states), by = kv -> kv[2]))
    percolation_data::Tuple{Float64,Dict{Int64,Float64}} = percolation_differences(sorted_dict,n)
    ntasks = nthreads()
    d::Array{Array{Float64}} = [zeros(Float64,n) for _ in 1:ntasks]
    task_size = cld(sample_size, ntasks)
    vs_active = [i for i in 1:sample_size]
    @info("Using "*string(ntasks)*" threads")
    @sync for (t, task_range) in enumerate(Iterators.partition(1:sample_size, task_size))
        Threads.@spawn for _ in @view(vs_active[task_range])
            s = sample(1:n)
            z = s
            while s == z
                z = sample(1:n)
            end
            if ramp(X[s],X[z]) > 0
                for v in 1:n
                    til_kappa::Float64 = ramp(X[s],X[z])/(percolation_data[1] - (ramp(X[s],X[z]) + ramp(X[z],X[s])))
                    kappa::Float64 = ramp(X[s],X[z])/percolation_data[2][v]
                    tmp_dv::Float64 =  kappa/til_kappa
                    d[t][v] = tmp_dv 
                end
            end
        end
    end
    @info("Sampling completed, reduction phase")
    final_d_v::Array{Float64} = zeros(Float64,n)
    task_size = cld(n, ntasks)
    vs_active = [i for i in 1:n]
    @sync for (t, task_range) in enumerate(Iterators.partition(1:n, task_size))
        Threads.@spawn for v in @view(vs_active[task_range])
            for i in 1:ntasks
                final_d_v[v] = max(final_d_v[v],d[i][v])
            end
        end
    end
    end_time::Float64 = time()-start_time
    @info("minimum d_v = "*string(minimum(final_d_v))*" maximum d_v = "*string(maximum(final_d_v))*" computed in "*string(end_time)*" seconds")
    return final_d_v,end_time
end

=#


@inline function partition_percolation_states(percolation_states::Array{Float64},alpha::Float64 = 2.0)::Dict{Int64,Array{Int64}}
    @info("Partitioning the percolation states")
    partition::Dict{Int64,Array{Int64}} = Dict{Float64,Array{Int64}}()
    n::Int64 = lastindex(percolation_states)
    #k::Int64 = 0
    #range::Tuple{Float64,Float64} = Tuple{Float64,Float64}(0.,0.)
    #eps::Float64 = 0.000000001
    #partition_index::Int64 = 0
    state_partition_index::Int64 = 0
    min_inv_state::Float64 = 0
    partition_number::Int64 = 0
    for i in 1:n
        min_inv_state = min(1. /percolation_states[i],n)
        state_partition_index = trunc(Int64, log(alpha,min_inv_state)+1)

        if haskey(partition,state_partition_index)
            push!(partition[state_partition_index],i)
        else
            partition[state_partition_index] = [i]
            partition_number += 1
        end
    end
    @info("Number of partitions "*string(partition_number))
    return partition
end



