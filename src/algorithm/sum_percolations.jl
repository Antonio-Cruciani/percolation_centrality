function percolation_differences_dep(percolation_states::Array{Float64},n::Int64)::Tuple{Float64,Array{Float64}}
    summation::Float64 = 0
    minus_sum::Array{Float64} = zeros(Float64,n)
    svp::Array{Float64} = zeros(Float64,n+1)
    for i in 2:n
        svp[i] = svp[i-1] + percolation_states[i-1]
        summation += (i-1)*percolation_states[i] - svp[i]
    end
    svp[n+1] = svp[n] + percolation_states[n]
    for i in 1:n
        minus_sum[i] = summation - percolation_states[i] * (2*i-n-2) - svp[n+1] + 2*svp[i]
    end
    return summation,minus_sum
end

function ramp(x::Float64,y::Float64)::Float64
    if x-y>0
        return (x-y)
    else
        return 0
    end
end


function random_percolations(n::Int64)::Array{Float64}
    return Array{Float64}([rand() for _ in 1:n])
end

function custom_percolation_randomized(n::Int64,rnd_size::Int64,eps_size::Int64,zero_size::Int64,one_size::Int64,epsilon::Float64)::Array{Float64}
    @assert n == rnd_size + eps_size + zero_size + one_size "the sum of the ranges must sum to n"
    percolations::Array{Float64} = zeros(Float64,n)
    indexes::Dict{Int64,Int16} = Dict(i => 0 for i in 1:n)
    if rnd_size > 0
        rnd_indexes::Array{Int64}  = sample(collect(keys(indexes)),rnd_size,replace = false)
        for e in rnd_indexes
            percolations[e] = rand()
            delete!(indexes, e)
        end
    end
    if eps_size > 0
        eps_indexes::Array{Int64} = sample(collect(keys(indexes)),eps_size,replace = false)
        for e in eps_indexes
            percolations[e] = epsilon
            delete!(indexes, e)
        end
    end
    # Can be removed 
    if zero_size > 0
        zero_indexes::Array{Int64} = sample(collect(keys(indexes)),zero_size,replace = false)
        for e in zero_indexes
            percolations[e] = 0.0
            delete!(indexes, e)
        end
    end
    if one_size > 0
        one_indexes::Array{Int64} = sample(collect(keys(indexes)),one_size,replace = false)
        for e in one_indexes
            percolations[e] = 1.0
            delete!(indexes, e)
        end
    end
    return percolations
end


function custom_percolations(n::Int64,target,ranges::Array{Int64})::Array{Float64}
    @assert sum(ranges) == n "the sum of the ranges must sum to n"
    percolations::Array{Float64} = zeros(Float64,n)
    start = 1
    end_v = ranges[1]
    k = 1
    for i in 1:lastindex(ranges)
        if i > 1
            start = k
            end_v = k + ranges[i] -1
        end
        for j in start:end_v
            if target[i] == "r"
                percolations[j] = rand()
            else
                percolations[j] = target[i]
            end
            k+=1
        end
    end
    return percolations
end

function percolation_differences(percolation_states,n::Int64)::Tuple{Float64,Dict{Int64,Float64}}
    sum::Float64 = 0.0
    minus_sum::Dict{Int64,Float64} = Dict(v => 0.0 for v in 1:n)
    svp::Array{Float64} = zeros(Float64, n + 1)

    j::Int64 = 0
    k::Float64 = 0.0
    for (key, value) in percolation_states
        if j > 0
            svp[j + 1] = svp[j] + k
            sum += j * value - svp[j + 1]
        end
        k = value
        j += 1
    end
    svp[end] = svp[n] + k

    j = 0
    for (key, value) in percolation_states
        minus_sum[key] = sum - value * (2 * j - n) - svp[end] + 2 * svp[j + 1]
        j += 1
    end

    return sum,minus_sum
end