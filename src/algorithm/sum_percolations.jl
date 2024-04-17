function percolation_differences(percolation_states::Array{Float64},n::Int64)::Tuple{Float64,Array{Float64}}
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