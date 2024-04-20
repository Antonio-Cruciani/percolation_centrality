function supremum_deviation(x::Array{Float64},y::Array{Float64})::Float64
    @assert lastindex(x) == lastindex(y) "x and y must have the same length"
    sd::Float64 = 0.0
    for i in 1:lastindex(x)
        if sd < abs(x[i] - y[i])
            sd = abs(x[i] - y[i])
        end
    end
    return sd
end