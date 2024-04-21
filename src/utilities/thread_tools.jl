function _reduce_data!(u::Int64,tn::Int64,src_x::Array{Array{Float64}},src_y::Array{Array{Float64}},src_z::Array{Array{Int64}},dst_x::Array{Float64},dst_y::Array{Float64},dst_z::Array{Int64})

    for t in 1:tn
        dst_x[u]+=src_x[t][u]
        dst_y[u]+=src_y[t][u]
        dst_z[u]+=src_z[t][u]
    end

    return nothing
end


function _reduce_data_a!(u::Int64,tn::Int64,src_x::Array{Array{Float64}},dst_x::Array{Float64})

    for t in 1:tn
        dst_x[u]+=src_x[t][u]
    end

    return nothing
end