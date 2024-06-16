function save_percolation_array(nn::String,x::Array{Float64})
    mkpath("percolation_states/")
    f::IOStream = open("percolation_states/" * nn *".txt", "w")
    for u in 1:lastindex(x)
        write(f, string(x[u]) * "\n")
    end
    close(f)
end

function read_percolation(file_name::String)
    @assert isfile(file_name) "The percolation values file does not exist"
    f = open(file_name, "r")
    p::Array{Float64} = []
    for line in eachline(f)
        l = parse(Float64, line)
        if l > 1.0 || l < 0.0
            println("Error, percolation state "*string(l)*" is not a valid state!")
            return nothing
        end 
        push!(p,parse(Float64, line))
    end
    close(f)
    return p
end

function save_centrality_values(file_name::String, centrality::Array{Float64})::Nothing
    f::IOStream = open(file_name, "w")
    for u in 1:lastindex(centrality)
        write(f, string(centrality[u]) * "\n")
    end
    close(f)
end

function append_centrality_values(file_name::String, centrality::Array{Float64})::Nothing
    f::IOStream = open(file_name, "a")
    for u in 1:lastindex(centrality)
        write(f, string(centrality[u]) * "\n")
    end
    close(f)
end


function read_time(file_name::String)
    @assert isfile(file_name) "The time value file does not exist: "*string(file_name)
    f = open(file_name, "r")
    t::Array{Float64} = []
    for line in eachline(f)
        split_line::Vector{String} = split(line, " ")
        push!(t,parse(Float64, split_line[1]))
    end
    close(f)
    if length(t) == 1
        return t[1]
    else
        return t
    end 
end

function read_sample_size(file_name::String)
    @assert isfile(file_name) "The time value file does not exist"
    f = open(file_name, "r")
    t::Array{Float64} = []
    for line in eachline(f)
        split_line::Vector{String} = split(line, " ")
        push!(t,parse(Float64, split_line[2]))
    end
    close(f)
    if length(t) == 1
        return t[1]
    else
        return t
    end 
end

function read_xi(file_name::String)
    @assert isfile(file_name) "The time value file does not exist"
    f = open(file_name, "r")
    t::Array{Float64} = []
    for line in eachline(f)
        split_line::Vector{String} = split(line, " ")
        push!(t,parse(Float64, split_line[3]))
    end
    close(f)
    if length(t) == 1
        return t[1]
    else
        return t
    end 
end
function save_time(nn::String,algo::String,tt::Float64)::nothing
    mkpath("times/" * nn * "/")
    f = open("times/" * nn *"/"*algo*".txt","a")
    write(f, string(round(tt; digits=4)))
    close(f)
end

function save_diameter(nn::String,diam::Int64,avg_dist::Float64,tt::Float64)
    mkpath("distance_metrics/")
    f = open("distance_metrics/" * nn *".txt","a")
    write(f, string(diam)*" "*string(avg_dist)*" "*string(round(tt; digits=4)))
    close(f)
end

function save_xi(nn::String,algo::String,xi::Float64)::nothing
    mkpath("xis/" * nn * "/")
    f = open("xis/" * nn *"/"*algo*".txt","a")
    write(f, string(round(xi; digits=4)))
    close(f)
end

function save_sample_size(nn::String,algo::String,ss::Int64)::nothing
    mkpath("sample_sizes/" * nn * "/")
    f = open("sample_sizes/" * nn *"/"*algo*".txt","a")
    write(f, string(ss))
    close(f)
end

function save_results_progressive_sampling(nn::String,cn::String, c::Array{Float64}, ss::Int64, t::Float64,starting_ss::Int64,xi::Float64 = -1.0)
    if (length(c) > 0)
        mkpath("scores/" * nn * "/")
        append_centrality_values("scores/" * nn * "/"*cn*"_"*string(starting_ss)*".txt", c)
        mkpath("times/" * nn * "/")
        f = open("times/" * nn * "/"*cn*"_"*string(starting_ss)*".txt", "a")
        write(f, string(round(t; digits=4)) * " " * string(ss) *" "*string(xi) *"\n")
        close(f)
    else
        mkpath("times/" * nn * "/")
        f = open("times/" * nn * "/"*cn*"_"*string(starting_ss)*".txt", "a")
        write(f, "-1.0 -1.0 -1.0 -1.0,-1.0")
        close(f)
    end
end

function save_results(nn::String, cn::String, c::Array{Float64}, t::Float64)
    if (length(c) > 0)
        mkpath("scores/" * nn * "/")
        save_centrality_values("scores/" * nn * "/" * cn * ".txt", c)
        mkpath("times/" * nn * "/")
        f = open("times/" * nn * "/time_" * cn * ".txt", "w")
        write(f, string(t))
        close(f)
    else
        mkpath("times/" * nn * "/")
        f = open("times/" * nn * "/time_" * cn * ".txt", "w")
        write(f, "-1.0\n")
        close(f)
    end
end

function read_centrality_values(file_name::String)::Array{Float64}
    @assert isfile(file_name) "The centrality value file does not exist "*file_name
    f::IOStream = open(file_name, "r")
    centrality::Array{Float64} = []
    value::Float64 = 0.0
    for l in eachline(f)
        value = parse(Float64, l)
        if (value < -0.1)
            println("ERROR. There are negative values with absolute big values")
            return Array{Float64}([])
        end
        if (value < 0)
            value = 0
        end
        push!(centrality, value)
    end
    close(f)
    return centrality
end

function read_distance_measures(nn::String)
    
    res = []

    file_name = "distance_metrics/"*nn*".txt"
    @assert isfile(file_name) "The diameter does not exists"
    f = open(file_name, "r")
    res = [parse(Float64,x) for x in split(readline(f)," ")]
    close(f)
  
    return res
end