
struct static_graph
    adjacency::Array{Int64,Array{Int64}}
    incidency::Array{Int64,Array{Int64}}
    function static_graph(n::Int64)
        return new(Array{Int64,Array{Int64}}(undef,n),Array{Int64,Array{Int64}}(undef,n))
    end
end


function load_graph(file_name::String,directed::Boolean = true,sep::String = " ")
    @assert isfile(file_name) "The edge list file " * file_name * " does not exist"
    start_time = time()
    file_id_to_graph_id::Dict{Int64,Int64} = Dict{Int64,Int64}()
    current_node_id::Int64 = 1
    f::IOStream = open(file_name, "r")
	edges::Set{Tuple{Int64,Int64}} = Set{Tuple{Int64,Int64}}()
    nodes::Set{Int64} = Set{Int64}()
    for line in eachline(f)
        split_line::Vector{String} = split(line, sep)
        @assert length(split_line) == 2 "Bad line format: " * line
        u = parse(Int64,split_line[1])
		v = parse(Int64,split_line[2])
        if u != v
            if (!haskey(file_id_to_graph_id, u))
                file_id_to_graph_id[u] = current_node_id
                push!(file_id, u)
                current_node_id = current_node_id + 1
            end
            if (!haskey(file_id_to_graph_id, v))
                file_id_to_graph_id[v] = current_node_id
                push!(file_id, v)
                current_node_id = current_node_id + 1
            end 
        end
        if directed
            push!(edges,(u,v))
        else
            push!(edges,(min(u,v),max(u,v)))
        end
        push!(nodes,u)
        push!(nodes,v)
    end
    close(f)

    g = SimpleGraph()
    n_vertices::Int64 = maximum(collect(keys(file_id_to_graph_id)))
    if directed
        g = SimpleDiGraph(n_vertices)
    else
        g = SimpleGraph(n_vertices)
    end
    for e in collect(edges)
        res = add_edge!(g,file_id_to_graph_id[e[1]],file_id_to_graph_id[e[2]])
        if !res 
			println("Edge not added! "*string(e[1])*" , "*string(e[2]))
			println("Edge not added! "*string(file_id_to_graph_id[e[1]])*" , "*string(file_id_to_graph_id[e[2]]))
		end
    end
    loading_time = time() - start_time
    println("Loaded Graph ")
	println("#Nodes "*string(nv(g)))
	println("#Edges "*string(ne(g)))
	println("Directed ? "*string(directed))
    println("Loading time "*string(loading_time)*" seconds")
	return g
end


function print_stats(g; graph_name="anonymous")
    println("====================================================")
    println("Network: " * graph_name)
    println("====================================================")
    println("Number of nodes " * string(nv(g))))
    println("Number of edges " * string(length(ne(g))))
    println("Is the graph directed? " * string(is_directed(g)))
    println("====================================================")
end


function adjacency_list(g)::Array{Int64,Array{Int64}}
    adj_list::Array{Int64,Array{Int64}} = Array{Int64,Array{Int64}}([ for _ in 1:nv(g)])
    for u in 1:nv(g)
        push!(adj_list[u],outneighbors(g,u))
    end
    return adj_list
end

function incidency_list(g)::Array{Int64,Array{Int64}}
    in_list::Array{Int64,Array{Int64}} = Array{Int64,Array{Int64}}([ for _ in 1:nv(g)])
    for u in 1:nv(g)
        push!(in_list[u],inneighbors(g,u))
    end
    return in_list
end