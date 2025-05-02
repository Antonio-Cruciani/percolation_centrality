
struct static_graph
    adjacency::Array{Array{Int64}}
    incidency::Array{Array{Int64}}
    degrees_adj::Array{Int64}
    degrees_idj::Array{Int64}
    function static_graph(adj::Array{Array{Int64}},idj::Array{Array{Int64}})
        return new(adj,idj,[lastindex(adj[u]) for u in 1:lastindex(adj)],[lastindex(idj[u]) for u in 1:lastindex(idj)])
    end
end


function load_graph(file_name::String,directed::Bool = true,sep::String = " ")
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
                #push!(file_id, u)
                current_node_id = current_node_id + 1
            end
            if (!haskey(file_id_to_graph_id, v))
                file_id_to_graph_id[v] = current_node_id
                #push!(file_id, v)
                current_node_id = current_node_id + 1
            end 
            
            if directed
                push!(edges,(u,v))
            else
                push!(edges,(min(u,v),max(u,v)))
            end
            push!(nodes,u)
            push!(nodes,v)
        end
    end
    close(f)

    g = SimpleGraph()
    #n_vertices::Int64 = maximum(collect(keys(file_id_to_graph_id)))
    n_vertices::Int64 = maximum(collect(values(file_id_to_graph_id)))
    if directed
        g = SimpleDiGraph(n_vertices)
    else
        g = SimpleGraph(n_vertices)
    end
    for e in collect(edges)
        res = add_edge!(g,file_id_to_graph_id[e[1]],file_id_to_graph_id[e[2]])
        if !res 
			@info("Edge not added! "*string(e[1])*" , "*string(e[2]))
			@info("Edge not added! "*string(file_id_to_graph_id[e[1]])*" , "*string(file_id_to_graph_id[e[2]]))
		end
    end
    loading_time = time() - start_time
    @info("Loaded Graph ")
	@info("#Nodes "*string(nv(g)))
	@info("#Edges "*string(ne(g)))
	@info("Directed ? "*string(directed))
    @info("Loading time "*string(loading_time)*" seconds")
    flush(stderr)
	return g
end


function print_stats(g; graph_name="anonymous")
    @info("====================================================")
    @info("Network: " * graph_name)
    @info("====================================================")
    @info("Number of nodes " * string(nv(g)))
    @info("Number of edges " * string(length(ne(g))))
    @info("Is the graph directed? " * string(is_directed(g)))
    @info("====================================================")
    flush(stderr)
end


function adjacency_list(g)::Array{Array{Int64}}
    adj_list::Array{Array{Int64}} = Array{Array{Int64}}([])
    for u in 1:nv(g)
        push!(adj_list,outneighbors(g,u))
    end
    return adj_list
end

function incidency_list(g)::Array{Array{Int64}}
    in_list::Array{Array{Int64}} = Array{Array{Int64}}([])
    for u in 1:nv(g)
        push!(in_list,inneighbors(g,u))
    end
    return in_list
end


function save_graph(g,nn,sep::String = " ")
    mkpath("components/")
    f = open("components/" * nn * ".txt", "w")
    for e in collect(edges(g))
        write(f,string(src(e))*sep*string(dst(e))*"\n")
    end
    close(f)
end


function load_and_normalize_edgelist(path::String)
    # Containers
    edges = Tuple{String, String}[]
    label_to_id = Dict{String, Int}()
    id_to_label = String[]
    
    # Read file
    open(path, "r") do file
        for line in eachline(file)
            if isempty(strip(line)) || startswith(line, "#")
                continue  # skip comments and empty lines
            end
            tokens = split(strip(line))
            if length(tokens) < 2
                continue  # skip malformed lines
            end
            push!(edges, (tokens[1], tokens[2]))
        end
    end

    # Map labels to consecutive IDs
    next_id = 0
    for (u, v) in edges
        for node in (u, v)
            if !haskey(label_to_id, node)
                label_to_id[node] = next_id
                push!(id_to_label, node)
                next_id += 1
            end
        end
    end

    # Normalize edges
    normalized_edges = [(label_to_id[u], label_to_id[v]) for (u, v) in edges]

    return normalized_edges, label_to_id, id_to_label
end


function save_edge_list(el,nn,sep = " ")
    mkpath("normalized/")
    f = open("normalized/" * nn * ".txt", "w")
    for e in el
        write(f,string(e[1])*sep*string(e[2])*"\n")
    end
    close(f)

end