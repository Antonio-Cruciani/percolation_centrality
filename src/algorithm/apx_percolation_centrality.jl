
function estimate_percolation_centrality(g,percolation_states::Array{Float64},epsilon::Float64,delta::Float64, alpha_sampling::Float64 = 1.0)
    n::Int64 = nv(g)
    m::Int64 = ne(g)
    mc_trials::Int64 = 25
    directed::Bool = is_directed(g)
    sg::static_graph = static_graph(adjacency_list(g),incidency_list(g))
    #sg.adjacency = adjacency_list(g)
    #sg.in_list = incidency_list(g)
    q::Array{Int64} = zeros(Int64,n)
    ball::Array{Int16} = zeros(Int16,n)
    n_paths::Array{Int64} = zeros(Int64,n)
    dist::Array{Int64} = zeros(Int64,n)
    pred::Array{Array{Int64}} = [Array{Int64}([]) for _ in 1:n]
    num_paths::Array{Int64} = zeros(Int64,n)
    percolation_centrality::Array{Float64} = zeros(Float64,n)
    wimpy_variance::Array{Float64} = zeros(Float64,n)
    shortest_path_length::Array{Int64} = zeros(Int64,n+1)
    percolated_path_length::Array{Float64} = zeros(Float64,n+1)
    mcrade::Array{Float64} = zeros(Float64,(n+1)*mc_trials)
    tmp_perc_states::Array{Float64} = copy(percolation_states)
    percolation_data::Tuple{Float64,Array{Float64}} = percolation_differences(sort(tmp_perc_states),n)
    betweenness::Array{Float64} = zeros(Float64,n)

    #= Test
    test_couples = []
    for i in 1:n
        for j in i:n
            if i!=j
                push!(test_couples,(i,j))
            end
        end
    end
    println(" number of couples ",length(test_couples))
    #println(test_couples)
    for tc in test_couples
        _random_path!(tc[1],tc[2],betweenness,sg,n,q,ball,n_paths,dist,pred,num_paths,percolation_centrality,wimpy_variance,percolation_states,percolation_data,shortest_path_length,percolated_path_length,mcrade,mc_trials,alpha_sampling)
    end
    #println(betweenness)
    
    return (betweenness)
    =#
    _random_path!(sg,n,q,ball,n_paths,dist,pred,num_paths,percolation_centrality,wimpy_variance,percolation_states,percolation_data,shortest_path_length,percolated_path_length,mcrade,mc_trials,betweenness,alpha_sampling)

    println("---------------")
    println(percolation_centrality)
    println(wimpy_variance)
    println("-------------------------")
    println(shortest_path_length)
    println(percolated_path_length)

end