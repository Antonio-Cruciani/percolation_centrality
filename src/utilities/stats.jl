
#using Conda
#using PyCall
#Conda.add("scipy")

function supremum_deviation(x,y)::Float64
    @assert lastindex(x) == lastindex(y) "x and y must have the same length"
    sd::Float64 = 0.0
    for i in 1:lastindex(x)
        if sd < abs(x[i] - y[i])
            sd = abs(x[i] - y[i])
        end
    end
    return sd
end



#=
py"""
from scipy import stats
def ktau(x, y):
    return stats.kendalltau(x, y)
def wktau(x, y):
    return stats.weightedtau(x, y)
def spear(x,y):
    return stats.spearmanr(x, y)
def pear(x,y):
    return stats.pearsonr(x, y)
"""

function compute_correlations(x, y, verbose::Bool)
    sp = py"spear"(x, y)
    if (verbose)
        log("    Spearman computed")
    end
    kt = py"ktau"(x, y)
    if (verbose)
        log("    Kendall tau computed")
    end
    wkt = py"wktau"(x, y)
    if (verbose)
        log("    Weighted Kendall tau computed")
    end
    p =  py"pear"(x, y)
    if (verbose)
        log("    Pearson computed")
    end
    return sp, kt, wkt,p
end

function normalize_centrality(x::Array{Float64})::Array{Float64}
    n::Int64 = length(x)
    for i in 1:n
        x[i] = x[i]/(n*(n-1))
    end
    return x
end



=#

function mean_sqaured_error(x,y)
    MSE = 0
    for i in 1:lastindex(x)
        MSE += (x[i]-y[i])^2
    end
    return (1/lastindex(x)) * MSE
end


function normalize_centrality(x::Array{Float64})::Array{Float64}
    n::Int64 = length(x)
    for i in 1:n
        x[i] = x[i]/(n*(n-1))
    end
    return x
end


function get_errors(starting_sample::Int64,target_epsilon::Float64,datasets::Array{String},algo::String ="cmcera",trials::Int64 = 10,directed::Bool = true)
    results = []
    
    for graph in datasets
        g = load_graph("graphs/"*graph, directed,"\t")
        gn = split(graph,".txt")[1]
        #println("Analyzing "*gn*" ε ",target_epsilon, " Algorithm ",algo," Progressive Sampler ",prog_sampler, " Path Optimality ",method, " Upper bound samples via ",upper_bound_samples)
        exact = normalize_centrality(read_centrality_values("scores/"*gn*"/exact.txt"))
        #if prog_sampler == "cm"
        apx_cc = read_centrality_values("scores/"*gn*"/"*algo*"_"*string(starting_sample)*".txt")
        #else
        #apx_cc = read_centrality_values("scores/"*gn*"/"*prog_sampler*"_"*algo*"_"*method*"_"*upper_bound_samples*"_"*string(starting_sample)*".txt")
        #end
        apx_bc::Array{Float64} = zeros(Float64,nv(g))
        mse_list = []
        sd_list = []
        for k in 1:trials
            tmp_scores = []
            for i in 1:nv(g)
                push!(tmp_scores,apx_cc[i+(nv(g)*(k-1))])
            end
            push!(mse_list,mean_sqaured_error(exact,tmp_scores))
            push!(sd_list,supremum_deviation(exact,tmp_scores))
        end
        for k in 1:trials
            for i in 1:nv(g)
                apx_bc[i] += apx_cc[i+(nv(g)*(k-1))]
            end
        end
        apx_bc = apx_bc .* [1/trials]
        for i in 1:nv(g)
            if apx_bc[i]> 1
                apx_bc[i] = 0.0
            end
        end
        SD = supremum_deviation(exact,apx_bc)
        MSE = mean_sqaured_error(exact,apx_bc)
        push!(results,[gn,algo,SD,MSE,mean(sd_list),std(sd_list),mean(mse_list),std(mse_list)])
    end
    mkpath("analysis/")
    header = false
    if !isfile("analysis/errors.txt")
        header = true
    end
    open("analysis/errors.txt","a") do file
        if header
            write(file,"Graph,Epsilon,Algorithm,SD,MSE,SDavg,SDstd,MSEavg,MSEstd\n")            
        end
        for res in results
            write(file,res[1]*","*string(target_epsilon)*","*res[2]*","*string(res[3])*","*string(res[4])*","*string(res[5])*","*string(res[6])*","*string(res[7])*","*string(res[8])*"\n")
        end
    end
end

function get_times(starting_sample::Int64,target_epsilon::Float64,datasets::Array{String},algo::String ="cmcera")
    results = []
    
    for graph in datasets
        gn = split(graph,".txt")[1]
        #println("Analyzing "*gn*" ε ",target_epsilon, " Algorithm ",algo," Progressive Sampler ",prog_sampler, " Path Optimality ",method, " Upper bound samples via ",upper_bound_samples)
        exact = read_time("times/"*gn*"/time_exact.txt")
        #if prog_sampler == "cm"
        apx_path = "times/"*gn*"/"*algo*"_"*string(starting_sample)*".txt"
        #else
        #apx_path = "times/"*gn*"/"*prog_sampler*"_"*algo*"_"*method*"_"*upper_bound_samples*"_"*string(starting_sample)*".txt"
        #end
        times = read_time(apx_path)
        samples =read_sample_size(apx_path)
        xi = read_xi(apx_path)
        push!(results,[gn,algo,exact,mean(times),std(times),mean(samples),std(samples),mean(xi),std(xi)])
    end
    mkpath("analysis/")
    header = false
    if !isfile("analysis/times.txt")
        header = true
    end
    open("analysis/times.txt","a") do file
        if header
            write(file,"Graph,Epsilon,Algorithm,ExactTime,ApxTimeMean,ApxTimeStd,SampleSizeMean,SampleSizeStd,XiMean,XiStd\n")            
        end
        for res in results
            write(file,res[1]*","*string(target_epsilon)*","*res[2]*","*string(res[3])*","*string(res[4])*","*string(res[5])*","*string(res[6])*","*string(res[7])*","*string(res[8])*","*string(res[9])*"\n")
        end
    end
end



function distance_metrics_stats(datasets,directed = true)
 
    all_results = []
    for graph in datasets
        gn = string(split(graph,".txt")[1])
        exact = read_distance_measures(gn)
        # write the procedure
        push!(all_results,append!([gn,directed],exact))
    end
    mkpath("analysis/")
    header = false
    if !isfile("analysis/distance_metrics.txt")
        header = true
    end
    open("analysis/distance_metrics.txt","a") do file

        if header
            write(file,"Graph,Directed,Diameter,Rho,ExactTime\n")            
        end
        for res in all_results
            write(file,res[1]*","*string(res[2])*","*string(res[3])*","*string(res[4])*","*string(res[5])*"\n")
        end
    end
end


function graph_properties(datasets,directed = true)
    results = []
    
    for graph in datasets
        g = load_graph("graphs/"*graph, directed,"\t")
        gn = split(graph,".txt")[1]
        #println("Analyzing "*gn*" ε ",target_epsilon, " Algorithm ",algo," Progressive Sampler ",prog_sampler, " Path Optimality ",method, " Upper bound samples via ",upper_bound_samples)
        @info("Stats of "*gn)
        exact = normalize_centrality(read_centrality_values("scores/"*gn*"/exact.txt"))
        push!(results,[gn,directed,nv(g),ne(g),sum(exact),maximum(exact)])
    end
    mkpath("analysis/")
    header = false
    if !isfile("analysis/graph_stats.txt")
        header = true
    end
    open("analysis/graph_stats.txt","a") do file

        if header
            write(file,"Graph,Directed,Nodes,Edges,SumPercCentr,MaxPercCentr\n")            
        end
        for res in results
            write(file,res[1]*","*string(res[2])*","*string(res[3])*","*string(res[4])*","*string(res[5])*","*string(res[6])*"\n")
        end
    end

end


function get_ranking_intersections(tk::Int64,starting_sample::Int64,target_epsilon::Float64,datasets::Array{String},algo::String ="cmcera",trials::Int64 = 10,directed::Bool = true)
    results = []
    for graph in datasets
        g = load_graph("graphs/"*graph,directed, "\t")
        gn = split(graph,".txt")[1]
        #println("Analyzing "*gn*" ε ",target_epsilon, " Algorithm ",algo," Progressive Sampler ",prog_sampler, " Path Optimality ",method, " Upper bound samples via ",upper_bound_samples)
        exact = normalize_centrality(read_centrality_values("scores/"*gn*"/exact.txt"))
        #if prog_sampler == "cm"
        apx_cc = read_centrality_values("scores/"*gn*"/"*algo*"_"*string(starting_sample)*".txt")
        #else
        #apx_cc = read_centrality_values("scores/"*gn*"/"*prog_sampler*"_"*algo*"_"*method*"_"*upper_bound_samples*"_"*string(starting_sample)*".txt")
        #end
        top_k_exact::Array{Tuple{Int64,Float64}} = Array{Tuple{Int64,Float64}}([])
        for u in 1:nv(g)
            push!(top_k_exact,(u,exact[u]))
        end
        sort!(top_k_exact, by=top_k_exact->(-top_k_exact[2],top_k_exact[1]))
        exact_rank = [u[1] for u in top_k_exact]
        apx_bc = zeros(nv(g))
        intersection_list = []
        min_h_list = []
        jaccard_list = []
        for k in 1:trials
            tmp_top_k_apx = []
            tmp_tbc = Array{Float64}([])
            for i in 1:nv(g)
                push!(tmp_top_k_apx,(i,apx_cc[i+(nv(g)*(k-1))]))
                push!(tmp_tbc,apx_cc[i+(nv(g)*(k-1))])
            end
            sort!(tmp_top_k_apx, by=tmp_top_k_apx->(-tmp_top_k_apx[2],tmp_top_k_apx[1]))
            tmp_top_k_rank = [u[1] for u in tmp_top_k_apx]
            push!(intersection_list,length(intersect(Set(exact_rank[1:tk]), Set(tmp_top_k_rank[1:tk]))))
            push!(min_h_list,min_h_k(exact,tmp_tbc,tk))
            push!(jaccard_list,jaccard(exact,tmp_tbc,tk)[end])
        end
        for k in 1:trials
            for i in 1:nv(g)
                apx_bc[i] += apx_cc[i+(nv(g)*(k-1))]
            end
        end
        apx_bc = apx_bc .* [1/trials]
        tmp_top_k_apx = []
        for i in 1:nv(g)
            push!(tmp_top_k_apx,(i,apx_bc[i]))
        end
        sort!(tmp_top_k_apx, by=tmp_top_k_apx->(-tmp_top_k_apx[2],tmp_top_k_apx[1]))
        tmp_top_k_rank = [u[1] for u in tmp_top_k_apx]
        push!(results,[gn,algo,tk,length(intersect(Set(exact_rank[1:tk]), Set(tmp_top_k_rank[1:tk]))),mean(intersection_list),std(intersection_list),min_h_k(exact,apx_bc,tk),mean(min_h_list),std(min_h_list),jaccard(exact,apx_bc,tk)[end],mean(jaccard_list),std(jaccard_list)])
    end
    mkpath("analysis/")
    header = false
    if !isfile("analysis/intersections.txt")
        header = true
    end
    open("analysis/intersections.txt","a") do file
        if header
            write(file,"Graph,Algorithm,Epsilon,k,Intersection,IntersectionMean,IntersectionStd,h,hMean,hStd,Jaccard,JaccardMean,JaccardStd\n")            
        end
        for res in results
            write(file,res[1]*","*res[2]*","*string(target_epsilon)*","*string(res[3])*","*string(res[4])*","*string(res[5])*","*string(res[6])*","*string(res[7])*","*string(res[8])*","*string(res[9])*","*string(res[10])*","*string(res[11])*","*string(res[12])*"\n")
        end
    end
end



function min_h_k(a::Array{Float64}, b::Array{Float64}, k::Int64)::Int64
    @assert length(a) == length(b) "The two rankings have different number of elements"
    if (k > length(a))
        k = length(a)
    end
    ai::Vector{Int64} = sortperm(a, rev=true)
    bi::Vector{Int64} = sortperm(b, rev=true)
    min_h_b_k_a::Int64 = 0
    for a_j in 1:k
        b_j::Int64 = 1
        while (bi[b_j] != ai[a_j])
            b_j = b_j + 1
        end
        if (b_j > min_h_b_k_a)
            min_h_b_k_a = b_j
        end
    end
    return min_h_b_k_a
end

function intersection(a::Array{Float64}, b::Array{Float64}, max_k::Int64)::Array{Float64}
    @assert length(a) == length(b) "The two rankings have different number of elements"
    if (max_k > length(a))
        max_k = length(a)
    end
    ai::Vector{Int64} = sortperm(a, rev=true)
    bi::Vector{Int64} = sortperm(b, rev=true)
    inter::Array{Float64} = zeros(max_k)
    for k in 1:max_k
        inter[k] = length(intersect(Set(ai[1:k]), Set(bi[1:k])))
    end
    return inter
end

function jaccard(a::Array{Float64}, b::Array{Float64}, max_k::Int64)::Array{Float64}
    @assert length(a) == length(b) "The two rankings have different number of elements"
    if (max_k > length(a))
        max_k = length(a)
    end
    ai::Vector{Int64} = sortperm(a, rev=true)
    bi::Vector{Int64} = sortperm(b, rev=true)
    jac::Array{Float64} = zeros(max_k)
    for k in 1:max_k
        jac[k] = length(intersect(Set(ai[1:k]), Set(bi[1:k]))) / length(union(Set(ai[1:k]), Set(bi[1:k])))
    end
    return jac
end