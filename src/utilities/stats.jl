
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
        gn = split(graph,".txt")[1]
        exact = read_distance_measures(gn,path_optimality)
        # write the procedure
        push!(all_results,append!([gn,string(directed)],exact))
    end
    return all_results
end