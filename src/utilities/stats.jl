
#using Conda
using PyCall
#Conda.add("scipy")

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









