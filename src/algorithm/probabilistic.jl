function upper_bound_average_diameter(delta::Float64,diam::Int64,tdd::Array{Int64},sample_size::Int64,verbose::Bool=false)::Float64
    avg_dist::Float64 = 0.0
    # Diam is the vertex diameter it's already diam +1 
    for i in 1:(diam+1)
        avg_dist += (i-1) * tdd[i]
    end
    avg_dist = avg_dist/sample_size
    # Upper bound using Bernstein bound
    log_term_avg_dist::Float64 = log(1. /delta)
    c_term_avg_dist::Float64 = (diam - 2)*log_term_avg_dist/sample_size
    average_diam_ub_b = avg_dist + c_term_avg_dist + sqrt(2*c_term_avg_dist*diam + c_term_avg_dist^2)
    var_estimate_diam::Float64 = 0.0
    # Upper bound using Empirical Bernstein bound
    for i in 1:diam
        var_estimate_diam+= (tdd[i] - avg_dist)^2
    end
    var_estimate_diam = var_estimate_diam/(sample_size-1)
    #=
    for i in 1:diam
        for j in (i+1):diam
            var_estimate_diam += ((i-1)-(j-1))^2 * tdd[i]/sample_size * tdd[j]/(sample_size-1) * norm
        end
    end
    =#
    log_term_avg_dist = log(2/delta)
    average_diam_ub_eb::Float64 = avg_dist + 7/3 * (diam -2) * log_term_avg_dist/sample_size + sqrt(2*var_estimate_diam*log_term_avg_dist / sample_size)
    avg_diam_upperbound = min(average_diam_ub_b,average_diam_ub_eb)
    if verbose
        @info("Average diameter "*string(avg_dist))
        @info("Average diameter UB (Bernstein) "*string(average_diam_ub_b))
        @info("Average diameter UB (Emp-Bernstein) "*string(average_diam_ub_eb))
        @info("Variance estimate average diameter "*string(var_estimate_diam))
        flush(stderr)
    end
    if diam -2 >0
        return min(avg_diam_upperbound,diam-2)
    else
        return avg_diam_upperbound
    end

end

function upper_bound_top_1_bc(top1_est_bc::Float64,delta::Float64,sample_size::Int64)::Float64
    log_term_top_1_bc::Float64 = log(1/delta)
    const_term_top_1_bc::Float64 =  log_term_top_1_bc/sample_size
    top1bc_upperbound::Float64  = top1_est_bc +  const_term_top_1_bc + sqrt(2*const_term_top_1_bc*top1_est_bc+const_term_top_1_bc^2)
    return min(1.0,top1bc_upperbound)
end

function number_samples_bound_function(x::Float64,rho::Float64,delta::Float64,eps::Float64)::Float64
    v_x::Float64 = x*(1-x)
    arg_h::Float64 = eps/v_x
    denom::Float64 = v_x * ((1+arg_h) * log(1+arg_h)-arg_h)
    return log(2*rho / (x*delta))/denom

end
function upper_bound_samples(max_bc::Float64,max_var::Float64, avg_dist::Float64,eps::Float64, delta_bound::Float64,debug::Bool = false)::Float64
    x_hat::Float64 = 0.0
    x_hat_l::Float64 = 0.5-sqrt(eps/3)
    x_hat_l = max(x_hat_l,0.0)
    x_hat_h::Float64 = 0.5 
    v_x::Float64 = 0.0
    arg_h::Float64 = 0.0
    f_val::Float64 = 0.0
    while x_hat_h - x_hat_l > 0.0001
        x_hat = (x_hat_h + x_hat_l) /2.0
        if debug
            println(" X_hat "*string(x_hat))
        end
        v_x = x_hat * (1-x_hat)
        arg_h = eps/v_x
        f_val = v_x * ((1+arg_h) * log(1+arg_h)-arg_h)
        if f_val <= 2*eps^2
            x_hat_h = x_hat
        else
            x_hat_l = x_hat
        end
    end
    x_hat = x_hat_h 
    x_h::Float64 = min(x_hat,max_bc)
    x_h_var::Float64 = 1.0
    if max_var < 0.25
        x_h_var = 0.5 - sqrt(0.25-max_var)
        x_h = min(x_h,x_h_var)
    end
    x_l::Float64 = x_h
    step::Float64 = x_h /1000.
    num_samples_bound_high::Float64 = number_samples_bound_function(x_h,avg_dist,delta_bound,eps) 
    num_samples_bound::Float64 =num_samples_bound_high+1
    while num_samples_bound > num_samples_bound_high
        x_l = x_h - step
        if x_l > 0
            num_samples_bound = number_samples_bound_function(x_l,avg_dist,delta_bound,eps)
            if num_samples_bound > num_samples_bound_high
                x_h = x_l
                num_samples_bound_high = num_samples_bound
            end
        else
            num_samples_bound = num_samples_bound_high -1
        end
    end
    return num_samples_bound_high
end


function _check_stopping_condition!(percolation::Array{Float64},wv::Array{Float64},last_stopping_samples::Int64,num_samples::Int64,eps::Float64,delta::Float64,iteration::Int64,second_phase::Bool,diam::Int64,tdd::Array{Int64},sample_size::Int64,mc_trials::Int64,partition_index::Array{Int64},partitions_ids_map::Dict{Int64,Int64},mcrade::Array{Float64},number_of_non_empty_partitions::Int64,omega::Vector{Float64},has_to_stop::Vector{Bool})
    n::Int64 = lastindex(percolation)
    num_samples_d::Float64 = num_samples
    delta_for_progressive_bound::Float64 = delta/2^iteration
    #println("Checking stopping condition at iteration "*string(iteration)*" sample size "*string(num_samples)*" δ = "*string(delta_for_progressive_bound))
    
    if second_phase
        avg_diam_upperbound::Float64 = upper_bound_average_diameter(delta_for_progressive_bound,diam,tdd,sample_size,false)
        top1_est_bc::Float64 = maximum(percolation)/num_samples_d
        top1bc_upperbound::Float64 = upper_bound_top_1_bc(top1_est_bc,delta_for_progressive_bound,sample_size)
        wimpy_var_upper_bound::Float64 = upper_bound_top_1_bc(maximum(wv)/num_samples_d,delta_for_progressive_bound,sample_size)
        max_num_samples::Int64 = trunc(Int,upper_bound_samples(top1bc_upperbound,wimpy_var_upper_bound,avg_diam_upperbound,eps,delta_for_progressive_bound))
        if last_stopping_samples > max_num_samples
            last_stopping_samples = max_num_samples
            
            omega[1] = last_stopping_samples
            @info("New stopping condition update, last stopping samples "*string(last_stopping_samples))
            flush(stderr)
        end
        if max_num_samples <= num_samples
            @info("New stopping condition TRUE")
            flush(stderr)
        end
    end
    
    
    sup_pcest_partition::Array{Float64} = zeros(n)
    sup_empwvar_partition::Array{Float64} = zeros(n)
    epsilon_partition::Array{Float64} = ones(n)
    max_mcera_partition::Array{Float64} = [-num_samples for i in 1:(mc_trials*n) ]
    # Update MCERA
    for i in 1:n
        v_rade_idx = i*mc_trials
        node_partition_idx = partition_index[i]
        mapped_partition_index = partitions_ids_map[node_partition_idx]
        sup_pcest_partition[mapped_partition_index] = max(sup_pcest_partition[mapped_partition_index],percolation[i] )
        sup_empwvar_partition[mapped_partition_index] = max(sup_empwvar_partition[mapped_partition_index],wv[i])
        mcera_partition_index = mc_trials*mapped_partition_index
        for j in 1:mc_trials
            max_mcera_partition[j+mcera_partition_index] = max(max_mcera_partition[j+mcera_partition_index] ,mcrade[v_rade_idx+j])
        end
    end
    mcera_partition_avg::Array{Float64} = zeros(number_of_non_empty_partitions)
    mcera_avg::Float64 = 0.0
    mcera_partition_index::Int64 = 0.0 
    delta_each_partition::Float64 = delta_for_progressive_bound/number_of_non_empty_partitions
    for i in 1:number_of_non_empty_partitions
        mcera_avg = 0.0
        mcera_partition_index = mc_trials *i
        for j in 1:mc_trials
            mcera_avg+=max_mcera_partition[mcera_partition_index+j]/mc_trials
        end
        mcera_avg = mcera_avg/num_samples_d
        mcera_partition_avg[i] = mcera_avg
        sup_emp_wimpy_var = sup_empwvar_partition[i]/num_samples_d
        current_eps = epsilon_mcrade(sup_emp_wimpy_var,mcera_avg,delta_each_partition,num_samples_d,mc_trials)
        epsilon_partition[i] = current_eps
    end
    sup_eps::Float64 = 0.0
    for i in 1:number_of_non_empty_partitions
        sup_eps = max(sup_eps,epsilon_partition[i])
    end
    if sup_eps <= eps
        @info("c-Monte Carlo R.A. STOPS with ξ : "*string(sup_eps))
        @info("c-Monte Carlo R.A. STOPS at iteration : "*string(iteration))
        @info("c-Monte Carlo R.A. STOPS at sample size : "*string(num_samples))
        flush(stderr)
    else
        @info("c-Monte Carlo R.A. ξ "*string(sup_eps)*" target ε  "*string(eps) )
        flush(stderr)
    end
    has_to_stop[1]= (sup_eps <= eps)
    return nothing
end

function epsilon_mcrade(sup_emp_wimpy_var,mcera,delta,num_samples,mc_trials)
    mcera = max(mcera,0.0)
    log_term_mcrade::Float64 = log(5/delta)/num_samples
    var_ub::Float64 = sup_emp_wimpy_var +log_term_mcrade+sqrt(log_term_mcrade^2 +2 * log_term_mcrade*sup_emp_wimpy_var) 
    era_ub::Float64 = mcera + sqrt(4*sup_emp_wimpy_var *log_term_mcrade / mc_trials)
    ra_ub::Float64 = era_ub + log_term_mcrade + sqrt(log_term_mcrade^2 + 2*log_term_mcrade * era_ub)
    eps_ub::Float64 = 2* ra_ub+ sqrt(2*log_term_mcrade*(var_ub+4*ra_ub))
    return eps_ub
end 

function get_next_stopping_sample(ss::Float64,iteration_index::Int64)
    ss = ss * 1.2
    iteration_index +=1
    return ss,iteration_index
end


function compute_xi(B,r)
    myf(x)= 1/x * log(sum([exp(x^2 * B[i]/(2*r^2)) for i in 1:lastindex(B)])) 
    model = Model(NLopt.Optimizer)
    set_optimizer_attribute(model, "algorithm", :LN_COBYLA)
    local_optimizer = NLopt.Opt(:LD_LBFGS, 1)
    set_optimizer_attribute(model, "local_optimizer", :LD_LBFGS)
    local_optimizer.lower_bounds = [floatmin(Float64)]
    local_optimizer.xtol_abs = 1e-4
    local_optimizer.ftol_abs =1e-6
    set_optimizer_attribute(model, "local_optimizer", local_optimizer)
    #set_optimizer_attribute(model, "algorithm", :LN_COBYLA)
    register(model,:myf, 1, myf, autodiff=true)
    @variable(model, x >= 0)
    @constraint(model,x >= 0)
    @NLobjective(model, Min,myf(x))
    JuMP.optimize!(model)
    return(objective_value(model))
end


function empirical_variance(tilde_b::Array{Float64}, sample_size::Int64, v::Int64)::Float64
    n::Int64 = div(length(tilde_b), sample_size)
    variance::Float64 = 0
    for i in 1:sample_size
        for j in (i+1):sample_size
            variance += (((tilde_b[(i-1)*n+v] - tilde_b[(j-1)*n+v]))^2)
        end
    end
    return variance / (sample_size * (sample_size - 1))
end

function theoretical_error_bound(tilde_b::Array{Float64},tb::Array{Float64}, sample_size::Int64, eta::Float64)::Float64
    n::Int64 = lastindex(tb)
    errors::Array{Float64} = zeros(n)
    Base.Threads.@threads for u in 1:n
        variance::Float64 = 0.0
        for i in 1:sample_size
            variance += (tilde_b[(i-1)*n+u] - (tb[u]/sample_size))^2
            #variance+= empirical_variance(tilde_b,sample_size,u)
        end
        variance = variance/(sample_size-1)
       # variance = variance
        errors[u] = sqrt(2 * variance * log(4 * n / eta) / sample_size) + 7 * log(4 * n / eta) / (3 * (sample_size - 1))
    end
    return maximum(errors)
end

function sample_couples(nodes::Array{Int64},X::Array{Float64},normalizer::Float64)::Tuple{Int64,Int64}
    s::Int64 = 0



end


@inline function weighted_sample_kappa_source(X::Array{Float64})::Int64
    n::Int64 = length(X)

    # Sort X and keep track of indices
    sorted_indices::Array{Int64} = sortperm(X)
    sorted_X::Array{Float64} = X[sorted_indices]

    # Compute cumulative sums of sorted X
    prefix_sum::Float64 = cumsum(sorted_X)
    s::Int64 = -1
    # Compute weights using efficient max(0, X[s] - X[i]) 
    weights::Array{Float64} = zeros(Float64,n)
    for rank in 1:n
        s = sorted_indices[rank]
        # Contribution of values less than X[s]
        weights[s] += (rank - 1) * sorted_X[rank] - (rank > 1 ? prefix_sum[rank - 1] : 0)
        # Contribution of values greater than X[s]
        weights[s] += (prefix_sum[end] - prefix_sum[rank]) - (n - rank) * sorted_X[rank]
    end

    # Normalize weights to probabilities
    total_weight::Float64 = sum(weights)
    probabilities::Array{Float64} = weights / total_weight

    # Use weighted sampling
    return sample(1:n, Weights(probabilities))
end


@inline function weighted_sample_kappa(X::Array{Float64})::Tuple{Int64,Int64}
    n::Int64 = length(X)

    # Sort X and keep track of indices
    sorted_indices::Array{Int64} = sortperm(X)
    sorted_X::Array{Float64} = X[sorted_indices]

    # Compute cumulative sums of sorted X
    prefix_sum::Array{Float64} = cumsum(sorted_X)

    # Compute weights using efficient max(0, X[s] - X[i]) logic
    weights::Array{Float64} = zeros(Float64,n)
    for rank in 1:n
        s = sorted_indices[rank]
        # Contribution of values less than X[s]
        weights[s] += (rank - 1) * sorted_X[rank] - (rank > 1 ? prefix_sum[rank - 1] : 0)
        # Contribution of values greater than X[s]
        weights[s] += (prefix_sum[end] - prefix_sum[rank]) - (n - rank) * sorted_X[rank]
    end

    # Normalize weights to probabilities
    total_weight::Float64 = sum(weights)
    probabilities::Array{Float64} = weights / total_weight

    # Use weighted sampling to select s
    s::Int64 = sample(1:n, Weights(probabilities))

    # Compute conditional probabilities for z given s
    conditional_weights::Array{Float64} = zeros(Float64,n)
    for j in 1:n
        if j != s
            conditional_weights[j] = max(0, X[s] - X[j])
        end
    end

    # Normalize conditional weights
    total_conditional_weight::Float64 = sum(conditional_weights)
    conditional_probabilities::Array{Float64} = Array{Float64}([])
    if total_conditional_weight > 0
        conditional_probabilities = conditional_weights / total_conditional_weight
    else
        conditional_probabilities = fill(1.0 / (n - 1), n)  # Handle edge case where all weights are zero
    end

    # Use weighted sampling to select z
    z::Int64 = sample(1:n, Weights(conditional_probabilities))

    return s, z
end
