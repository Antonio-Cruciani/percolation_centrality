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
        println("Average diameter "*string(avg_dist))
        println("Average diameter UB (Bernstein) "*string(average_diam_ub_b))
        println("Average diameter UB (Emp-Bernstein) "*string(average_diam_ub_eb))
        println("Variance estimate average diameter "*string(var_estimate_diam))
        flush(stdout)
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
