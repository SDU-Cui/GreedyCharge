using JLD2

Base.@kwdef mutable struct GreedyConfig
    p_max::Float64 = 6000
    phi_max::Float64 = 4e-2
end

struct Pi
    indices::Vector{Int}
    variables::Matrix{Int}
    power::Vector{Float64}
end


function reshape_X(X::AbstractMatrix{INT})
    T, K = size(X)
    pi_A = env.phases[1]
    var_A = Matrix{INT}(undef, 0, 0)
    pi_B = env.phases[2]
    var_B = Matrix{INT}(undef, 0, 0)
    pi_C = env.phases[3]
    var_C = Matrix{INT}(undef, 0, 0)
    for k in env.phases[1]
        var_A = vcat(var_A, X[:, k])
    end

    for k in env.phases[2]
        var_B = vcat(var_B, X[:, k])
    end  

    for k in env.phases[3]
        var_C = vcat(var_C, X[:, k])
    end

    power = compute_power(X, env)       # power::Matrix (T,3)
    power_A = power[:, 1] |> vec
    power_B = power[:, 2] |> vec
    power_C = power[:, 3] |> vec

    return Dict(:A => Pi(pi_A, var_A, power_A), :B => Pi(pi_B, var_B, power_B), :C => Pi(pi_C, var_C, power_C))
end


function repair!(X::AbstractMatrix{INT}, env::Environment, ti::INT)
    # pi = reshape_X(X)                          # distinguish all EVs into three phase
    p = compute_power(X, env)                    # compute all time steps A,B,C gross load 
    # pt = compute_power(X[ti, :], env)          # compute time step t A,B,C gross phase
    xi = X[ti, :]
    pt = p[ti, :]                                # compute time step t A,B,C gross phase
    p = sum(pt)
    max_pi = which_max_phase(pt)

    while p > env.P_max
        rndstop!(xi, pt, env, max_pi)
        p = sum(pt)
        max_pi = which_max_phase(pt)
    end
    
    imb = compute_imbalance(pt)

    while imb > env.pi_max
        rndstop!(xi, pt, env, max_pi)
        imb = compute_imbalance(pt)
        max_pi = which_max_phase(pt)
    end

    X[ti, :] = xi
end


function repair!(h::AVI, xt::AbstractVector{INT}, env::Environment, ti::INT)
    pt = compute_power(xt, env, ti)                    # compute time step t, A,B,C gross load 
    p = sum(pt)
    max_pi = which_max_phase(pt)

    while p > env.P_max
        rndstop!(h, xt, pt, env, max_pi)
        p = sum(pt)
        max_pi = which_max_phase(pt)
    end
    
    imb = compute_imbalance(pt)

    while imb > env.pi_max
        rndstop!(h, xt, pt, env, max_pi)
        imb = compute_imbalance(pt)
        max_pi = which_max_phase(pt)
    end
end


function create_individual_greedy!(h::AVI, env::Environment, ti::Int)
    # h = env.tasks.h             # remaining number of steps still required to fullfil charging
    l = env.tasks.l             # remaining duration of stay in charging station allowed for charging
    EVs = get_num_EVs(env)
    xt = zeros(INT, EVs)
    
    for v in 1:EVs
        if h[v] == 0
            continue
        else
            if ti <= l[v]
                xt[v] = 1
                h[v] -= 1
            end
        end
    end

    return xt
end


function compute_power(arr::AbstractMatrix{INT}, env::Environment)
    # compute gross load at all time steps; power::(T, 3)
    base_load = env.base_load
    P = env.tasks.P
    pi_A = env.phases[1]
    pi_B = env.phases[2]
    pi_C = env.phases[3]
    power_list = arr .* P'
    power_A = sum(power_list[:, pi_A[:]], dims=2) |> vec
    power_B = sum(power_list[:, pi_B[:]], dims=2) |> vec
    power_C = sum(power_list[:, pi_C[:]], dims=2) |> vec
    EVspower = hcat(power_A, power_B, power_C)
    power = base_load + EVspower

    return power
end


function compute_power(arr::AbstractVector{INT}, env::Environment, ti::Int)
    # compute gross load at time step t; power::(1, 3)
    base_load = env.base_load[ti, :]             # step t base load::(1, 3)
    P = env.tasks.P
    pi_A = env.phases[1]
    pi_B = env.phases[2]
    pi_C = env.phases[3]
    power_list = P .* arr
    power_A = sum(power_list[pi_A])
    power_B = sum(power_list[pi_B])
    power_C = sum(power_list[pi_C])
    EVspower = Vector{FLOAT}([power_A, power_B, power_C])
    power = base_load + EVspower

    return power
end


function compute_imbalance(power::AVF)
    minp, maxp = extrema(power)
    avgp = sum(power) / 3

    return (maxp - minp) / avgp
end


function compute_imbalance(power::AMF)
    maxp = mapslices(maximum, power, dims=2)
    minp = mapslices(minimum, power, dims=2)
    avgp = sum(power, dims=2) / 3

    return ((maxp - minp) ./ avgp) |> vec
end


function which_max_phase(power::AVF)
    max_phase = argmax(power)

    return max_phase
end


function rndstop!(h::AVI, xi::Vector{INT}, pt::Vector{FLOAT}, env::Environment, max_pi::Int)
    max_phases = env.phases[max_pi]
    powers = env.tasks.P
    # ones_indices = findall(x -> x == 1, xi[max_phases])
    ones_indices = filter(i -> xi[i] == 1, max_phases)
    rnd_v = rand(ones_indices)

    xi[rnd_v] = 0
    h[rnd_v] += 1
    pt[max_pi] -= powers[rnd_v]
end


function minhstop!(h::AVI, xi::AVI, pt::AVF, env::Environment, max_pi::Int)
    max_phases = env.phases[max_pi]
    powers = env.tasks.P
    ones_indices = filter(i -> xi[i] == 1, max_phases)
    minh_index = sortperm(h[ones_indices])[1]
    minh_v = ones_indices[minh_index]

    xi[minh_v] = 0
    h[minh_v] += 1
    pt[max_pi] -= powers[minh_v]
end


function rndstart!(h::AVI, X::AMI, env::Environment)
    
end


function power_violation(X::AMI, env::Environment)
    power = compute_power(X, env)
    gross_power = sum(power, dims=2) |> vec

    return gross_power / env.P_max
end


function imbalance_violation(X::AMI, env::Environment)
    power = compute_power(X, env)
    imb = compute_imbalance(power)

    return imb / env.pi_max
end


function repairSOC!(h::AVI, X::AMI, env::Environment)
    err_indicies = findall(x -> x > 0, h)
    for i in err_indicies
        power_vio = power_violation(X, env)
        imb_vio = imbalance_violation(X, env)
        vio_sort = sortperm(imb_vio)
        for t in vio_sort
            if h[i] == 0
                break
            elseif t > env.tasks.l[i]
                continue
            elseif X[t, i] == 0
                if power_vio[t] + env.tasks.P[i] / env.P_max > 1
                    continue
                else 
                    h[i] -= 1
                    X[t, i] = 1
                end
            else
                continue
            end
        end
    end
end