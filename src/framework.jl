function get_value_greedy(X::AMI, env::Environment)
    P = copy(env.tasks.P)
    price = copy(env.price)
    value = sum(X * P * env.Δt .* price)
    
    return value
end

function run_greedy(env::Environment; jld2_file::String="")
    price = env.price           # ::Vector{FLOAT}
    sorted_p = sortperm(price)  # order of price from high to low
    X = zeros(INT, env.T, size(env.tasks.P)[1])
    h = copy(env.tasks.h)       # expected SOC
    
    start_time = time()
    for ti in sorted_p
        xt = create_individual_greedy!(h, env, ti)
        repair!(h, xt, env, ti)
        X[ti, :] = xt
    end
    end_time = time()
    solve_time = end_time - start_time
    
    return X, solve_time
end