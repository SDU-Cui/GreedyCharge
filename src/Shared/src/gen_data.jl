"""
    gen_tasks(T::INT, K::INT)::StructArray{CTask}

Generate `K` tasks supposing the scheduling horizon is `T`.
"""
function gen_tasks(T::INT, K::INT)::StructArray{CTask}
    dist = Normal(FLOAT(7.0), FLOAT(2.0))
    P = rand(dist, K)
    clamp!(P, 4.0, 10.0)
    # ϕ = rand(Int(1):Int(3), K)
    n = div(K, 3) 
    ϕ = vcat(fill(1, n), fill(2, n), fill(3, n))
    if K % 3 != 0
        ϕ = vcat(ϕ, 1:(K % 3))
    end

    dist = Normal(FLOAT(0.95), 0.03)
    η = rand(dist, K)
    clamp!(η, 0.9, 0.99)
    l = rand(T÷3:T, K)
    clamp!(l, 1, T)
    h = zeros(INT, K)
    for k in 1:K
        h[k] = rand(ceil(INT, l[k]*0.2):ceil(INT, l[k]*0.8))
        h[k] = clamp(h[k], 1, l[k])
    end
    return StructArray{CTask}((ϕ, P, η, h, l))
end