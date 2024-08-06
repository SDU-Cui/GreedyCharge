"""
    Base.@kwdef mutable struct CTask

Represent a charging task.
"""
Base.@kwdef mutable struct CTask
    ϕ::INT          # phase (1, 2, 3) for ('A', 'B', 'C')
    P::FLOAT        # charging power, kW
    η::FLOAT        # charging efficiency
    h::INT          # remaining number of steps still required to fullfil charging
    l::INT          # remaining duration of stay in charging station allowed for charging
end

const CTasks = StructArray{CTask}

to_struct_array(task_list::AbstractVector{CTask})::CTasks = StructArray(task_list)


struct Environment
    T::INT                      # horizon length 
    Δt::FLOAT                   # unit: hour
    P_max::FLOAT                # total power limit
    pi_max::FLOAT               # phase imbalance limit
    # exogenous state variables, each vector of length `T`
    base_load::Matrix{FLOAT}    # T × 3
    price::Vector{FLOAT}

    tasks::CTasks
    phases::Vector{Vector{Int}} # phases[1] stores indices of EVs in phase A 
    # powers::Vector{FLOAT}       # charging power of all EVs
end 


function Environment(T::Integer, Δt::AbstractFloat, P_max::AbstractFloat,
        pi_max::AbstractFloat, base_load::AbstractMatrix{<:AbstractFloat}, 
        price::AbstractVector{<:AbstractFloat}, tasks::CTasks)
    evs = collect(INT, 1:length(tasks))
    phases = Vector{INT}[]
    for p in 1:3
        push!(phases, evs[tasks.ϕ .== p])
    end
    return Environment(T, Δt, P_max, pi_max, base_load, price, tasks, phases)
end


get_num_EVs(env::Environment) = length(env.tasks)
get_phase(env::Environment, k::Integer) = env.tasks.ϕ[k]
get_power(env::Environment, k::Integer) = env.tasks.P[k]


function validate_solution_feasibility(x::AMI, env::Environment)
    K = get_num_EVs(env)
    # time constraint
    for k in 1:K
        l = env.tasks.l[k]
        @assert all(e == 0 for e in @view x[l+1:end, k])
    end
    # SOC constraint
    h = sum(x; dims=1) |> vec
    @assert all(h .== env.tasks.h)
    # power limit
    P_ch = x .* env.tasks.P'
    P_phase = copy(env.base_load)
    for p in 1:3
        @views P_phase[:, p] .+= sum(P_ch[:, env.tasks.ϕ .== p]; dims=2)
    end
    P_total = sum(P_phase; dims=2)
    @assert all(P_total .<= env.P_max + 1e-2)
    # phase imbalance constraint
    pim = compute_phase_imbalance.(eachrow(P_phase))
    @assert all(pim .<= env.pi_max + 1e-3)
    println("👍The solution is feasible!")
end
