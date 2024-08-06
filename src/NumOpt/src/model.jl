# set up the optimization model in JuMP

function build_constraints!(model::JuMP.Model, env::Environment)
    T = env.T
    K = get_num_EVs(env)
    @variable(model, x[1:T, 1:K], Bin)
    # each EV is allowed at most l steps for charging
    for k in 1:K
        l = env.tasks.l[k]
        fix.(x[l+1:T, k], 0, force = true)
    end
    # within the l steps, exactly h steps must be charging, i.e., of value 1
    for k in 1:K
        h = env.tasks.h[k]
        @constraint(model, sum(@view x[:, k]) == h)
    end
    # max power limit and phase imbalance constraint
    @expression(model, P_ch, x .* env.tasks.P')    # T×K
    @expression(model, P_phase, AffExpr.(env.base_load))   # T×3, now contain base load only
    # add the charging power via phase
    for p in 1:3
        @views add_to_expression!.(P_phase[:, p], sum(P_ch[:, env.tasks.ϕ .== p]; dims=2))
    end
    # max power limit
    @expression(model, P_total, sum(P_phase; dims=2))   # T×1
    @constraint(model, P_total .<= env.P_max)
    # phase imbalance
    @expression(model, P_phase_avg, sum(P_phase; dims=2) ./ 3)   # T×1
    @expression(model, P_phase_diff, zeros(AffExpr, T, 3))
    for (i, (a, b)) in enumerate([(1, 2), (2, 3), (3, 1)])
        P_phase_diff[:, i] .= @views P_phase[:, a] .- P_phase[:, b]
    end
    # we want |P_phase_diff| / P_phase_avg <= pi_max
    @constraint(model, .-P_phase_avg .* env.pi_max .<= P_phase_diff) 
    @constraint(model, P_phase_diff .<= P_phase_avg .* env.pi_max)

end


function set_objective!(model::JuMP.Model, env::Environment)
    P_ch = model[:P_ch]         # T×K
    @expression(model, cost, P_ch .* env.Δt .* env.price)   # T×K
    @expression(model, obj, sum(cost))
    @objective(model, Min, obj)
end

"""
    solve!(env::Environment; optimizer::Symbol=:Cbc, opt_options::AbstractDict=Dict())

Given a problem specified by `env`, optimize it with a solver and its options.
"""
function solve!(env::Environment; optimizer::Symbol=:Cbc, opt_options::AbstractDict=Dict())
    if optimizer === :Cbc
        model = JuMP.Model(optimizer_with_attributes(Cbc.Optimizer, opt_options...))
    elseif optimizer === :GLPK
        model = JuMP.Model(optimizer_with_attributes(GLPK.Optimizer, opt_options...))
    elseif optimizer ===:Gurobi
        model = JuMP.Model(optimizer_with_attributes(Gurobi.Optimizer, opt_options...))
    else
        error("Unrecognized optimizer: $optimizer")
    end
    build_constraints!(model, env)
    set_objective!(model, env)
    optimize!(model)
    solve_time = JuMP.solve_time(model)
    model[:found_solution] = true
    model[:opt_status] = "OPTIMAL"
    if termination_status(model) != JuMP.OPTIMAL 
        @warn "Global optimality failed: " termination_status(model)
        if termination_status(model) == JuMP.TIME_LIMIT
            # check whether the result is feasible
            x = value.(model[:x])
            rep = primal_feasibility_report(model, Dict(model[:x] .=> x); atol=1e-5)
            if isempty(rep)
                @info "Feasible solution is found depsite TIME_LIMIT"
                model[:opt_status] = "TIME_LIMIT_FEASIBLE"
            else
                @info "Feasibility check report:" rep
                model[:found_solution] = false
                model[:opt_status] = "TIME_LIMIT_AND_INFEASIBLE"
                @error "The candidate solution is infeasible in this TIME_LIMIT case."
            end
        else
            model[:found_solution] = false
            model[:opt_status] = termination_status(model) |> string
            @error "Bad status: $(termination_status(model))"
        end
    end
    return model, solve_time
end