module NumOpt

using JuMP, Cbc, GLPK, Gurobi
using Revise, Shared

include("model.jl")

get_value(T::Type, x) = JuMP.value.(x) .|> T 
get_value(x) = JuMP.value.(x)

get_decision_matrix(model::JuMP.Model) = INT.(value.(model[:x]) .> 0.5)

export solve!, get_value, get_decision_matrix

end # module NumOpt
