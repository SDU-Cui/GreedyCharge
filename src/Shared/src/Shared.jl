module Shared

using StructArrays
using Distributions


const INT = Int64
const FLOAT = Float64
const AVF = AbstractVector{FLOAT}
const AMF = AbstractMatrix{FLOAT}
const AVI = AbstractVector{INT}
const AMI = AbstractMatrix{INT}

include("struct.jl")
include("gen_data.jl")

"""
    compute_phase_imbalance(Pt::AVF)

Given three-phase power `Pt`, compute the phase imbalance.
"""
function compute_phase_imbalance(Pt::AVF)
    min_v, max_v = extrema(Pt)
    avg_v = sum(Pt) / 3
    pim = (max_v - min_v) / avg_v
    return pim
end

export compute_phase_imbalance, INT, FLOAT, AVF, AMF, AVI, AMI, CTask, CTasks, Environment,
    get_num_EVs, get_power, get_phase, to_struct_array, gen_tasks, validate_solution_feasibility

end # module Shared
