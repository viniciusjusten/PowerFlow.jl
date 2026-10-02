# Every branch and DC line flow below is the power leaving the AC bus at that terminal, the same
# convention PowerModels uses. A branch carrying power from its from-bus to its to-bus therefore has
# a positive `_from` value and a negative `_to` value, and the sum of the two is the branch loss. A DC
# line works the same way, so its `_to` value is minus the power the line delivers to the to-bus.
struct ACSolution
    voltage::Vector{Float64}
    voltage_angle::Vector{Float64}
    active_power_generation::Vector{Float64}
    reactive_power_generation::Vector{Float64}
    branch_active_power_from::Vector{Float64}
    branch_active_power_to::Vector{Float64}
    branch_reactive_power_from::Vector{Float64}
    branch_reactive_power_to::Vector{Float64}
    dcline_active_power_from::Vector{Float64}
    dcline_active_power_to::Vector{Float64}
    dcline_reactive_power_from::Vector{Float64}
    dcline_reactive_power_to::Vector{Float64}
    objective_value::Float64
    termination_status::Any
    # TODO - add tap solution and duals
end

struct DCSolution
    voltage_angle::Vector{Float64}
    active_power_generation::Vector{Float64}
    branch_active_power_from::Vector{Float64}
    branch_active_power_to::Vector{Float64}
    dcline_active_power_from::Vector{Float64}
    dcline_active_power_to::Vector{Float64}
    objective_value::Float64
    termination_status::Any
    # TODO - add duals
end

function extract_ac_solution(model::JuMP.Model)
    return ACSolution(
        JuMP.value.(model[:V]),
        JuMP.value.(model[:angle]),
        JuMP.value.(model[:Pg]),
        JuMP.value.(model[:Qg]),
        JuMP.value.(model[:P_inj_from]),
        JuMP.value.(model[:P_inj_to]),
        JuMP.value.(model[:Q_inj_from]),
        JuMP.value.(model[:Q_inj_to]),
        JuMP.value.(model[:Pf]),
        JuMP.value.(model[:Pt]),
        JuMP.value.(model[:Qf]),
        JuMP.value.(model[:Qt]),
        JuMP.objective_value(model),
        JuMP.termination_status(model),
    )
end

# the DC branch flows are lossless, so the to-bus side is the from-bus flow mirrored
function extract_dc_solution(model::JuMP.Model)
    flow = JuMP.value.(model[:flow])
    return DCSolution(
        JuMP.value.(model[:angle]),
        JuMP.value.(model[:Pg]),
        flow,
        -flow,
        JuMP.value.(model[:Pf]),
        JuMP.value.(model[:Pt]),
        JuMP.objective_value(model),
        JuMP.termination_status(model),
    )
end

"""
    ac_branch_flows(input, V, angle)

Evaluate the active and reactive power leaving each terminal of every branch for a given voltage
state, using the fixed tap ratio of each branch. Returns the vectors `(P_from, P_to, Q_from, Q_to)`.
"""
function ac_branch_flows(input, V, angle)
    n_branches = num_branches(input)
    P_from = zeros(n_branches)
    P_to = zeros(n_branches)
    Q_from = zeros(n_branches)
    Q_to = zeros(n_branches)

    for k in 1:n_branches
        P_from[k] = active_injection(input, k, V, angle, :from)
        P_to[k] = active_injection(input, k, V, angle, :to)
        Q_from[k] = reactive_injection(input, k, V, angle, :from)
        Q_to[k] = reactive_injection(input, k, V, angle, :to)
    end

    return P_from, P_to, Q_from, Q_to
end

"""
    dcline_scheduled_flows(input)

Return the active and reactive power leaving each terminal of every DC line according to its
scheduled data, as the vectors `(P_from, P_to, Q_from, Q_to)`. This is what the power flows that do
not dispatch the DC lines use.
"""
function dcline_scheduled_flows(input)
    n_dc = num_dclines(input)
    P_from = [dcline_active_power(input, k) for k in 1:n_dc]
    P_to = [-dcline_to_bus_active_power(input, k) for k in 1:n_dc]
    Q_from = [dcline_reactive_power_from(input, k) for k in 1:n_dc]
    Q_to = [dcline_reactive_power_to(input, k) for k in 1:n_dc]
    return P_from, P_to, Q_from, Q_to
end

"""
    dc_branch_flows(input, angle)

Evaluate the active power flow of every branch under the textbook DC approximation, the same one
used by `dc_susceptance_matrix` and by the DC optimal power flow model.
"""
function dc_branch_flows(input, angle)
    return [
        (angle[branch_from_bus_idx(input, k)] - angle[branch_to_bus_idx(input, k)]) / branch_reactance(input, k)
        for k in 1:num_branches(input)
    ]
end
