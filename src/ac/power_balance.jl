function branch_injections!(input, model, formulation::ACFormulation, ::Type{ModelBuild})
    V, angle = model[:V], model[:angle]
    n_branches = PSD.num_branches(input)

    terms = [branch_admittance_terms(input, k, branch_tap(input, model, formulation, k)) for k in 1:n_branches]
    from = [PSD.branch_from_bus_idx(input, k) for k in 1:n_branches]
    to = [PSD.branch_to_bus_idx(input, k) for k in 1:n_branches]

    # the flows are variables tied to the voltages by equality constraints, as in PowerModels, so the
    # power balance stays linear in them and the thermal limits are quadratics in the flows instead
    # of in the trigonometric expressions
    @variable(model, P_inj_from[k in 1:n_branches])
    @variable(model, P_inj_to[k in 1:n_branches])
    @variable(model, Q_inj_from[k in 1:n_branches])
    @variable(model, Q_inj_to[k in 1:n_branches])

    @constraint(
        model,
        P_inj_from_definition[k in 1:n_branches],
        P_inj_from[k] == active_injection(V, angle, from[k], to[k], :from, terms[k]),
    )
    @constraint(
        model,
        P_inj_to_definition[k in 1:n_branches],
        P_inj_to[k] == active_injection(V, angle, from[k], to[k], :to, terms[k]),
    )
    @constraint(
        model,
        Q_inj_from_definition[k in 1:n_branches],
        Q_inj_from[k] == reactive_injection(V, angle, from[k], to[k], :from, terms[k]),
    )
    @constraint(
        model,
        Q_inj_to_definition[k in 1:n_branches],
        Q_inj_to[k] == reactive_injection(V, angle, from[k], to[k], :to, terms[k]),
    )

    return nothing
end

function branch_flow_limit!(input, model, formulation::ACFormulation, ::Type{ModelBuild})
    if !allows_generation_redispatch(formulation)
        return nothing
    end

    P_inj_from, P_inj_to = model[:P_inj_from], model[:P_inj_to]
    Q_inj_from, Q_inj_to = model[:Q_inj_from], model[:Q_inj_to]

    for k in 1:PSD.num_branches(input)
        max_flow = PSD.branch_max_flow(input, k)
        if !isfinite(max_flow)
            continue
        end

        @constraint(model, P_inj_from[k]^2 + Q_inj_from[k]^2 <= max_flow^2)
        @constraint(model, P_inj_to[k]^2 + Q_inj_to[k]^2 <= max_flow^2)

        # redundant with the limit above at the solution, but Ipopt keeps every iterate inside
        # variable bounds, while it only enforces the constraint at convergence
        for flow in (P_inj_from[k], P_inj_to[k], Q_inj_from[k], Q_inj_to[k])
            JuMP.set_lower_bound(flow, -max_flow)
            JuMP.set_upper_bound(flow, max_flow)
        end
    end

    return nothing
end

# the angle of the from-bus minus the angle of the to-bus of every branch stays within the branch's
# limits. Like the thermal limits, this only applies to the formulations that redispatch, so that a
# plain power flow solves whatever the scheduled data imposes.
function branch_angle_difference_limit!(input, model, formulation::ACFormulation, ::Type{ModelBuild})
    if !allows_generation_redispatch(formulation)
        return nothing
    end

    angle = model[:angle]
    n_branches = PSD.num_branches(input)
    from = [PSD.branch_from_bus_idx(input, k) for k in 1:n_branches]
    to = [PSD.branch_to_bus_idx(input, k) for k in 1:n_branches]

    @constraint(
        model,
        angle_difference_min[k in 1:n_branches],
        angle[from[k]] - angle[to[k]] >= PSD.branch_min_angle_difference(input, k),
    )
    @constraint(
        model,
        angle_difference_max[k in 1:n_branches],
        angle[from[k]] - angle[to[k]] <= PSD.branch_max_angle_difference(input, k),
    )

    return nothing
end

function power_balance!(input, model, formulation::ACFormulation, ::Type{ModelBuild})
    n = PSD.num_buses(input)
    n_g = PSD.num_generators(input)
    n_l = PSD.num_loads(input)
    n_branches = PSD.num_branches(input)
    n_dc = PSD.num_dclines(input)

    V, Pg, Qg = model[:V], model[:Pg], model[:Qg]
    P_inj_from, P_inj_to = model[:P_inj_from], model[:P_inj_to]
    Q_inj_from, Q_inj_to = model[:Q_inj_from], model[:Q_inj_to]
    Pf, Pt, Qf, Qt = model[:Pf], model[:Pt], model[:Qf], model[:Qt]

    # the fictitious reactive injection of each bus, which is zero everywhere except at the buses
    # that get reactive balance slacks
    Q_violation = zeros(JuMP.AffExpr, n)
    if minimizes_reactive_balance_violation(formulation)
        Q_deficit, Q_surplus = model[:Q_deficit], model[:Q_surplus]
        for i in variable_tap_buses(input)
            Q_violation[i] = Q_deficit[i] - Q_surplus[i]
        end
    end

    @constraint(
        model,
        active_power_balance[i in 1:n],
        sum(Pg[g] for g in 1:n_g if PSD.generator_bus_idx(input, g) == i; init = 0.0) -
        sum(PSD.load_active_power(input, l) for l in 1:n_l if PSD.load_bus_idx(input, l) == i; init = 0.0) ==
        sum(Pf[k] for k in 1:n_dc if PSD.dcline_from_bus_idx(input, k) == i; init = 0.0) +
        sum(Pt[k] for k in 1:n_dc if PSD.dcline_to_bus_idx(input, k) == i; init = 0.0) +
        sum(P_inj_from[k] for k in 1:n_branches if PSD.branch_from_bus_idx(input, k) == i; init = 0.0) +
        sum(P_inj_to[k] for k in 1:n_branches if PSD.branch_to_bus_idx(input, k) == i; init = 0.0)
    )

    @constraint(
        model,
        reactive_power_balance[i in 1:n],
        sum(Qg[g] for g in 1:n_g if PSD.generator_bus_idx(input, g) == i; init = 0.0) -
        sum(PSD.load_reactive_power(input, l) for l in 1:n_l if PSD.load_bus_idx(input, l) == i; init = 0.0) +
        Q_violation[i] ==
        sum(Qf[k] for k in 1:n_dc if PSD.dcline_from_bus_idx(input, k) == i; init = 0.0) +
        sum(Qt[k] for k in 1:n_dc if PSD.dcline_to_bus_idx(input, k) == i; init = 0.0) +
        sum(Q_inj_from[k] for k in 1:n_branches if PSD.branch_from_bus_idx(input, k) == i; init = 0.0) +
        sum(Q_inj_to[k] for k in 1:n_branches if PSD.branch_to_bus_idx(input, k) == i; init = 0.0) -
        PSD.bus_shunt_susceptance(input, i) * V[i]^2
    )

    return nothing
end

function voltage_setpoint!(input, model, formulation::ACFormulation, ::Type{ModelBuild})
    if !has_voltage_setpoint(formulation)
        return nothing
    end

    V = model[:V]
    for g in 1:PSD.num_generators(input)
        bus_idx = PSD.generator_bus_idx(input, g)
        if PSD.bus_type(input, bus_idx) in (PSD.BusType.PV, PSD.BusType.Slack)
            JuMP.fix(V[bus_idx], PSD.generator_voltage_setpoint(input, g); force = true)
        end
    end

    return nothing
end
