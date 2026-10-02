function branch_injections!(input, model, formulation::ACFormulation, ::Type{ModelBuild})
    V, angle = model[:V], model[:angle]
    n_branches = num_branches(input)

    terms = [branch_admittance_terms(input, k, branch_tap(input, model, formulation, k)) for k in 1:n_branches]
    from = [branch_from_bus_idx(input, k) for k in 1:n_branches]
    to = [branch_to_bus_idx(input, k) for k in 1:n_branches]

    @expression(model, P_inj_from[k in 1:n_branches], active_injection(V, angle, from[k], to[k], :from, terms[k]))
    @expression(model, P_inj_to[k in 1:n_branches], active_injection(V, angle, from[k], to[k], :to, terms[k]))
    @expression(model, Q_inj_from[k in 1:n_branches], reactive_injection(V, angle, from[k], to[k], :from, terms[k]))
    @expression(model, Q_inj_to[k in 1:n_branches], reactive_injection(V, angle, from[k], to[k], :to, terms[k]))

    return nothing
end

function branch_flow_limit!(input, model, formulation::ACFormulation, ::Type{ModelBuild})
    if !allows_generation_redispatch(formulation)
        return nothing
    end

    P_inj_from, P_inj_to = model[:P_inj_from], model[:P_inj_to]
    Q_inj_from, Q_inj_to = model[:Q_inj_from], model[:Q_inj_to]

    for k in 1:num_branches(input)
        max_flow = branch_max_flow(input, k)
        if !isfinite(max_flow)
            continue
        end

        @constraint(model, P_inj_from[k]^2 + Q_inj_from[k]^2 <= max_flow^2)
        @constraint(model, P_inj_to[k]^2 + Q_inj_to[k]^2 <= max_flow^2)
    end

    return nothing
end

function power_balance!(input, model, ::ACFormulation, ::Type{ModelBuild})
    n = num_buses(input)
    n_g = num_generators(input)
    n_l = num_loads(input)
    n_branches = num_branches(input)
    n_dc = num_dclines(input)

    V, Pg, Qg = model[:V], model[:Pg], model[:Qg]
    P_inj_from, P_inj_to = model[:P_inj_from], model[:P_inj_to]
    Q_inj_from, Q_inj_to = model[:Q_inj_from], model[:Q_inj_to]
    Pf, Pt, Qf, Qt = model[:Pf], model[:Pt], model[:Qf], model[:Qt]

    @constraint(
        model,
        active_power_balance[i in 1:n],
        sum(Pg[g] for g in 1:n_g if generator_bus_idx(input, g) == i; init = 0.0) -
        sum(load_active_power(input, l) for l in 1:n_l if load_bus_idx(input, l) == i; init = 0.0) ==
        sum(Pf[k] for k in 1:n_dc if dcline_from_bus_idx(input, k) == i; init = 0.0) +
        sum(Pt[k] for k in 1:n_dc if dcline_to_bus_idx(input, k) == i; init = 0.0) +
        sum(P_inj_from[k] for k in 1:n_branches if branch_from_bus_idx(input, k) == i; init = 0.0) +
        sum(P_inj_to[k] for k in 1:n_branches if branch_to_bus_idx(input, k) == i; init = 0.0)
    )

    @constraint(
        model,
        reactive_power_balance[i in 1:n],
        sum(Qg[g] for g in 1:n_g if generator_bus_idx(input, g) == i; init = 0.0) -
        sum(load_reactive_power(input, l) for l in 1:n_l if load_bus_idx(input, l) == i; init = 0.0) ==
        sum(Qf[k] for k in 1:n_dc if dcline_from_bus_idx(input, k) == i; init = 0.0) +
        sum(Qt[k] for k in 1:n_dc if dcline_to_bus_idx(input, k) == i; init = 0.0) +
        sum(Q_inj_from[k] for k in 1:n_branches if branch_from_bus_idx(input, k) == i; init = 0.0) +
        sum(Q_inj_to[k] for k in 1:n_branches if branch_to_bus_idx(input, k) == i; init = 0.0) -
        bus_shunt_susceptance(input, i) * V[i]^2
    )

    return nothing
end

function voltage_setpoint!(input, model, formulation::ACFormulation, ::Type{ModelBuild})
    if !has_voltage_setpoint(formulation)
        return nothing
    end

    V = model[:V]
    for g in 1:num_generators(input)
        bus_idx = generator_bus_idx(input, g)
        if bus_type(input, bus_idx) in (BusType.PV, BusType.Slack)
            JuMP.fix(V[bus_idx], generator_voltage_setpoint(input, g); force = true)
        end
    end

    return nothing
end
