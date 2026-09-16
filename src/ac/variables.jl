function voltage_variables!(input, model, ::ACFormulation, ::Type{ModelBuild})
    n = num_buses(input)

    @variable(
        model,
        V[i in 1:n],
        lower_bound = bus_min_voltage(input, i),
        upper_bound = bus_max_voltage(input, i),
        start = bus_voltage(input, i),
    )
    @variable(model, angle[i in 1:n], start = bus_voltage_angle(input, i))

    for i in 1:n
        if bus_type(input, i) == BusType.Slack
            JuMP.fix(angle[i], bus_voltage_angle(input, i); force = true)
        end
    end

    return nothing
end

function generator_variables!(input, model, formulation::ACFormulation, ::Type{ModelBuild})
    n_g = num_generators(input)
    redispatch = allows_generation_redispatch(formulation)

    @variable(
        model,
        Pg[g in 1:n_g],
        lower_bound = generator_min_active_power(input, g),
        upper_bound = generator_max_active_power(input, g),
        # start = generator_active_power(input, g),
    )
    @variable(
        model,
        Qg[g in 1:n_g],
        lower_bound = generator_min_reactive_power(input, g),
        upper_bound = generator_max_reactive_power(input, g),
        # start = generator_reactive_power(input, g),
    )

    for g in 1:n_g
        bus_idx = generator_bus_idx(input, g)

        # the slack bus's active power is never "fixed data" — it's whatever the network balance
        # equation requires, even in a non-redispatch (standard power-flow) formulation
        if !redispatch && bus_type(input, bus_idx) != BusType.Slack
            JuMP.fix(Pg[g], generator_active_power(input, g); force = true)
        end

        if !redispatch && !(bus_type(input, bus_idx) in (BusType.PV, BusType.Slack))
            JuMP.fix(Qg[g], generator_reactive_power(input, g); force = true)
        end
    end

    return nothing
end

function tap_variables!(input, model, formulation::ACFormulation, ::Type{ModelBuild})
    if !optimizes_tap(formulation)
        return nothing
    end

    n_branches = num_branches(input)
    @variable(
        model,
        tap_var[
            k in 1:n_branches;
            branch_is_transformer(input, k) && branch_min_tap_ratio(input, k) != branch_max_tap_ratio(input, k)
        ],
        lower_bound = branch_min_tap_ratio(input, k),
        upper_bound = branch_max_tap_ratio(input, k),
        # start = branch_tap_ratio(input, k),
    )

    @expression(
        model,
        tap[k in 1:n_branches],
        if branch_is_transformer(input, k) && branch_min_tap_ratio(input, k) != branch_max_tap_ratio(input, k)
            tap_var[k]
        else
            branch_tap_ratio(input, k)
        end,
    )

    return nothing
end

# the tap ratio to use for branch k: a JuMP variable when the formulation optimizes tap, otherwise
# the branch's fixed data
function branch_tap(input, model, formulation::ACFormulation, k::Int)
    if optimizes_tap(formulation)
        return model[:tap][k]
    else
        return branch_tap_ratio(input, k)
    end
end

# Pf is the flow leaving the from-bus terminal (before losses); Pt is what arrives at the to-bus
# terminal, after the linear loss model (constant + proportional to Pf) is applied. Qf/Qt are
# independent reactive power injections at each terminal, bounded by each terminal's own limits.
function dcline_variables!(input, model, formulation::ACFormulation, ::Type{ModelBuild})
    n_dc = num_dclines(input)
    redispatch = allows_generation_redispatch(formulation)

    @variable(
        model,
        Pf[k in 1:n_dc],
        lower_bound = -dcline_max_flow_to_from(input, k),
        upper_bound = dcline_max_flow_from_to(input, k),
        # start = dcline_active_power(input, k),
    )
    @variable(
        model,
        Qf[k in 1:n_dc],
        lower_bound = dcline_min_reactive_power_from(input, k),
        upper_bound = dcline_max_reactive_power_from(input, k),
        # start = dcline_reactive_power_from(input, k),
    )
    @variable(
        model,
        Qt[k in 1:n_dc],
        lower_bound = dcline_min_reactive_power_to(input, k),
        upper_bound = dcline_max_reactive_power_to(input, k),
        # start = dcline_reactive_power_to(input, k),
    )

    @expression(
        model,
        Pt[k in 1:n_dc],
        Pf[k] - (dcline_constant_loss(input, k) + dcline_loss_coefficient(input, k) * Pf[k])
    )

    if !redispatch
        for k in 1:n_dc
            JuMP.fix(Pf[k], dcline_active_power(input, k); force = true)
            JuMP.fix(Qf[k], dcline_reactive_power_from(input, k); force = true)
            JuMP.fix(Qt[k], dcline_reactive_power_to(input, k); force = true)
        end
    end

    return nothing
end
