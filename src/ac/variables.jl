# the buses whose voltage limits become soft under a formulation that minimizes the voltage
# violation: the terminals of every variable-tap transformer. Empty for every other formulation, so
# there the voltage limits stay hard bounds everywhere.
function voltage_violation_buses(input, formulation::ACFormulation)
    if !minimizes_voltage_violation(formulation)
        return Int[]
    end

    buses = Int[]
    for k in 1:num_branches(input)
        if branch_has_variable_tap(input, k)
            push!(buses, branch_from_bus_idx(input, k), branch_to_bus_idx(input, k))
        end
    end

    return sort!(unique!(buses))
end

function voltage_variables!(input, model, formulation::ACFormulation, ::Type{ModelBuild})
    n = num_buses(input)

    # V is created free and gets its limits below, because the buses whose violation is penalized
    # must be able to leave the limits in the first place
    @variable(model, V[i in 1:n], start = 1) # bus_voltage(input, i))
    @variable(model, angle[i in 1:n], start = 0) # bus_voltage_angle(input, i))

    violation_buses = voltage_violation_buses(input, formulation)

    for i in 1:n
        if !(i in violation_buses)
            JuMP.set_lower_bound(V[i], bus_min_voltage(input, i))
            JuMP.set_upper_bound(V[i], bus_max_voltage(input, i))
        end

        if bus_type(input, i) == BusType.Slack
            JuMP.fix(angle[i], bus_voltage_angle(input, i); force = true)
        end
    end

    # on the penalized buses the same limits come back as constraints with a slack measuring how far
    # past each one the voltage goes; the objective is what pushes those slacks to zero
    @variable(model, V_over[i in violation_buses] >= 0)
    @variable(model, V_under[i in violation_buses] >= 0)
    @constraint(
        model,
        V_over_constraint[i in violation_buses],
        V[i] <= bus_max_voltage(input, i) + V_over[i],
    )
    @constraint(
        model,
        V_under_constraint[i in violation_buses],
        V[i] >= bus_min_voltage(input, i) - V_under[i],
    )

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
        tap_var[k in 1:n_branches; branch_has_variable_tap(input, k)],
        lower_bound = branch_min_tap_ratio(input, k),
        upper_bound = branch_max_tap_ratio(input, k),
        # start = branch_tap_ratio(input, k),
    )

    @expression(
        model,
        tap[k in 1:n_branches],
        if branch_has_variable_tap(input, k)
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

# Every DC line variable is the power leaving the AC bus and entering the DC line at its terminal,
# the PowerModels convention that the DC line data also follows. Pf is the flow leaving the
# from-bus, and Pt is minus what arrives at the to-bus after the linear loss model (constant plus
# proportional to Pf) is applied. Qf and Qt are independent at each terminal, bounded by each
# terminal's own limits.
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
        -Pf[k] + (dcline_constant_loss(input, k) + dcline_loss_coefficient(input, k) * Pf[k])
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
