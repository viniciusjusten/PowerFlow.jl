function voltage_variables!(input, model, ::ACFormulation, ::Type{ModelBuild})
    n = num_buses(input)

    @variable(
        model,
        bus_min_voltage(input, i) <= V[i in 1:n] <= bus_max_voltage(input, i),
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
        generator_min_active_power(input, g) <= Pg[g in 1:n_g] <= generator_max_active_power(input, g),
        start = generator_active_power(input, g),
    )
    @variable(
        model,
        generator_min_reactive_power(input, g) <= Qg[g in 1:n_g] <= generator_max_reactive_power(input, g),
        start = generator_reactive_power(input, g),
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
