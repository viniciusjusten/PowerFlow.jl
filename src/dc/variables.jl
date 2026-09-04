function angle_variables!(input, model, ::Type{ModelBuild})
    n = num_buses(input)

    @variable(model, angle[i in 1:n], start = bus_voltage_angle(input, i))

    for i in 1:n
        if bus_type(input, i) == BusType.Slack
            JuMP.fix(angle[i], bus_voltage_angle(input, i); force = true)
        end
    end

    return nothing
end

function generator_variables!(input, model, ::Type{ModelBuild})
    n_g = num_generators(input)

    @variable(
        model,
        generator_min_active_power(input, g) <= Pg[g in 1:n_g] <= generator_max_active_power(input, g),
        start = generator_active_power(input, g),
    )

    return nothing
end
