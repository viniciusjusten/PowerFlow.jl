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

# reactive power has no meaning in this linearized (textbook DC) formulation, so only the DC
# line's active-power flow and linear loss model (constant + proportional to Pf) are represented
function dcline_variables!(input, model, ::Type{ModelBuild})
    n_dc = num_dclines(input)

    @variable(
        model,
        Pf[k in 1:n_dc],
        lower_bound = -dcline_max_flow_to_from(input, k),
        upper_bound = dcline_max_flow_from_to(input, k),
        start = dcline_active_power(input, k),
    )

    @expression(
        model,
        Pt[k in 1:n_dc],
        Pf[k] - (dcline_constant_loss(input, k) + dcline_loss_coefficient(input, k) * Pf[k])
    )

    return nothing
end
