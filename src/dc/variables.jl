function angle_variables!(input, model, ::Type{ModelBuild})
    n = PSD.num_buses(input)

    @variable(model, angle[i in 1:n], start = PSD.bus_voltage_angle(input, i))

    for i in 1:n
        if PSD.bus_type(input, i) == PSD.BusType.Slack
            JuMP.fix(angle[i], PSD.bus_voltage_angle(input, i); force = true)
        end
    end

    return nothing
end

function generator_variables!(input, model, ::Type{ModelBuild})
    n_g = PSD.num_generators(input)

    @variable(
        model,
        PSD.generator_min_active_power(input, g) <= Pg[g in 1:n_g] <= PSD.generator_max_active_power(input, g),
        # start = generator_active_power(input, g),
    )

    return nothing
end

# reactive power has no meaning in this linearized (textbook DC) formulation, so only the DC
# line's active-power flow and linear loss model (constant + proportional to Pf) are represented.
# As in the AC model, Pf and Pt are the power leaving the AC bus at each terminal, so Pt is minus
# what arrives at the to-bus.
function dcline_variables!(input, model, ::Type{ModelBuild})
    n_dc = PSD.num_dclines(input)

    @variable(
        model,
        Pf[k in 1:n_dc],
        lower_bound = -PSD.dcline_max_flow_to_from(input, k),
        upper_bound = PSD.dcline_max_flow_from_to(input, k),
        # start = dcline_active_power(input, k),
    )

    @expression(
        model,
        Pt[k in 1:n_dc],
        -Pf[k] + (PSD.dcline_constant_loss(input, k) + PSD.dcline_loss_coefficient(input, k) * Pf[k])
    )

    return nothing
end
