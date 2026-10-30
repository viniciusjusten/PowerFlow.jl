function objective_function!(input, model, ::Union{ACOPFFormulation, ACOPFTapFormulation}, ::Type{ModelBuild})
    Pg = model[:Pg]
    @objective(model, Min, sum(PSD.generator_cost(input, g) * Pg[g] for g in 1:PSD.num_generators(input)))
    return nothing
end

function objective_function!(input, model, ::ACOPFMinRedispatchFormulation, ::Type{ModelBuild})
    Pg = model[:Pg]
    @objective(
        model,
        Min,
        sum((Pg[g] - PSD.generator_active_power(input, g))^2 for g in 1:PSD.num_generators(input))
    )
    return nothing
end

function objective_function!(input, model, ::ACPFFeasibilityFormulation, ::Type{ModelBuild})
    @objective(model, Min, 0.0)
    return nothing
end

function objective_function!(input, model, ::ACPFPmisVmisFormulation, ::Type{ModelBuild})
    V, Pg = model[:V], model[:Pg]
    @objective(
        model,
        Min,
        sum(
            (Pg[g] - PSD.generator_active_power(input, g))^2 +
            (V[PSD.generator_bus_idx(input, g)] - PSD.generator_voltage_setpoint(input, g))^2 for
            g in 1:PSD.num_generators(input)
        )
    )
    return nothing
end

# nothing but the voltage violation: the generators redispatch freely inside their own limits, and
# the taps stay at their input data (only ACOPFTapFormulation makes them variables)
function objective_function!(
    input,
    model,
    ::ACOPFMinVoltageViolationFormulation,
    ::Type{ModelBuild},
)
    V_over, V_under = model[:V_over], model[:V_under]
    @objective(
        model,
        Min,
        sum(V_over[i] + V_under[i] for i in variable_tap_buses(input); init = 0.0)
    )
    return nothing
end

# nothing but the reactive balance violation at the terminals of the variable-tap transformers. As in
# the voltage violation formulation, the generators redispatch freely inside their own limits and the
# taps stay at their input data.
function objective_function!(
    input,
    model,
    ::ACOPFMinReactiveBalanceViolationFormulation,
    ::Type{ModelBuild},
)
    Q_deficit, Q_surplus = model[:Q_deficit], model[:Q_surplus]
    @objective(
        model,
        Min,
        sum(Q_deficit[i] + Q_surplus[i] for i in variable_tap_buses(input); init = 0.0)
    )
    return nothing
end
