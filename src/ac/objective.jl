function objective_function!(input, model, ::Union{ACOPFFormulation, ACOPFTapFormulation}, ::Type{ModelBuild})
    Pg = model[:Pg]
    @objective(model, Min, sum(generator_cost(input, g) * Pg[g] for g in 1:num_generators(input)))
    return nothing
end

function objective_function!(input, model, ::ACOPFMinRedispatchFormulation, ::Type{ModelBuild})
    Pg = model[:Pg]
    @objective(
        model,
        Min,
        sum((Pg[g] - generator_active_power(input, g))^2 for g in 1:num_generators(input))
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
            (Pg[g] - generator_active_power(input, g))^2 +
            (V[generator_bus_idx(input, g)] - generator_voltage_setpoint(input, g))^2 for
            g in 1:num_generators(input)
        )
    )
    return nothing
end

# nothing but the voltage violation: the generators redispatch freely inside their own limits, and
# the taps stay at their input data (only ACOPFTapFormulation makes them variables)
function objective_function!(
    input,
    model,
    formulation::ACOPFMinVoltageViolationFormulation,
    ::Type{ModelBuild},
)
    V_over, V_under = model[:V_over], model[:V_under]
    @objective(
        model,
        Min,
        sum(V_over[i] + V_under[i] for i in voltage_violation_buses(input, formulation); init = 0.0)
    )
    return nothing
end
