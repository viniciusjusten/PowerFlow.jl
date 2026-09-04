function objective_function!(input, model, ::ACOPFFormulation, ::Type{ModelBuild})
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
