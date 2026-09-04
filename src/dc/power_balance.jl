function power_balance!(input, model, ::Type{ModelBuild})
    n = num_buses(input)
    n_g = num_generators(input)
    n_l = num_loads(input)
    n_branches = num_branches(input)

    B = dc_susceptance_matrix(input)
    angle, Pg = model[:angle], model[:Pg]

    @expression(
        model,
        flow[k in 1:n_branches],
        (angle[branch_from_bus_idx(input, k)] - angle[branch_to_bus_idx(input, k)]) / branch_reactance(input, k)
    )

    @constraint(
        model,
        active_power_balance[i in 1:n],
        sum(Pg[g] for g in 1:n_g if generator_bus_idx(input, g) == i; init = 0.0) -
        sum(load_active_power(input, l) for l in 1:n_l if load_bus_idx(input, l) == i; init = 0.0) ==
        sum(B[i, j] * (angle[i] - angle[j]) for j in 1:n if j != i)
    )

    for k in 1:n_branches
        max_flow = branch_max_flow(input, k)
        if isfinite(max_flow)
            @constraint(model, flow[k] <= max_flow)
            @constraint(model, flow[k] >= -max_flow)
        end
    end

    return nothing
end
