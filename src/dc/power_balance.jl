function power_balance!(input, model, ::Type{ModelBuild})
    n = PSD.num_buses(input)
    n_g = PSD.num_generators(input)
    n_l = PSD.num_loads(input)
    n_branches = PSD.num_branches(input)
    n_dc = PSD.num_dclines(input)

    B = dc_susceptance_matrix(input)
    angle, Pg = model[:angle], model[:Pg]
    Pf, Pt = model[:Pf], model[:Pt]

    @expression(
        model,
        flow[k in 1:n_branches],
        (angle[PSD.branch_from_bus_idx(input, k)] - angle[PSD.branch_to_bus_idx(input, k)]) /
        PSD.branch_reactance(input, k)
    )

    @constraint(
        model,
        active_power_balance[i in 1:n],
        sum(Pg[g] for g in 1:n_g if PSD.generator_bus_idx(input, g) == i; init = 0.0) -
        sum(PSD.load_active_power(input, l) for l in 1:n_l if PSD.load_bus_idx(input, l) == i; init = 0.0) -
        sum(Pf[k] for k in 1:n_dc if PSD.dcline_from_bus_idx(input, k) == i; init = 0.0) -
        sum(Pt[k] for k in 1:n_dc if PSD.dcline_to_bus_idx(input, k) == i; init = 0.0) ==
        sum(B[i, j] * (angle[i] - angle[j]) for j in 1:n if j != i)
    )

    for k in 1:n_branches
        max_flow = PSD.branch_max_flow(input, k)
        if isfinite(max_flow)
            @constraint(model, flow[k] <= max_flow)
            @constraint(model, flow[k] >= -max_flow)
        end
    end

    return nothing
end
