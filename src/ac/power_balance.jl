function power_balance!(input, model, ::ACFormulation, ::Type{ModelBuild})
    n = num_buses(input)
    n_g = num_generators(input)
    n_l = num_loads(input)

    Ybus = PowerSystemData.ybus(input)
    G, B = real.(Ybus), imag.(Ybus)

    V, angle, Pg, Qg = model[:V], model[:angle], model[:Pg], model[:Qg]

    @constraint(
        model,
        active_power_balance[i in 1:n],
        sum(Pg[g] for g in 1:n_g if generator_bus_idx(input, g) == i; init = 0.0) -
        sum(load_active_power(input, l) for l in 1:n_l if load_bus_idx(input, l) == i; init = 0.0) ==
        active_power_calc(i, V, angle, G, B)
    )

    @constraint(
        model,
        reactive_power_balance[i in 1:n],
        sum(Qg[g] for g in 1:n_g if generator_bus_idx(input, g) == i; init = 0.0) -
        sum(load_reactive_power(input, l) for l in 1:n_l if load_bus_idx(input, l) == i; init = 0.0) ==
        reactive_power_calc(i, V, angle, G, B)
    )

    return nothing
end

function voltage_setpoint!(input, model, formulation::ACFormulation, ::Type{ModelBuild})
    if !has_voltage_setpoint(formulation)
        return nothing
    end

    V = model[:V]
    for g in 1:num_generators(input)
        bus_idx = generator_bus_idx(input, g)
        if bus_type(input, bus_idx) in (BusType.PV, BusType.Slack)
            JuMP.fix(V[bus_idx], generator_voltage_setpoint(input, g); force = true)
        end
    end

    return nothing
end

# thermal (rated apparent power) limits are only enforced for the redispatchable OPF-style
# formulations — a plain power-flow solve (ACPFFeasibilityFormulation, or the NonlinearSolve-based
# ACPFFormulation) reports whatever flows result, it doesn't reject/redispatch around an overload
function branch_flow_limit!(input, model, formulation::ACFormulation, ::Type{ModelBuild})
    if !allows_generation_redispatch(formulation)
        return nothing
    end

    V, angle = model[:V], model[:angle]
    n_branches = num_branches(input)

    # named subexpressions (not variables — no extra columns in the KKT system) so any other
    # constraint/objective built later can reference model[:Pkm][k, side] instead of re-deriving it,
    # and JuMP/MOI's nonlinear evaluator computes each one once per iteration and reuses it everywhere
    @expression(
        model,
        Pkm[k in 1:n_branches, side in (:from, :to)],
        branch_active_power_flow(input, k, V, angle, side)
    )
    @expression(
        model,
        Qkm[k in 1:n_branches, side in (:from, :to)],
        branch_reactive_power_flow(input, k, V, angle, side)
    )

    for k in 1:n_branches
        max_flow = branch_max_flow(input, k)
        if !isfinite(max_flow)
            continue
        end

        for side in (:from, :to)
            @constraint(model, Pkm[k, side]^2 + Qkm[k, side]^2 <= max_flow^2)
        end
    end

    return nothing
end
