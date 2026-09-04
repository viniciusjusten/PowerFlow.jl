# TODO: Ybus (and therefore G, B here) already bakes in each branch's tap_ratio/phase_shift as a
# fixed constant from the data. That's fine as long as taps are given, fixed inputs. If we ever want
# to optimize the tap ratio itself (e.g. OLTC control), these formulas won't work as-is. Ybus would
# need to be rebuilt from a JuMP tap variable inside the model instead of being precomputed once from
# PowerSystemData. Revisit in a future discussion.
function active_power_calc(i::Int, V, angle, G, B)
    n = length(V)
    return V[i] * sum(V[k] * (G[i, k] * cos(angle[i] - angle[k]) + B[i, k] * sin(angle[i] - angle[k])) for k in 1:n)
end

function reactive_power_calc(i::Int, V, angle, G, B)
    n = length(V)
    return V[i] * sum(V[k] * (G[i, k] * sin(angle[i] - angle[k]) - B[i, k] * cos(angle[i] - angle[k])) for k in 1:n)
end

# per-branch admittance terms feeding the branch flow equations below — kept separate from the
# aggregate Ybus so branch flows are correct even with parallel branches between the same bus pair
function branch_flow_terms(input, k::Int)
    y = 1 / (branch_resistance(input, k) + im * branch_reactance(input, k))
    tap = branch_tap_ratio(input, k)
    shift = branch_phase_shift(input, k)
    g_sh = branch_shunt_conductance(input, k)
    b_sh = branch_shunt_susceptance(input, k)

    y_from = y * tap * cis(shift)
    y_to = y * tap * cis(-shift)

    G_from = real(y) * tap^2 + g_sh
    B_from = imag(y) * tap^2 + b_sh
    G_to = real(y) + g_sh
    B_to = imag(y) + b_sh

    return (G_from, B_from, real(y_from), imag(y_from)), (G_to, B_to, real(y_to), imag(y_to))
end

function branch_active_power_flow(input, k::Int, V, angle, side::Symbol)
    from = branch_from_bus_idx(input, k)
    to = branch_to_bus_idx(input, k)
    (G_from, _, Gt_from, Bt_from), (G_to, _, Gt_to, Bt_to) = branch_flow_terms(input, k)

    if side == :from
        return V[from]^2 * G_from -
               V[from] * V[to] * (Gt_from * cos(angle[from] - angle[to]) + Bt_from * sin(angle[from] - angle[to]))
    else
        return V[to]^2 * G_to -
               V[to] * V[from] * (Gt_to * cos(angle[to] - angle[from]) + Bt_to * sin(angle[to] - angle[from]))
    end
end

function branch_reactive_power_flow(input, k::Int, V, angle, side::Symbol)
    from = branch_from_bus_idx(input, k)
    to = branch_to_bus_idx(input, k)
    (_, B_from, Gt_from, Bt_from), (_, B_to, Gt_to, Bt_to) = branch_flow_terms(input, k)

    if side == :from
        return -V[from]^2 * B_from -
               V[from] * V[to] * (Gt_from * sin(angle[from] - angle[to]) - Bt_from * cos(angle[from] - angle[to]))
    else
        return -V[to]^2 * B_to -
               V[to] * V[from] * (Gt_to * sin(angle[to] - angle[from]) - Bt_to * cos(angle[to] - angle[from]))
    end
end

function bus_active_power_injection(input, i::Int)
    total = 0.0
    for g in 1:num_generators(input)
        if generator_bus_idx(input, g) == i
            total += generator_active_power(input, g)
        end
    end
    for l in 1:num_loads(input)
        if load_bus_idx(input, l) == i
            total -= load_active_power(input, l)
        end
    end
    return total
end

function bus_reactive_power_injection(input, i::Int)
    total = 0.0
    for g in 1:num_generators(input)
        if generator_bus_idx(input, g) == i
            total += generator_reactive_power(input, g)
        end
    end
    for l in 1:num_loads(input)
        if load_bus_idx(input, l) == i
            total -= load_reactive_power(input, l)
        end
    end
    return total
end

# textbook DC power flow approximation: resistance, shunt elements, tap ratio and phase shift are ignored
function dc_susceptance_matrix(input)
    n = num_buses(input)
    B = zeros(n, n)

    for i in 1:num_branches(input)
        from = branch_from_bus_idx(input, i)
        to = branch_to_bus_idx(input, i)
        b = 1 / branch_reactance(input, i)

        B[from, from] += b
        B[to, to] += b
        B[from, to] -= b
        B[to, from] -= b
    end

    return B
end
