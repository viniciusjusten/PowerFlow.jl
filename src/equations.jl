# used by the NonlinearSolve.jl-based standard power flow (acpf_standard.jl), which always takes tap
# as fixed data baked into a precomputed Ybus. The JuMP model path (power_balance.jl) instead uses
# active_injection/reactive_injection below, which support tap as either fixed data or a JuMP variable.
function active_power_calc(i::Int, V, angle, G, B)
    n = length(V)
    return V[i] * sum(V[k] * (G[i, k] * cos(angle[i] - angle[k]) + B[i, k] * sin(angle[i] - angle[k])) for k in 1:n)
end

function reactive_power_calc(i::Int, V, angle, G, B)
    n = length(V)
    return V[i] * sum(V[k] * (G[i, k] * sin(angle[i] - angle[k]) - B[i, k] * cos(angle[i] - angle[k])) for k in 1:n)
end

# a branch whose tap ratio is a degree of freedom rather than fixed data: a transformer whose
# min/max tap ratios leave it some room. Used both to declare the JuMP tap variables and to pick the
# transformers that get tap sensitivities.
function branch_has_variable_tap(input, k::Int)
    return PSD.branch_is_transformer(input, k) &&
           PSD.branch_min_tap_ratio(input, k) != PSD.branch_max_tap_ratio(input, k)
end

# per-branch admittance terms shared by active_injection/reactive_injection at both sides of a branch —
# computed once per branch (not once per call) since `y = 1/(r + im*x)` would otherwise be redone
# redundantly for every one of the 4 (side, active/reactive) combinations that need it.
# `tap` is taken as an argument (rather than read from `input`) so it can be either the branch's fixed
# tap_ratio data or a JuMP tap variable, depending on the formulation.
function branch_admittance_terms(input, k::Int, tap)
    y = 1 / (PSD.branch_resistance(input, k) + im * PSD.branch_reactance(input, k))
    g, b = real(y), imag(y)
    shift = PSD.branch_phase_shift(input, k)
    cos_shift, sin_shift = cos(shift), sin(shift)
    g_sh = PSD.branch_shunt_conductance(input, k)
    b_sh = PSD.branch_shunt_susceptance(input, k)

    # Matpower/PowerModels convention (see Branch.tap_ratio's docstring): tap^2 divides the whole
    # from-bus self admittance, shunt included, and the to-bus side is unscaled.
    G_from = (g + g_sh) / tap^2
    B_from = (b + b_sh) / tap^2
    G_to = g + g_sh
    B_to = b + b_sh

    # real/imag parts of y / tap * cis(shift) and y / tap * cis(-shift), expanded term-by-term to
    # avoid Complex arithmetic over `tap` (which breaks when `tap` is a JuMP variable/expression)
    Gt_from = (g * cos_shift - b * sin_shift) / tap
    Bt_from = (g * sin_shift + b * cos_shift) / tap
    Gt_to = (g * cos_shift + b * sin_shift) / tap
    Bt_to = (b * cos_shift - g * sin_shift) / tap

    return (G_from, B_from, Gt_from, Bt_from), (G_to, B_to, Gt_to, Bt_to)
end

# active power injected into a branch at `side`, from precomputed branch_admittance_terms
function active_injection(V, angle, from::Int, to::Int, side::Symbol, terms)
    (G_from, _, Gt_from, Bt_from), (G_to, _, Gt_to, Bt_to) = terms

    if side == :from
        return V[from]^2 * G_from -
               V[from] * V[to] * (Gt_from * cos(angle[from] - angle[to]) + Bt_from * sin(angle[from] - angle[to]))
    else
        return V[to]^2 * G_to -
               V[to] * V[from] * (Gt_to * cos(angle[to] - angle[from]) + Bt_to * sin(angle[to] - angle[from]))
    end
end

# reactive power injected into a branch at `side`, from precomputed branch_admittance_terms
function reactive_injection(V, angle, from::Int, to::Int, side::Symbol, terms)
    (_, B_from, Gt_from, Bt_from), (_, B_to, Gt_to, Bt_to) = terms

    if side == :from
        return -V[from]^2 * B_from -
               V[from] * V[to] * (Gt_from * sin(angle[from] - angle[to]) - Bt_from * cos(angle[from] - angle[to]))
    else
        return -V[to]^2 * B_to -
               V[to] * V[from] * (Gt_to * sin(angle[to] - angle[from]) - Bt_to * cos(angle[to] - angle[from]))
    end
end

# convenience forms that compute branch_admittance_terms internally, for standalone/one-off use
# (e.g. tests) where reusing terms across multiple calls doesn't matter
function active_injection(input, k::Int, V, angle, side::Symbol, tap = PSD.branch_tap_ratio(input, k))
    terms = branch_admittance_terms(input, k, tap)
    return active_injection(V, angle, PSD.branch_from_bus_idx(input, k), PSD.branch_to_bus_idx(input, k), side, terms)
end

function reactive_injection(input, k::Int, V, angle, side::Symbol, tap = PSD.branch_tap_ratio(input, k))
    terms = branch_admittance_terms(input, k, tap)
    return reactive_injection(V, angle, PSD.branch_from_bus_idx(input, k), PSD.branch_to_bus_idx(input, k), side, terms)
end

function bus_active_power_injection(input, i::Int)
    total = 0.0
    for g in 1:PSD.num_generators(input)
        if PSD.generator_bus_idx(input, g) == i
            total += PSD.generator_active_power(input, g)
        end
    end
    for l in 1:PSD.num_loads(input)
        if PSD.load_bus_idx(input, l) == i
            total -= PSD.load_active_power(input, l)
        end
    end
    for k in 1:PSD.num_dclines(input)
        if PSD.dcline_from_bus_idx(input, k) == i
            total -= PSD.dcline_active_power(input, k)
        end
        if PSD.dcline_to_bus_idx(input, k) == i
            total += dcline_to_bus_active_power(input, k)
        end
    end
    return total
end

function bus_reactive_power_injection(input, i::Int)
    total = 0.0
    for g in 1:PSD.num_generators(input)
        if PSD.generator_bus_idx(input, g) == i
            total += PSD.generator_reactive_power(input, g)
        end
    end
    for l in 1:PSD.num_loads(input)
        if PSD.load_bus_idx(input, l) == i
            total -= PSD.load_reactive_power(input, l)
        end
    end
    # the DC line data is the reactive power leaving the bus into the line at each terminal
    for k in 1:PSD.num_dclines(input)
        if PSD.dcline_from_bus_idx(input, k) == i
            total -= PSD.dcline_reactive_power_from(input, k)
        end
        if PSD.dcline_to_bus_idx(input, k) == i
            total -= PSD.dcline_reactive_power_to(input, k)
        end
    end
    return total
end

# active power arriving at the to-bus terminal of a DC line, after the linear loss model
# (constant + proportional to the from-bus flow) is applied to the scheduled from-bus flow
function dcline_to_bus_active_power(input, k::Int)
    Pf = PSD.dcline_active_power(input, k)
    return Pf - (PSD.dcline_constant_loss(input, k) + PSD.dcline_loss_coefficient(input, k) * Pf)
end

# textbook DC power flow approximation: resistance, shunt elements, tap ratio and phase shift are ignored
function dc_susceptance_matrix(input)
    n = PSD.num_buses(input)
    B = zeros(n, n)

    for i in 1:PSD.num_branches(input)
        from = PSD.branch_from_bus_idx(input, i)
        to = PSD.branch_to_bus_idx(input, i)
        b = 1 / PSD.branch_reactance(input, i)

        B[from, from] += b
        B[to, to] += b
        B[from, to] -= b
        B[to, from] -= b
    end

    return B
end
