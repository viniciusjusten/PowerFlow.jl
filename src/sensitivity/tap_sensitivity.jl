"""
    TapSensitivity

Sensitivities of one variable-tap transformer, referenced by `branch_idx` (the branch that owns the
tap) and by the pair of buses its terminals connect. All derivatives are taken at a solved operating
point; `from`/`to` in the field names always refer to that branch's terminals.

The fields are filled in stages, each by its own function, so a partially filled struct is the normal
intermediate state:

  - [`build_tap_sensitivities`](@ref) fills the operating point and the dQ/dV block,
  - [`reactive_power_tap_derivatives!`](@ref) fills `dQ_from_dtap` / `dQ_to_dtap`,
  - [`voltage_tap_derivatives!`](@ref) fills `dV_from_dtap` / `dV_to_dtap`.

Fields that have not been filled yet hold `NaN`.
"""
Base.@kwdef mutable struct TapSensitivity
    branch_idx::Int
    branch_name::String
    from_bus_idx::Int
    to_bus_idx::Int
    tap_ratio::Float64
    # operating point the derivatives are taken at
    voltage_from::Float64
    voltage_to::Float64
    angle_difference::Float64 # rad, angle[from] - angle[to]
    # bus reactive power injection vs. bus voltage magnitude, for the four combinations of the
    # transformer's two terminals: dQ_<numerator bus>_dV_<denominator bus>
    dQ_from_dV_from::Float64 = NaN
    dQ_from_dV_to::Float64 = NaN
    dQ_to_dV_from::Float64 = NaN
    dQ_to_dV_to::Float64 = NaN
    # bus reactive power injection vs. the branch's tap ratio
    dQ_from_dtap::Float64 = NaN
    dQ_to_dtap::Float64 = NaN
    # bus voltage magnitude vs. the branch's tap ratio, from the chain rule over the two above
    dV_from_dtap::Float64 = NaN
    dV_to_dtap::Float64 = NaN
end

# the dQ/dV block as a 2x2 matrix, rows indexed by the bus whose reactive injection varies and
# columns by the bus whose voltage magnitude varies, both ordered (from, to)
function dQ_dV_matrix(s::TapSensitivity)
    return [s.dQ_from_dV_from s.dQ_from_dV_to; s.dQ_to_dV_from s.dQ_to_dV_to]
end

dQ_dtap_vector(s::TapSensitivity) = [s.dQ_from_dtap, s.dQ_to_dtap]

# d(bus i reactive power injection)/d(bus j voltage magnitude), differentiating the injection
# Q_i = V_i * sum_k V_k * (G_ik * sin(angle_ik) - B_ik * cos(angle_ik)) used by reactive_power_calc,
# holding every angle and every other voltage magnitude fixed
function reactive_power_voltage_derivative(i::Int, j::Int, V, angle, G, B)
    if i == j
        return reactive_power_calc(i, V, angle, G, B) / V[i] - B[i, i] * V[i]
    end
    return V[i] * (G[i, j] * sin(angle[i] - angle[j]) - B[i, j] * cos(angle[i] - angle[j]))
end

"""
    build_tap_sensitivities(input, solution::ACSolution)

Build one [`TapSensitivity`](@ref) per variable-tap transformer of `input`, with the operating point
and the dQ/dV block filled in. Nothing is solved here: the voltages and angles come from `solution`
and `input` only supplies the network model, so `input` must already carry the tap ratios `solution`
was obtained at.

Only the four dQ/dV combinations over the transformer's own two terminals are computed — the rest of
the network's reaction to a voltage change is not included.
"""
function build_tap_sensitivities(input, solution::ACSolution)
    problem = build_acpf_standard_problem(input)
    G, B = problem.G, problem.B
    V, angle = solution.voltage, solution.voltage_angle

    sensitivities = TapSensitivity[]
    for k in 1:num_branches(input)
        if !branch_has_variable_tap(input, k)
            continue
        end

        from = branch_from_bus_idx(input, k)
        to = branch_to_bus_idx(input, k)

        push!(
            sensitivities,
            TapSensitivity(;
                branch_idx = k,
                branch_name = branch_name(input, k),
                from_bus_idx = from,
                to_bus_idx = to,
                tap_ratio = branch_tap_ratio(input, k),
                voltage_from = V[from],
                voltage_to = V[to],
                angle_difference = angle[from] - angle[to],
                dQ_from_dV_from = reactive_power_voltage_derivative(from, from, V, angle, G, B),
                dQ_from_dV_to = reactive_power_voltage_derivative(from, to, V, angle, G, B),
                dQ_to_dV_from = reactive_power_voltage_derivative(to, from, V, angle, G, B),
                dQ_to_dV_to = reactive_power_voltage_derivative(to, to, V, angle, G, B),
            ),
        )
    end

    return sensitivities
end

"""
    reactive_power_tap_derivatives!(sensitivities, input)

Fill `dQ_from_dtap` / `dQ_to_dtap` on each element of `sensitivities` with the analytic derivative of
the branch's reactive power injection with respect to its tap ratio, evaluated at the operating point
already stored on the struct. Returns `sensitivities`.

With `g_s + im * b_s` the series admittance, `b_c / 2` the branch's per-side charging susceptance,
`tap` the tap ratio and `dtheta = angle_from - angle_to`:

    dQ_from/dtap = 2 * (V_from^2 / tap^3) * (b_s + b_c / 2)
                   + (V_from * V_to / tap^2) * (g_s * sin(dtheta) - b_s * cos(dtheta))
    dQ_to/dtap   = -(V_from * V_to / tap^2) * (g_s * sin(dtheta) + b_s * cos(dtheta))

The `b_c / 2` inside the first term comes from the Matpower form of the from-side self admittance,
`(y_s + im * b_c / 2) / tap^2`, which is what `branch_admittance_terms` and `PowerSystemData.ybus`
build — so these are the exact derivatives of the network model being solved here. They do assume a
zero phase shift, the only case `branch_admittance_terms` and these formulas agree on.
"""
function reactive_power_tap_derivatives!(sensitivities, input)
    for s in sensitivities
        k = s.branch_idx
        y = 1 / (branch_resistance(input, k) + im * branch_reactance(input, k))
        g_s, b_s = real(y), imag(y)
        b_charging = branch_shunt_susceptance(input, k) # already the per-side b_c / 2

        tap = s.tap_ratio
        V_from, V_to = s.voltage_from, s.voltage_to
        sin_diff, cos_diff = sin(s.angle_difference), cos(s.angle_difference)

        s.dQ_from_dtap =
            2 * (V_from^2 / tap^3) * (b_s + b_charging) +
            (V_from * V_to / tap^2) * (g_s * sin_diff - b_s * cos_diff)
        s.dQ_to_dtap = -(V_from * V_to / tap^2) * (g_s * sin_diff + b_s * cos_diff)
    end

    return sensitivities
end

"""
    voltage_tap_derivatives!(sensitivities, input)

Fill `dV_from_dtap` / `dV_to_dtap` on each element of `sensitivities` by chaining the dQ/dV block and
the dQ/dtap vector already stored on the struct. Returns `sensitivities`.

The specified reactive injections at the two terminals do not move with the tap, so the total
derivative of each injection vanishes:

    dQ/dtap + (dQ/dV) * (dV/dtap) = 0   =>   dV/dtap = -(dQ/dV) \\ (dQ/dtap)

restricted to the terminals whose voltage is actually free. A PV or Slack terminal is held at its
voltage setpoint, so its dV/dtap is zero and it is dropped from the linear system — there the tap
moves the terminal's reactive generation instead of its voltage. Requires
[`reactive_power_tap_derivatives!`](@ref) to have run first.
"""
function voltage_tap_derivatives!(sensitivities, input)
    for s in sensitivities
        terminals = (s.from_bus_idx, s.to_bus_idx)
        free = [pos for (pos, bus) in enumerate(terminals) if bus_type(input, bus) == BusType.PQ]

        dV_dtap = zeros(2)
        if !isempty(free)
            dV_dtap[free] = -(dQ_dV_matrix(s)[free, free] \ dQ_dtap_vector(s)[free])
        end

        s.dV_from_dtap, s.dV_to_dtap = dV_dtap[1], dV_dtap[2]
    end

    return sensitivities
end

"""
    run_tap_sensitivities(input)
    run_tap_sensitivities(input, solution::ACSolution)

Full pipeline for the variable-tap transformers of a system: build the dQ/dV block at a solved
operating point, add the analytic dQ/dtap derivatives, then chain the two into dV/dtap. Returns the
vector of [`TapSensitivity`](@ref).

Everything is analytic from `solution` onwards, so the two-argument form solves nothing — pass the
solution of an already run case together with the inputs it was run on, taps included. The
single-argument form has no solution to start from and runs the standard power flow itself.
"""
function run_tap_sensitivities(input, solution::ACSolution)
    sensitivities = build_tap_sensitivities(input, solution)
    reactive_power_tap_derivatives!(sensitivities, input)
    voltage_tap_derivatives!(sensitivities, input)
    return sensitivities
end

function run_tap_sensitivities(input)
    return run_tap_sensitivities(input, run_acpf(input, ACPFFormulation()))
end
