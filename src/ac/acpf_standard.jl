struct ACPFStandardProblem
    n::Int
    G::Matrix{Float64}
    B::Matrix{Float64}
    slack_idx::Int
    slack_angle::Float64
    non_slack_buses::Vector{Int}
    pq_buses::Vector{Int}
    voltage_setpoint::Vector{Float64} # meaningful only at PV/Slack buses
    Pspec::Vector{Float64} # fixed generator data minus load, per bus
    Qspec::Vector{Float64}
end

function build_acpf_standard_problem(input)
    n = num_buses(input)
    Ybus = PowerSystemData.ybus(input)
    G, B = real.(Ybus), imag.(Ybus)

    slack_idx = findfirst(i -> bus_type(input, i) == BusType.Slack, 1:n)
    non_slack_buses = [i for i in 1:n if i != slack_idx]
    pq_buses = [i for i in 1:n if bus_type(input, i) == BusType.PQ]

    voltage_setpoint = [bus_voltage(input, i) for i in 1:n]
    for g in 1:num_generators(input)
        bus_idx = generator_bus_idx(input, g)
        if bus_type(input, bus_idx) in (BusType.PV, BusType.Slack)
            voltage_setpoint[bus_idx] = generator_voltage_setpoint(input, g)
        end
    end

    Pspec = [bus_active_power_injection(input, i) for i in 1:n]
    Qspec = [bus_reactive_power_injection(input, i) for i in 1:n]

    return ACPFStandardProblem(
        n,
        G,
        B,
        slack_idx,
        bus_voltage_angle(input, slack_idx),
        non_slack_buses,
        pq_buses,
        voltage_setpoint,
        Pspec,
        Qspec,
    )
end

function unpack_state(x, p::ACPFStandardProblem)
    T = eltype(x)
    angle = zeros(T, p.n)
    V = Vector{T}(p.voltage_setpoint)

    angle[p.slack_idx] = p.slack_angle
    n_non_slack = length(p.non_slack_buses)
    for (k, i) in enumerate(p.non_slack_buses)
        angle[i] = x[k]
    end
    for (k, i) in enumerate(p.pq_buses)
        V[i] = x[n_non_slack + k]
    end

    return V, angle
end

function power_flow_residual!(res, x, p::ACPFStandardProblem)
    V, angle = unpack_state(x, p)
    n_non_slack = length(p.non_slack_buses)

    for (k, i) in enumerate(p.non_slack_buses)
        res[k] = active_power_calc(i, V, angle, p.G, p.B) - p.Pspec[i]
    end
    for (k, i) in enumerate(p.pq_buses)
        res[n_non_slack + k] = reactive_power_calc(i, V, angle, p.G, p.B) - p.Qspec[i]
    end

    return nothing
end

# assumes at most one generator per PV/Slack bus — a documented simplification for now
function generator_reactive_power_output(input, g::Int, V, angle, G, B)
    bus_idx = generator_bus_idx(input, g)
    if bus_type(input, bus_idx) == BusType.PQ
        return generator_reactive_power(input, g)
    end

    load_q =
        sum(load_reactive_power(input, l) for l in 1:num_loads(input) if load_bus_idx(input, l) == bus_idx; init = 0.0)
    return reactive_power_calc(bus_idx, V, angle, G, B) + load_q
end

function run_acpf(input, ::ACPFFormulation)
    p = build_acpf_standard_problem(input)
    n_non_slack = length(p.non_slack_buses)
    n_pq = length(p.pq_buses)

    x0 = zeros(n_non_slack + n_pq)
    for (k, i) in enumerate(p.non_slack_buses)
        x0[k] = bus_voltage_angle(input, i)
    end
    for (k, i) in enumerate(p.pq_buses)
        x0[n_non_slack + k] = bus_voltage(input, i)
    end

    f = NonlinearSolve.NonlinearFunction(
        power_flow_residual!;
        sparsity = SparseConnectivityTracer.TracerSparsityDetector(),
    )
    prob = NonlinearSolve.NonlinearProblem(f, x0, p)
    sol = NonlinearSolve.solve(prob, NonlinearSolve.NewtonRaphson(; linsolve = LinearSolve.KrylovJL_GMRES())) # , show_trace = Val(true), trace_level = TraceAll(), store_trace = Val(true)

    V, angle = unpack_state(sol.u, p)
    Pg = [generator_active_power(input, g) for g in 1:num_generators(input)]
    Qg = [generator_reactive_power_output(input, g, V, angle, p.G, p.B) for g in 1:num_generators(input)]

    P_from, P_to, Q_from, Q_to = ac_branch_flows(input, V, angle)
    Pf, Pt, Qf, Qt = dcline_scheduled_flows(input)

    return ACSolution(V, angle, Pg, Qg, P_from, P_to, Q_from, Q_to, Pf, Pt, Qf, Qt, NaN, sol.retcode)
end
