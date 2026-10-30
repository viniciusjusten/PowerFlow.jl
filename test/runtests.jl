using PowerSystemData
using PowerFlow
using JuMP
using Ipopt
using HiGHS
using Test

function build_test_system(;
    l13_max_flow::Float64 = 5.0,
    l12_min_tap_ratio::Float64 = 1.0,
    l12_max_tap_ratio::Float64 = 1.0,
    l23_is_transformer::Bool = false,
    l23_shunt_susceptance::Float64 = 0.0,
    l23_tap_ratio::Float64 = 1.0,
    l23_min_tap_ratio::Float64 = l23_tap_ratio,
    l23_max_tap_ratio::Float64 = l23_tap_ratio,
    pq_min_voltage::Float64 = 0.9,
    pq_max_voltage::Float64 = 1.1,
)
    input = PowerSystemData.Input()

    add_bus!(input; name = "Slack", type = BusType.Slack, voltage = 1.0, min_voltage = 0.9, max_voltage = 1.1)
    add_bus!(input; name = "PV", type = BusType.PV, voltage = 1.0, min_voltage = 0.9, max_voltage = 1.1)
    add_bus!(input; name = "PQ", type = BusType.PQ, min_voltage = pq_min_voltage, max_voltage = pq_max_voltage)

    add_branch!(
        input;
        name = "L12",
        from_bus_idx = 1,
        to_bus_idx = 2,
        resistance = 0.01,
        reactance = 0.1,
        max_flow = 5.0,
        min_tap_ratio = l12_min_tap_ratio,
        max_tap_ratio = l12_max_tap_ratio,
    )
    add_branch!(
        input;
        name = "L23",
        from_bus_idx = 2,
        to_bus_idx = 3,
        resistance = 0.01,
        reactance = 0.1,
        max_flow = 5.0,
        shunt_susceptance = l23_shunt_susceptance,
        is_transformer = l23_is_transformer,
        tap_ratio = l23_tap_ratio,
        min_tap_ratio = l23_min_tap_ratio,
        max_tap_ratio = l23_max_tap_ratio,
    )
    add_branch!(
        input;
        name = "L13",
        from_bus_idx = 1,
        to_bus_idx = 3,
        resistance = 0.02,
        reactance = 0.2,
        max_flow = l13_max_flow,
    )

    add_generator!(
        input;
        name = "G1",
        bus_idx = 1,
        active_power = 0.5,
        reactive_power = 0.1,
        voltage_setpoint = 1.0,
        min_active_power = 0.0,
        max_active_power = 2.0,
        min_reactive_power = -1.0,
        max_reactive_power = 1.0,
        cost = 20.0,
    )
    add_generator!(
        input;
        name = "G2",
        bus_idx = 2,
        active_power = 0.5,
        reactive_power = 0.1,
        voltage_setpoint = 1.02,
        min_active_power = 0.0,
        max_active_power = 2.0,
        min_reactive_power = -1.0,
        max_reactive_power = 1.0,
        cost = 25.0,
    )

    add_load!(input; name = "D3", bus_idx = 3, active_power = 0.9, reactive_power = 0.3)

    return input
end

function add_test_dcline!(input)
    add_dcline!(
        input;
        name = "DC13",
        from_bus_idx = 1,
        to_bus_idx = 3,
        active_power = 0.2,
        reactive_power_from = 0.05,
        reactive_power_to = -0.03,
        constant_loss = 0.01,
        loss_coefficient = 0.02,
        max_flow_from_to = 0.3,
        max_flow_to_from = 0.25,
        min_reactive_power_from = -0.1,
        max_reactive_power_from = 0.1,
        min_reactive_power_to = -0.1,
        max_reactive_power_to = 0.1,
    )
    return input
end

function reactive_balance_mismatch(input, sol, i)
    generation = sum(
        sol.reactive_power_generation[g] for g in 1:num_generators(input) if generator_bus_idx(input, g) == i;
        init = 0.0,
    )
    load = sum(load_reactive_power(input, l) for l in 1:num_loads(input) if load_bus_idx(input, l) == i; init = 0.0)
    outflow = sum(
        sol.branch_reactive_power_from[k] for k in 1:num_branches(input) if branch_from_bus_idx(input, k) == i;
        init = 0.0,
    )
    outflow += sum(
        sol.branch_reactive_power_to[k] for k in 1:num_branches(input) if branch_to_bus_idx(input, k) == i;
        init = 0.0,
    )
    outflow -= bus_shunt_susceptance(input, i) * sol.voltage[i]^2
    return outflow - (generation - load)
end

const SOLVED = JuMP.MOI.LOCALLY_SOLVED

@testset verbose = true "PowerFlow" begin
    @testset "ACOPF" begin
        input = build_test_system()
        sol = run_acopf(input)

        @test sol.termination_status == SOLVED
        # G1 is cheaper than G2, so it should cover essentially all of the load + losses
        @test sol.active_power_generation[1] > sol.active_power_generation[2]
        @test sum(sol.active_power_generation) ≈ 0.9 atol = 0.05
    end

    @testset "ACOPF min-redispatch" begin
        input = build_test_system()
        sol = run_acopf(input, ACOPFMinRedispatchFormulation())
        P_star = [generator_active_power(input, g) for g in 1:num_generators(input)]

        @test sol.termination_status == SOLVED
        @test all(abs.(sol.active_power_generation .- P_star) .< 0.1)
    end

    @testset "ACOPF with tap optimization" begin
        input = build_test_system(; l12_min_tap_ratio = 0.9, l12_max_tap_ratio = 1.1)
        sol = run_acopf(input, ACOPFTapFormulation())

        model = JuMP.Model(Ipopt.Optimizer)
        JuMP.set_silent(model)
        PowerFlow.ac_model!(input, model, ACOPFTapFormulation(), PowerFlow.ModelBuild)
        JuMP.optimize!(model)

        tap = JuMP.value(model[:tap][1])

        @test sol.termination_status == SOLVED
        @test 0.9 - 1e-6 <= tap <= 1.1 + 1e-6
    end

    @testset "ACOPF respects branch thermal limits" begin
        input = build_test_system(; l13_max_flow = 0.3)
        sol = run_acopf(input)

        model = JuMP.Model(Ipopt.Optimizer)
        JuMP.set_silent(model)
        PowerFlow.ac_model!(input, model, ACOPFFormulation(), PowerFlow.ModelBuild)
        JuMP.optimize!(model)

        P_inj_from, P_inj_to = model[:P_inj_from], model[:P_inj_to]
        Q_inj_from, Q_inj_to = model[:Q_inj_from], model[:Q_inj_to]
        S_from = sqrt(JuMP.value(P_inj_from[3])^2 + JuMP.value(Q_inj_from[3])^2)
        S_to = sqrt(JuMP.value(P_inj_to[3])^2 + JuMP.value(Q_inj_to[3])^2)

        @test sol.termination_status == SOLVED
        @test S_from <= 0.3 + 1e-4
        @test S_to <= 0.3 + 1e-4
        # P_inj/Q_inj are registered for every branch, not just rated ones, so they can be reused elsewhere
        @test length(P_inj_from) == 3
    end

    @testset "ACPF feasibility matches standard NonlinearSolve solve" begin
        input = build_test_system()

        sol_feasibility = run_acpf(input, ACPFFeasibilityFormulation())
        sol_standard = run_acpf(input)

        @test sol_feasibility.termination_status == SOLVED
        @test sol_feasibility.objective_value ≈ 0.0 atol = 1e-6
        # the PV generator's dispatch and voltage setpoint are fixed data in this formulation
        @test sol_feasibility.active_power_generation[2] ≈ generator_active_power(input, 2) atol = 1e-6
        @test sol_feasibility.voltage[1] ≈ 1.0 atol = 1e-6
        @test sol_feasibility.voltage[2] ≈ 1.02 atol = 1e-6

        @test sol_standard.voltage ≈ sol_feasibility.voltage atol = 1e-4
        @test sol_standard.voltage_angle ≈ sol_feasibility.voltage_angle atol = 1e-4
    end

    @testset "ACPF feasibility ignores branch thermal limits (plain power flow)" begin
        input = build_test_system(; l13_max_flow = 0.3)
        sol = run_acpf(input, ACPFFeasibilityFormulation())

        P_from = PowerFlow.active_injection(input, 3, sol.voltage, sol.voltage_angle, :from)
        Q_from = PowerFlow.reactive_injection(input, 3, sol.voltage, sol.voltage_angle, :from)
        S_from = sqrt(P_from^2 + Q_from^2)

        @test sol.termination_status == SOLVED
        @test S_from > 0.3 # the rating is exceeded because this formulation never redispatches around it
    end

    @testset "ACPF PmisVmis" begin
        input = build_test_system()
        sol = run_acpf(input, ACPFPmisVmisFormulation())

        @test sol.termination_status == SOLVED
    end

    @testset "DCOPF" begin
        input = build_test_system()
        sol = run_dcopf(input)

        @test sol.termination_status == JuMP.MOI.OPTIMAL
        # G1 is cheaper than G2, so DCOPF should dispatch it first
        @test sol.active_power_generation[1] > sol.active_power_generation[2]
    end

    @testset "DCPF" begin
        input = build_test_system()
        sol = run_dcpf(input)

        B = PowerFlow.dc_susceptance_matrix(input)
        non_slack = [2, 3]
        P = [PowerFlow.bus_active_power_injection(input, i) for i in non_slack]
        residual = B[non_slack, non_slack] * sol.voltage_angle[non_slack] - P

        @test all(abs.(residual) .< 1e-8)
        @test sol.active_power_generation == [generator_active_power(input, g) for g in 1:num_generators(input)]
    end

    @testset "ACOPF with DC line" begin
        input = add_test_dcline!(build_test_system())
        sol = run_acopf(input)

        model = JuMP.Model(Ipopt.Optimizer)
        JuMP.set_silent(model)
        PowerFlow.ac_model!(input, model, ACOPFFormulation(), PowerFlow.ModelBuild)
        JuMP.optimize!(model)

        Pf, Pt = JuMP.value(model[:Pf][1]), JuMP.value(model[:Pt][1])
        Qf, Qt = JuMP.value(model[:Qf][1]), JuMP.value(model[:Qt][1])

        @test sol.termination_status == SOLVED
        @test -Pt ≈ Pf * (1 - dcline_loss_coefficient(input, 1)) - dcline_constant_loss(input, 1) atol = 1e-6
        @test -dcline_max_flow_to_from(input, 1) - 1e-6 <= Pf <= dcline_max_flow_from_to(input, 1) + 1e-6
        @test dcline_min_reactive_power_from(input, 1) - 1e-6 <= Qf <= dcline_max_reactive_power_from(input, 1) + 1e-6
        @test dcline_min_reactive_power_to(input, 1) - 1e-6 <= Qt <= dcline_max_reactive_power_to(input, 1) + 1e-6
    end

    @testset "ACPF feasibility fixes the DC line to its scheduled data" begin
        input = add_test_dcline!(build_test_system())

        model = JuMP.Model(Ipopt.Optimizer)
        JuMP.set_silent(model)
        PowerFlow.ac_model!(input, model, ACPFFeasibilityFormulation(), PowerFlow.ModelBuild)
        JuMP.optimize!(model)

        @test JuMP.termination_status(model) == SOLVED
        @test JuMP.value(model[:Pf][1]) ≈ dcline_active_power(input, 1) atol = 1e-6
        @test JuMP.value(model[:Qf][1]) ≈ dcline_reactive_power_from(input, 1) atol = 1e-6
        @test JuMP.value(model[:Qt][1]) ≈ dcline_reactive_power_to(input, 1) atol = 1e-6
    end

    @testset "ACPF feasibility with a DC line matches standard NonlinearSolve solve" begin
        # the scheduled reactive power of the test DC line is nonzero at both terminals, so a sign
        # mismatch between the JuMP model and bus_reactive_power_injection would move the voltages
        input = add_test_dcline!(build_test_system())

        sol_feasibility = run_acpf(input, ACPFFeasibilityFormulation())
        sol_standard = run_acpf(input)

        @test sol_feasibility.termination_status == SOLVED
        @test sol_standard.voltage ≈ sol_feasibility.voltage atol = 1e-4
        @test sol_standard.voltage_angle ≈ sol_feasibility.voltage_angle atol = 1e-4
        @test sol_standard.branch_active_power_from ≈ sol_feasibility.branch_active_power_from atol = 1e-4
        @test sol_standard.branch_reactive_power_to ≈ sol_feasibility.branch_reactive_power_to atol = 1e-4
        @test sol_standard.dcline_active_power_to ≈ sol_feasibility.dcline_active_power_to atol = 1e-6
        @test sol_standard.dcline_reactive_power_from ≈ sol_feasibility.dcline_reactive_power_from atol = 1e-6
    end

    @testset "DCOPF with DC line" begin
        input = add_test_dcline!(build_test_system())
        sol = run_dcopf(input)

        model = JuMP.Model(HiGHS.Optimizer)
        JuMP.set_silent(model)
        PowerFlow.dc_model!(input, model, PowerFlow.ModelBuild)
        JuMP.optimize!(model)

        Pf, Pt = JuMP.value(model[:Pf][1]), JuMP.value(model[:Pt][1])

        @test sol.termination_status == JuMP.MOI.OPTIMAL
        @test -Pt ≈ Pf * (1 - dcline_loss_coefficient(input, 1)) - dcline_constant_loss(input, 1) atol = 1e-6
        @test -dcline_max_flow_to_from(input, 1) - 1e-6 <= Pf <= dcline_max_flow_from_to(input, 1) + 1e-6
    end

    @testset "DCPF picks up the DC line's scheduled flow" begin
        input = add_test_dcline!(build_test_system())
        sol = run_dcpf(input)

        B = PowerFlow.dc_susceptance_matrix(input)
        non_slack = [2, 3]
        P = [PowerFlow.bus_active_power_injection(input, i) for i in non_slack]
        residual = B[non_slack, non_slack] * sol.voltage_angle[non_slack] - P

        @test all(abs.(residual) .< 1e-8)
        @test PowerFlow.dcline_to_bus_active_power(input, 1) ≈
              dcline_active_power(input, 1) * (1 - dcline_loss_coefficient(input, 1)) -
              dcline_constant_loss(input, 1) atol = 1e-10
    end

    @testset "ACOPF minimizing the voltage violation at variable-tap terminals" begin
        # L23 is the only variable-tap transformer, so buses 2 and 3 are the penalized terminals.
        # Bus 3 is asked for a band above what buses 1 and 2 can hold it at, so the violation cannot
        # reach zero.
        input = build_test_system(;
            l23_is_transformer = true,
            l23_min_tap_ratio = 0.9,
            l23_max_tap_ratio = 1.1,
            pq_min_voltage = 1.15,
            pq_max_voltage = 1.2,
        )
        sol = run_acopf(input, ACOPFMinVoltageViolationFormulation())

        @test sol.termination_status == SOLVED

        # the objective is the violation itself, summed over the two terminals of L23 only
        violation(i) =
            max(0.0, sol.voltage[i] - bus_max_voltage(input, i)) +
            max(0.0, bus_min_voltage(input, i) - sol.voltage[i])
        @test sol.objective_value > 1e-3
        @test sol.objective_value ≈ violation(2) + violation(3) atol = 1e-6
        @test violation(3) > 1e-3

        # bus 1 is not a tap terminal, so its limits stayed hard bounds
        @test 0.9 - 1e-6 <= sol.voltage[1] <= 1.1 + 1e-6
        # the taps are data, not variables: minimizing the violation must not move them
        @test branch_tap_ratio(input, 2) == 1.0

        # the same system with room to spare drives the violation to zero
        feasible = build_test_system(;
            l23_is_transformer = true,
            l23_min_tap_ratio = 0.9,
            l23_max_tap_ratio = 1.1,
        )
        feasible_sol = run_acopf(feasible, ACOPFMinVoltageViolationFormulation())
        @test feasible_sol.termination_status == SOLVED
        @test feasible_sol.objective_value ≈ 0.0 atol = 1e-6
        @test 0.9 - 1e-6 <= feasible_sol.voltage[3] <= 1.1 + 1e-6
    end

    @testset "voltage violation slacks are absent without variable-tap transformers" begin
        # no variable-tap transformer means no penalized bus, so every voltage limit stays a hard
        # bound and the objective is empty
        input = build_test_system()
        sol = run_acopf(input, ACOPFMinVoltageViolationFormulation())

        @test sol.termination_status == SOLVED
        @test sol.objective_value ≈ 0.0 atol = 1e-8
        for i in 1:num_buses(input)
            @test bus_min_voltage(input, i) - 1e-6 <= sol.voltage[i] <= bus_max_voltage(input, i) + 1e-6
        end
    end

    @testset "ACOPF minimizing the reactive balance violation at variable-tap terminals" begin
        # The same system as the voltage violation test above, but now the voltage limits are hard,
        # and bus 3 can only reach its band through a fictitious reactive injection at buses 2 or 3.
        input = build_test_system(;
            l23_is_transformer = true,
            l23_min_tap_ratio = 0.9,
            l23_max_tap_ratio = 1.1,
            pq_min_voltage = 1.15,
            pq_max_voltage = 1.2,
        )
        sol = run_acopf(input, ACOPFMinReactiveBalanceViolationFormulation())

        @test sol.termination_status == SOLVED

        # every voltage stays inside its limits, since they were not relaxed
        for i in 1:num_buses(input)
            @test bus_min_voltage(input, i) - 1e-6 <= sol.voltage[i] <= bus_max_voltage(input, i) + 1e-6
        end

        # the objective is the fictitious reactive injection, summed over the two terminals of L23
        # only, and bus 1 is not a tap terminal, so its balance holds exactly
        @test sol.objective_value > 1e-3
        @test sol.objective_value ≈
              abs(reactive_balance_mismatch(input, sol, 2)) + abs(reactive_balance_mismatch(input, sol, 3)) atol = 1e-6
        @test reactive_balance_mismatch(input, sol, 1) ≈ 0.0 atol = 1e-6
        @test branch_tap_ratio(input, 2) == 1.0

        # the same system with room to spare needs no fictitious injection
        feasible = build_test_system(;
            l23_is_transformer = true,
            l23_min_tap_ratio = 0.9,
            l23_max_tap_ratio = 1.1,
        )
        feasible_sol = run_acopf(feasible, ACOPFMinReactiveBalanceViolationFormulation())
        @test feasible_sol.termination_status == SOLVED
        @test feasible_sol.objective_value ≈ 0.0 atol = 1e-6
        for i in 1:num_buses(feasible)
            @test reactive_balance_mismatch(feasible, feasible_sol, i) ≈ 0.0 atol = 1e-6
        end
    end

    @testset "reactive balance slacks are absent without variable-tap transformers" begin
        input = build_test_system()
        sol = run_acopf(input, ACOPFMinReactiveBalanceViolationFormulation())

        @test sol.termination_status == SOLVED
        @test sol.objective_value ≈ 0.0 atol = 1e-8
        for i in 1:num_buses(input)
            @test reactive_balance_mismatch(input, sol, i) ≈ 0.0 atol = 1e-6
        end
    end

    @testset "build_inputs_from_results loads every setpoint from the solution" begin
        input = build_test_system()
        sol = run_acopf(input)
        rebuilt = build_inputs_from_results(input, sol)

        for i in 1:num_buses(rebuilt)
            @test bus_voltage(rebuilt, i) ≈ sol.voltage[i]
            @test bus_voltage_angle(rebuilt, i) ≈ sol.voltage_angle[i]
        end
        for g in 1:num_generators(rebuilt)
            @test generator_active_power(rebuilt, g) ≈ sol.active_power_generation[g]
            @test generator_reactive_power(rebuilt, g) ≈ sol.reactive_power_generation[g]
            @test generator_voltage_setpoint(rebuilt, g) ≈ sol.voltage[generator_bus_idx(rebuilt, g)]
        end

        # the original is left untouched, and re-solving the rebuilt inputs reproduces the same point
        @test bus_voltage(input, 3) == 1.0
        resolved = run_acpf(rebuilt)
        @test resolved.voltage ≈ sol.voltage atol = 1e-6
        @test resolved.voltage_angle ≈ sol.voltage_angle atol = 1e-6
    end

    @testset "tap sensitivities" begin
        input = build_test_system(; l23_is_transformer = true, l23_tap_ratio = 1.05, l23_min_tap_ratio = 0.9, l23_max_tap_ratio = 1.1)
        sol = run_acopf(input)
        sensitivities = run_tap_sensitivities(input, sol)

        # only L23 has tap slack, and it connects the PV bus 2 to the PQ bus 3
        @test length(sensitivities) == 1
        s = sensitivities[1]
        @test s.branch_idx == 2
        @test s.branch_name == "L23"
        @test (s.from_bus_idx, s.to_bus_idx) == (2, 3)
        @test s.tap_ratio == 1.05

        # the sensitivities are taken at the solution itself, with no power flow re-solve in between
        V, angle = sol.voltage, sol.voltage_angle
        problem = PowerFlow.build_acpf_standard_problem(input)
        G, B = problem.G, problem.B
        @test s.voltage_from ≈ V[2]
        @test s.voltage_to ≈ V[3]
        @test s.angle_difference ≈ angle[2] - angle[3]

        # dQ/dV against a central difference of the bus reactive power injection
        h = 1e-6
        function dQ_dV_finite_difference(i, j)
            V_up, V_down = copy(V), copy(V)
            V_up[j] += h
            V_down[j] -= h
            up = PowerFlow.reactive_power_calc(i, V_up, angle, G, B)
            down = PowerFlow.reactive_power_calc(i, V_down, angle, G, B)
            return (up - down) / (2h)
        end
        @test s.dQ_from_dV_from ≈ dQ_dV_finite_difference(2, 2) atol = 1e-6
        @test s.dQ_from_dV_to ≈ dQ_dV_finite_difference(2, 3) atol = 1e-6
        @test s.dQ_to_dV_from ≈ dQ_dV_finite_difference(3, 2) atol = 1e-6
        @test s.dQ_to_dV_to ≈ dQ_dV_finite_difference(3, 3) atol = 1e-6

        # dQ/dtap against a central difference of the branch reactive injection at each terminal
        function dQ_dtap_finite_difference(side)
            up = PowerFlow.reactive_injection(input, 2, V, angle, side, s.tap_ratio + h)
            down = PowerFlow.reactive_injection(input, 2, V, angle, side, s.tap_ratio - h)
            return (up - down) / (2h)
        end
        @test s.dQ_from_dtap ≈ dQ_dtap_finite_difference(:from) atol = 1e-5
        @test s.dQ_to_dtap ≈ dQ_dtap_finite_difference(:to) atol = 1e-5

        # bus 2 is PV, so the tap moves its reactive generation rather than its voltage; bus 3 is the
        # only free terminal, so the chain rule collapses to a scalar solve
        @test s.dV_from_dtap == 0.0
        @test s.dV_to_dtap ≈ -s.dQ_to_dtap / s.dQ_to_dV_to

        # the sign agrees with re-solving the power flow at a perturbed tap
        perturbed = build_inputs_from_results(input, sol)
        perturbed.static.branches[2].tap_ratio += 0.01
        perturbed_sol = run_acpf(perturbed)
        @test sign(perturbed_sol.voltage[3] - V[3]) == sign(s.dV_to_dtap)
    end

    @testset "dQ/dtap holds for a transformer with charging susceptance" begin
        # the b_c / 2 term of dQ_from/dtap only shows up when the transformer has a shunt, and it
        # only matches the network model because ybus/branch_admittance_terms divide that shunt by
        # tap^2 the way makeYbus.m does
        input = build_test_system(;
            l23_is_transformer = true,
            l23_shunt_susceptance = 0.04,
            l23_tap_ratio = 1.05,
            l23_min_tap_ratio = 0.9,
            l23_max_tap_ratio = 1.1,
        )
        sol = run_acopf(input)
        s = run_tap_sensitivities(input, sol)[1]

        V, angle = sol.voltage, sol.voltage_angle

        h = 1e-6
        function dQ_dtap_finite_difference(side)
            up = PowerFlow.reactive_injection(input, 2, V, angle, side, s.tap_ratio + h)
            down = PowerFlow.reactive_injection(input, 2, V, angle, side, s.tap_ratio - h)
            return (up - down) / (2h)
        end
        @test s.dQ_from_dtap ≈ dQ_dtap_finite_difference(:from) atol = 1e-5
        @test s.dQ_to_dtap ≈ dQ_dtap_finite_difference(:to) atol = 1e-5

        # the shunt contributes to the from side only
        @test s.dQ_from_dtap != 0.0
    end

    @testset "tap sensitivities are empty without tap slack" begin
        input = build_test_system()
        @test isempty(run_tap_sensitivities(input))
    end

end
