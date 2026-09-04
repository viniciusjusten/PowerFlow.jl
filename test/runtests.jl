using PowerSystemData
using PowerFlow
using JuMP
using Ipopt
using Test

function build_test_system(; l13_max_flow::Float64 = 5.0)
    input = PowerSystemData.Input()

    add_bus!(input; name = "Slack", type = BusType.Slack, voltage = 1.0, min_voltage = 0.9, max_voltage = 1.1)
    add_bus!(input; name = "PV", type = BusType.PV, voltage = 1.0, min_voltage = 0.9, max_voltage = 1.1)
    add_bus!(input; name = "PQ", type = BusType.PQ, min_voltage = 0.9, max_voltage = 1.1)

    add_branch!(
        input;
        name = "L12",
        from_bus_idx = 1,
        to_bus_idx = 2,
        resistance = 0.01,
        reactance = 0.1,
        max_flow = 5.0,
    )
    add_branch!(
        input;
        name = "L23",
        from_bus_idx = 2,
        to_bus_idx = 3,
        resistance = 0.01,
        reactance = 0.1,
        max_flow = 5.0,
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

    @testset "ACOPF respects branch thermal limits" begin
        input = build_test_system(; l13_max_flow = 0.3)
        sol = run_acopf(input)

        model = JuMP.Model(Ipopt.Optimizer)
        JuMP.set_silent(model)
        PowerFlow.ac_model!(input, model, ACOPFFormulation(), PowerFlow.ModelBuild)
        JuMP.optimize!(model)

        Pkm, Qkm = model[:Pkm], model[:Qkm]
        S_from = sqrt(JuMP.value(Pkm[3, :from])^2 + JuMP.value(Qkm[3, :from])^2)
        S_to = sqrt(JuMP.value(Pkm[3, :to])^2 + JuMP.value(Qkm[3, :to])^2)

        @test sol.termination_status == SOLVED
        @test S_from <= 0.3 + 1e-4
        @test S_to <= 0.3 + 1e-4
        # Pkm/Qkm are registered for every branch, not just rated ones, so they can be reused elsewhere
        @test size(Pkm) == (3, 2)
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

        P_from = PowerFlow.branch_active_power_flow(input, 3, sol.voltage, sol.voltage_angle, :from)
        Q_from = PowerFlow.branch_reactive_power_flow(input, 3, sol.voltage, sol.voltage_angle, :from)
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
end
