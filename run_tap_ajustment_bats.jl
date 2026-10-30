using PowerModels
import PowerFlow as PF
import PowerSystemData as PSD
using JuMP, Ipopt, HiGHS, HSL_jll

function solve_min_reactive_violation(inputs)
    formulation = PF.ACOPFMinReactiveBalanceViolationFormulation()
    opt = JuMP.optimizer_with_attributes(
        Ipopt.Optimizer,
        "max_iter" => 300,
        "hsllib" => HSL_jll.libhsl_path,
        "linear_solver" => "ma57",
    )
    model = JuMP.Model(opt)
    JuMP.set_silent(model)
    set_attribute(model, "max_iter", 300)
    PF.ac_model!(inputs, model, formulation, PF.ModelBuild)
    JuMP.optimize!(model)

    solution = PF.extract_ac_solution(model)
    Q_deficit = Dict{Int, Float64}()
    Q_surplus = Dict{Int, Float64}()
    for i in PF.variable_tap_buses(inputs)
        Q_deficit[i] = JuMP.value(model[:Q_deficit][i])
        Q_surplus[i] = JuMP.value(model[:Q_surplus][i])
    end

    return solution, JuMP.termination_status(model), Q_deficit, Q_surplus
end

function terminal_reactive_sensitivity(s::PF.TapSensitivity, bus)
    if bus == s.to_bus_idx
        return s.dQ_to_dtap
    end
    return s.dQ_from_dtap
end

function newton_tap_steps(sensitivities, Q_deficit, Q_surplus, tolerance)
    steps = Dict{Int, Float64}()
    for bus in sort!(collect(keys(Q_deficit)))
        violation = Q_deficit[bus] - Q_surplus[bus]
        if abs(violation) < tolerance
            continue
        end

        chosen = nothing
        chosen_sensitivity = 0.0
        for s in sensitivities
            if bus != s.from_bus_idx && bus != s.to_bus_idx
                continue
            end
            sensitivity = terminal_reactive_sensitivity(s, bus)
            if isnothing(chosen) || abs(sensitivity) > abs(chosen_sensitivity)
                chosen = s
                chosen_sensitivity = sensitivity
            end
        end

        if isnothing(chosen) || chosen_sensitivity == 0.0
            continue
        end
        steps[chosen.branch_idx] = get(steps, chosen.branch_idx, 0.0) - violation / chosen_sensitivity
    end

    return steps
end

function variable_taps(inputs)
    taps = Dict{Int, Float64}()
    for k in 1:PSD.num_branches(inputs)
        if PF.branch_has_variable_tap(inputs, k)
            taps[k] = PSD.branch_tap_ratio(inputs, k)
        end
    end
    return taps
end

function restore_taps!(inputs, taps)
    for (k, tap) in taps
        inputs.static.branches[k].tap_ratio = tap
    end
    return nothing
end

function set_taps!(inputs, base_taps, steps, step_size)
    largest_change = 0.0
    for (k, base_tap) in base_taps
        branch = inputs.static.branches[k]
        new_tap = clamp(base_tap + step_size * get(steps, k, 0.0), branch.min_tap_ratio, branch.max_tap_ratio)
        largest_change = max(largest_change, abs(new_tap - base_tap))
        branch.tap_ratio = new_tap
    end
    return largest_change
end

function write_csv_headers(summary_path, violations_path, taps_path)
    open(summary_path, "w") do io
        println(io, "iteration,step_size,status,objective,largest_Q_deficit,largest_Q_surplus,accepted")
    end
    open(violations_path, "w") do io
        println(io, "iteration,bus_idx,bus_name,Q_deficit,Q_surplus")
    end
    open(taps_path, "w") do io
        println(io, "iteration,branch_idx,branch_name,from_bus_idx,to_bus_idx,tap_ratio")
    end
    return nothing
end

function append_iteration(violations_path, taps_path, iteration, inputs, Q_deficit, Q_surplus)
    open(violations_path, "a") do io
        for i in sort!(collect(keys(Q_deficit)))
            println(io, "$iteration,$i,$(PSD.bus_name(inputs, i)),$(Q_deficit[i]),$(Q_surplus[i])")
        end
    end
    open(taps_path, "a") do io
        for k in 1:PSD.num_branches(inputs)
            if !PowerFlow.branch_has_variable_tap(inputs, k)
                continue
            end
            println(
                io,
                "$iteration,$k,$(PSD.branch_name(inputs, k)),$(PSD.branch_from_bus_idx(inputs, k))," *
                "$(PSD.branch_to_bus_idx(inputs, k)),$(PSD.branch_tap_ratio(inputs, k))",
            )
        end
    end
    return nothing
end

function append_summary(summary_path, iteration, step_size, status, objective, Q_deficit, Q_surplus, accepted)
    largest_deficit = maximum(values(Q_deficit); init = 0.0)
    largest_surplus = maximum(values(Q_surplus); init = 0.0)
    open(summary_path, "a") do io
        println(io, "$iteration,$step_size,$status,$objective,$largest_deficit,$largest_surplus,$accepted")
    end
    println(
        "iteration $iteration: step size $step_size, status $status, objective $objective, " *
        "largest Q_deficit $largest_deficit, largest Q_surplus $largest_surplus, accepted $accepted",
    )
    return nothing
end

function is_solved(status)
    return status in (JuMP.MOI.LOCALLY_SOLVED, JuMP.MOI.ALMOST_LOCALLY_SOLVED)
end

function run_tap_adjustment!(
    inputs;
    tolerance = 1e-6,
    max_iterations = 50,
    min_step_size = 1e-3,
    stall_tolerance = 1e-3,
    stall_iterations = 3,
    output_dir = @__DIR__,
)
    summary_path = joinpath(output_dir, "tap_adjustment_summary.csv")
    violations_path = joinpath(output_dir, "tap_adjustment_violations.csv")
    taps_path = joinpath(output_dir, "tap_adjustment_taps.csv")
    write_csv_headers(summary_path, violations_path, taps_path)

    accepted_taps = variable_taps(inputs)
    accepted_solution, status, Q_deficit, Q_surplus = solve_min_reactive_violation(inputs)
    append_iteration(violations_path, taps_path, 0, inputs, Q_deficit, Q_surplus)
    append_summary(summary_path, 0, 0.0, status, accepted_solution.objective_value, Q_deficit, Q_surplus, true)
    if !is_solved(status)
        println("Stopping because the ACOPF did not solve at the initial taps.")
        return accepted_solution
    end

    # the Newton steps only change when a trial is accepted, so they are kept between rejections
    sensitivities = PF.build_tap_sensitivities(inputs, accepted_solution)
    PF.reactive_power_tap_derivatives!(sensitivities, inputs)
    steps = newton_tap_steps(sensitivities, Q_deficit, Q_surplus, tolerance)

    step_size = 1.0
    stalled_count = 0
    for iteration in 1:max_iterations
        if maximum(values(Q_deficit); init = 0.0) < tolerance && maximum(values(Q_surplus); init = 0.0) < tolerance
            println("Converged: every reactive balance violation is below $tolerance.")
            break
        end

        if set_taps!(inputs, accepted_taps, steps, step_size) == 0.0
            println("Stopping because every tap that should move is already at its limit.")
            break
        end

        trial_solution, trial_status, trial_deficit, trial_surplus = solve_min_reactive_violation(inputs)
        accepted = is_solved(trial_status) && trial_solution.objective_value < accepted_solution.objective_value
        append_iteration(violations_path, taps_path, iteration, inputs, trial_deficit, trial_surplus)
        append_summary(
            summary_path,
            iteration,
            step_size,
            trial_status,
            trial_solution.objective_value,
            trial_deficit,
            trial_surplus,
            accepted,
        )

        if accepted
            improvement =
                (accepted_solution.objective_value - trial_solution.objective_value) / accepted_solution.objective_value
            accepted_taps = variable_taps(inputs)
            accepted_solution, Q_deficit, Q_surplus = trial_solution, trial_deficit, trial_surplus
            # step_size = min(1.0, 2 * step_size)
            efficiency = (accepted_solution.objective_value - trial_solution.objective_value) /
                (step_size * accepted_solution.objective_value)
            step_size = clamp(0.8 / efficiency, 0.25, 5.0)

            if improvement < stall_tolerance
                stalled_count += 1
            else
                stalled_count = 0
            end
            if stalled_count >= stall_iterations
                println(
                    "Stopping because the solution stagnated: $stall_iterations accepted trials in a row " *
                    "improved the objective by less than $(100 * stall_tolerance)%.",
                )
                break
            end

            sensitivities = PF.build_tap_sensitivities(inputs, accepted_solution)
            PF.reactive_power_tap_derivatives!(sensitivities, inputs)
            steps = newton_tap_steps(sensitivities, Q_deficit, Q_surplus, tolerance)
        else
            step_size /= 2
            if step_size < min_step_size
                println("Stopping because the solution stagnated: the step size fell below $min_step_size.")
                break
            end
        end

        if iteration == max_iterations
            println("Stopping after $max_iterations trials without converging.")
        end
    end

    restore_taps!(inputs, accepted_taps)
    return accepted_solution
end

include("run_bats_tap.jl")

# Run DCOPF
model = JuMP.Model(HiGHS.Optimizer)
PF.dc_model!(inputs, model, PF.ModelBuild)
JuMP.optimize!(model)
dc_sol = PF.extract_dc_solution(model)

# Update generation bounds based on DCOPF solution
for (i, gen) in enumerate(inputs.static.generators)
    bus_idx = gen.bus_idx
    if inputs.static.buses[bus_idx].type == PSD.BusType.Slack
        continue
    end
    Pg = dc_sol.active_power_generation[i]
    gen.min_active_power = max(0.8 * Pg, gen.min_active_power)
    gen.max_active_power = min(1.2 * Pg, gen.max_active_power)
end

# Adjust the taps iteratively with the ACOPF that minimizes the reactive balance violation
ac_feasibility_sol = run_tap_adjustment!(inputs; tolerance = 1e-6)
