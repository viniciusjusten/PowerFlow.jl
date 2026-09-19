"""
    build_inputs_from_results(input, solution::ACSolution; tap_ratio = nothing)

Build a copy of `input` with every setpoint replaced by the corresponding value in `solution`:
bus voltage magnitude and angle, and generator active/reactive power. Each generator's voltage
setpoint is also moved to the solved voltage of its bus, so that re-solving the copy with a
formulation that regulates PV/Slack voltages (the standard power flow in `acpf_standard.jl`,
`ACPFFeasibilityFormulation`) reproduces the same operating point instead of pulling the system
back to the original data.

`solution` carries no tap ratios, so a solution that came out of `ACOPFTapFormulation` must pass
its optimized taps through `tap_ratio` (one entry per branch) for the copy to represent that
operating point.
"""
function build_inputs_from_results(input, solution::ACSolution; tap_ratio = nothing)
    new_input = deepcopy(input)

    for i in 1:num_buses(new_input)
        bus = new_input.static.buses[i]
        bus.voltage = solution.voltage[i]
        bus.voltage_angle = solution.voltage_angle[i]
    end

    for g in 1:num_generators(new_input)
        generator = new_input.static.generators[g]
        generator.active_power = solution.active_power_generation[g]
        generator.reactive_power = solution.reactive_power_generation[g]
        generator.voltage_setpoint = solution.voltage[generator.bus_idx]
    end

    if tap_ratio !== nothing
        for k in 1:num_branches(new_input)
            new_input.static.branches[k].tap_ratio = tap_ratio[k]
        end
    end

    return new_input
end
