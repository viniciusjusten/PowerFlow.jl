using PowerModels, PowerFlow, PowerFlow.PowerSystemData
using JuMP, Ipopt, Printf

const CASE_PATH = joinpath(@__DIR__, "case30_dcline.m")

"""
    solve_powermodels(pm_data)

Solve the PowerModels ACOPF (polar formulation) for the data with Ipopt and return the result
dictionary.
"""
function solve_powermodels(pm_data)
    solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
    return PowerModels.run_ac_opf(pm_data, solver)
end

"""
    sorted_values(component_dict, key)

Return the field `key` of every entry of a PowerModels component dictionary, ordered by the entry
index. This is the same order in which PowerSystemData creates its components, as long as every
component is in service.
"""
function sorted_values(component_dict, key)
    ids = sort(parse.(Int, collect(keys(component_dict))))
    return [component_dict[string(i)][key] for i in ids]
end

"""
    compare(pm_data, label)

Solve the case with both packages and print the termination status, the objective, the largest
difference in each bus, generator and branch quantity, and every DC line quantity side by side.
All values are in per unit and every flow is the power leaving the AC bus at that terminal.
"""
function compare(pm_data, label)
    input = PowerSystemData.from_powermodels(pm_data)
    pf_sol = PowerFlow.run_acopf(input, PowerFlow.ACOPFFormulation())
    pm_result = solve_powermodels(pm_data)
    pm_sol = pm_result["solution"]

    println("="^72)
    println(label)
    println("="^72)
    @printf("status      PF = %-20s PM = %s\n", pf_sol.termination_status, pm_result["termination_status"])
    pm_obj = pm_result["objective"]
    @printf(
        "objective   PF = %.6f   PM = %.6f   rel diff = %.2e\n",
        pf_sol.objective_value,
        pm_obj,
        abs(pf_sol.objective_value - pm_obj) / abs(pm_obj),
    )

    pairs = [
        ("bus vm", pf_sol.voltage, sorted_values(pm_sol["bus"], "vm")),
        ("bus va", pf_sol.voltage_angle, sorted_values(pm_sol["bus"], "va")),
        ("gen pg", pf_sol.active_power_generation, sorted_values(pm_sol["gen"], "pg")),
        ("gen qg", pf_sol.reactive_power_generation, sorted_values(pm_sol["gen"], "qg")),
        ("branch pf", pf_sol.branch_active_power_from, sorted_values(pm_sol["branch"], "pf")),
        ("branch pt", pf_sol.branch_active_power_to, sorted_values(pm_sol["branch"], "pt")),
        ("branch qf", pf_sol.branch_reactive_power_from, sorted_values(pm_sol["branch"], "qf")),
        ("branch qt", pf_sol.branch_reactive_power_to, sorted_values(pm_sol["branch"], "qt")),
    ]
    for (name, pf_values, pm_values) in pairs
        @printf("max |Δ %-9s| = %.2e\n", name, maximum(abs.(pf_values .- pm_values)))
    end

    println()
    println("dcline     PF            PM            Δ")
    dcline_pairs = [
        ("pf", pf_sol.dcline_active_power_from, "pf"),
        ("pt", pf_sol.dcline_active_power_to, "pt"),
        ("qf", pf_sol.dcline_reactive_power_from, "qf"),
        ("qt", pf_sol.dcline_reactive_power_to, "qt"),
    ]
    for (name, pf_values, pm_key) in dcline_pairs
        pm_values = sorted_values(pm_sol["dcline"], pm_key)
        for k in eachindex(pf_values)
            @printf(
                "%s[%d]     % .6f     % .6f     % .2e\n",
                name,
                k,
                pf_values[k],
                pm_values[k],
                pf_values[k] - pm_values[k]
            )
        end
    end
    println()

    return nothing
end

PowerModels.silence()
pm_data = PowerModels.parse_file(CASE_PATH)

compare(pm_data, "case30 with one DC line (loose reactive limits)")

# A one-sided reactive window at the to terminal. With loose limits a sign error in the reactive
# bounds of the DC line goes unnoticed, but here the two conventions put qt on opposite sides of
# zero, so the comparison shows which one each package follows.
tight_data = deepcopy(pm_data)
tight_data["dcline"]["1"]["qmint"] = 0.02
tight_data["dcline"]["1"]["qmaxt"] = 0.05
compare(tight_data, "case30 with one DC line (PowerModels qmint, qmaxt = 0.02, 0.05 p.u.)")
