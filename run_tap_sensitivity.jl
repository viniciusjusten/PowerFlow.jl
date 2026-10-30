using PowerModels, PowerFlow, PowerFlow.PowerSystemData
using JuMP, Ipopt, HiGHS

path = raw"D:\codes\PowerFlowExperiments\pglib_benchmark\pglib-opf\pglib_opf_case30_ieee.m"
PowerModels.silence()
pm_data = PowerModels.parse_file(path)

inputs = PowerSystemData.from_powermodels(pm_data)

transformers = [
    11,
    12,
    13,
    14,
    15,
    16,
    36,
]
feasible = [
    0.978,
    0.969,
    1.0,
    1.0,
    0.932,
    1.0,
    0.968,
]

for (i, k) in enumerate(transformers)
    inputs.static.branches[k].min_tap_ratio = 0.9
    inputs.static.branches[k].max_tap_ratio = 1.1
    inputs.static.branches[k].tap_ratio = 1.1
end

# Run DCOPF
model = JuMP.Model(HiGHS.Optimizer)
PowerFlow.dc_model!(inputs, model, PowerFlow.ModelBuild)
JuMP.optimize!(model)
dc_sol = PowerFlow.extract_dc_solution(model)

# Update generation bounds based on DCOPF solution
for (i, gen) in enumerate(inputs.static.generators)
    bus_idx = gen.bus_idx
    if inputs.static.buses[bus_idx].type == PowerFlow.BusType.Slack
        continue
    end
    Pg = dc_sol.active_power_generation[i]
    gen.min_active_power = max(0.8 * Pg, gen.min_active_power)
    gen.max_active_power = min(1.2 * Pg, gen.max_active_power)
end

# Run ACOPF 
formulation = PowerFlow.ACOPFMinReactiveBalanceViolationFormulation()
model = JuMP.Model(Ipopt.Optimizer)
set_attribute(model, "max_iter", 300)
PowerFlow.ac_model!(inputs, model, formulation, PowerFlow.ModelBuild)
JuMP.optimize!(model)
ac_feasibility_sol = PowerFlow.extract_ac_solution(model)
