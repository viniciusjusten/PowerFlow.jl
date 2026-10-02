using PowerModels, PowerFlow, PowerFlow.PowerSystemData

# path = raw"D:\codes\PowerFlowExperiments\pglib_benchmark\pglib-opf\pglib_opf_case30_ieee.m"
path = raw"D:\codes\PowerFlow.jl\case_hour1_scenario1.m"
# path = "case30_dcline.m"
PowerModels.silence()
pm_data = PowerModels.parse_file(path)

# keep only the largest connected component, where DC lines also count as connections, and then
# switch off every branch, generator, load and DC line attached to a bus that became inactive
PowerModels.select_largest_component!(pm_data)
PowerModels.propagate_topology_status!(pm_data)

inputs = PowerSystemData.from_powermodels(pm_data)
sol = PowerFlow.run_acopf(inputs, PowerFlow.ACOPFFormulation())
