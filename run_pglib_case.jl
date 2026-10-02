using PowerModels, PowerFlow, PowerFlow.PowerSystemData

# path = raw"D:\codes\PowerFlowExperiments\pglib_benchmark\pglib-opf\pglib_opf_case30_ieee.m"
path = raw"D:\codes\PowerFlow.jl\case.m"
# path = "case30_dcline.m"
PowerModels.silence()
pm_data = PowerModels.parse_file(path)

inputs = PowerSystemData.from_powermodels(pm_data)
sol = PowerFlow.run_acopf(inputs, PowerFlow.ACOPFFormulation())
