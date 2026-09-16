using PowerModels, PowerFlow, PowerFlow.PowerSystemData

path = "pglib_benchmark/pglib-opf/pglib_opf_case30_ieee.m"
PowerModels.silence()
pm_data = PowerModels.parse_file(path)

inputs = PowerSystemData.from_powermodels(pm_data)
sol = PowerFlow.run_acopf(inputs, PowerFlow.ACOPFFormulation())
