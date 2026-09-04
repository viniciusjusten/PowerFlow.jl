function ac_model!(args...)
    voltage_variables!(args...)
    generator_variables!(args...)
    power_balance!(args...)
    voltage_setpoint!(args...)
    branch_flow_limit!(args...)
    objective_function!(args...)
    return nothing
end

function run_acopf(input, formulation::ACFormulation = ACOPFFormulation())
    model = JuMP.Model(Ipopt.Optimizer)
    JuMP.set_silent(model)

    ac_model!(input, model, formulation, ModelBuild)
    JuMP.optimize!(model)

    return extract_ac_solution(model)
end
