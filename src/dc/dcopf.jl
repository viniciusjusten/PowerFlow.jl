function dc_model!(args...)
    angle_variables!(args...)
    generator_variables!(args...)
    power_balance!(args...)
    objective_function!(args...)
    return nothing
end

function run_dcopf(input)
    model = JuMP.Model(HiGHS.Optimizer)
    JuMP.set_silent(model)

    dc_model!(input, model, ModelBuild)
    JuMP.optimize!(model)

    return extract_dc_solution(model)
end
