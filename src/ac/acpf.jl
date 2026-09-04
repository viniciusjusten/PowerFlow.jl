function run_acpf(input, formulation::ACFormulation)
    model = JuMP.Model(Ipopt.Optimizer)
    JuMP.set_silent(model)

    ac_model!(input, model, formulation, ModelBuild)
    JuMP.optimize!(model)

    return extract_ac_solution(model)
end

run_acpf(input) = run_acpf(input, ACPFFormulation())
