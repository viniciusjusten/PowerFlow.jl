function objective_function!(input, model, ::Type{ModelBuild})
    Pg = model[:Pg]
    @objective(model, Min, sum(PSD.generator_cost(input, g) * Pg[g] for g in 1:PSD.num_generators(input)))
    return nothing
end
