struct ACSolution
    voltage::Vector{Float64}
    voltage_angle::Vector{Float64}
    active_power_generation::Vector{Float64}
    reactive_power_generation::Vector{Float64}
    objective_value::Float64
    termination_status::Any
    # TODO - add tap solution and duals
end

struct DCSolution
    voltage_angle::Vector{Float64}
    active_power_generation::Vector{Float64}
    objective_value::Float64
    termination_status::Any
    # TODO - add duals
end

function extract_ac_solution(model::JuMP.Model)
    return ACSolution(
        JuMP.value.(model[:V]),
        JuMP.value.(model[:angle]),
        JuMP.value.(model[:Pg]),
        JuMP.value.(model[:Qg]),
        JuMP.objective_value(model),
        JuMP.termination_status(model),
    )
end

function extract_dc_solution(model::JuMP.Model)
    return DCSolution(
        JuMP.value.(model[:angle]),
        JuMP.value.(model[:Pg]),
        JuMP.objective_value(model),
        JuMP.termination_status(model),
    )
end
