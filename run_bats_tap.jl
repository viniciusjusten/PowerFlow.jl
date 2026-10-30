using PowerModels
import PowerFlow as PF
import PowerSystemData as PSD

function read_tap_limits(path)
    lines = filter(!isempty, strip.(readlines(path)))
    header = strip.(split(lines[1], ','))
    from_col = findfirst(==("barra_from"), header)
    to_col = findfirst(==("barra_to"), header)
    min_col = findfirst(==("tmn"), header)
    max_col = findfirst(==("tmx"), header)
    if any(isnothing, (from_col, to_col, min_col, max_col))
        error("$path must have the columns barra_from, barra_to, tmn and tmx")
    end

    limits = []
    for line in lines[2:end]
        fields = strip.(split(line, ','))
        push!(limits, (
            from_bus_name = String(fields[from_col]),
            to_bus_name = String(fields[to_col]),
            min_tap_ratio = parse(Float64, fields[min_col]),
            max_tap_ratio = parse(Float64, fields[max_col]),
        ))
    end
    return limits
end

function transformers_by_bus_names(inputs)
    transformers = Dict{Tuple{String, String}, Vector{Int}}()
    for (k, branch) in enumerate(inputs.static.branches)
        if !branch.is_transformer
            continue
        end
        key = (
            inputs.static.buses[branch.from_bus_idx].name,
            inputs.static.buses[branch.to_bus_idx].name,
        )
        push!(get!(transformers, key, Int[]), k)
    end
    return transformers
end

function update_tap_limits_from_csv!(inputs, path)
    transformers = transformers_by_bus_names(inputs)
    rows_used = Dict{Tuple{String, String}, Int}()
    updated = Int[]
    unmatched = []
    for row in read_tap_limits(path)
        key = (row.from_bus_name, row.to_bus_name)
        indices = get(transformers, key, Int[])
        n = get(rows_used, key, 0) + 1
        if n > length(indices)
            push!(unmatched, row)
            continue
        end
        rows_used[key] = n
        k = indices[n]
        inputs.static.branches[k].min_tap_ratio = row.min_tap_ratio
        inputs.static.branches[k].max_tap_ratio = row.max_tap_ratio
        push!(updated, k)
    end
    return updated, unmatched
end

case_path = joinpath(@__DIR__, "case_hour1_scenario1_sem_shunts_toml.m")
tap_csv_path = joinpath(@__DIR__, "var_tap.csv")

PowerModels.silence()
pm_data = PowerModels.parse_file(case_path)

# keep only the largest connected component, where DC lines also count as connections, and then
# switch off every branch, generator, load and DC line attached to a bus that became inactive
PowerModels.select_largest_component!(pm_data)
PowerModels.propagate_topology_status!(pm_data)

inputs = PSD.from_powermodels(pm_data)

updated, unmatched = update_tap_limits_from_csv!(inputs, tap_csv_path)
println("Updated the tap limits of $(length(updated)) transformers")
println("$(length(unmatched)) rows of $(basename(tap_csv_path)) matched no transformer")
