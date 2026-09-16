# Joins the PowerFlow.jl and PowerModels partial result CSVs (written
# incrementally by run_powerflow.jl / run_powermodels.jl) into the final
# dcopf.csv / acopf.csv comparison tables. Uses only Base, no environment
# needed. Safe to re-run at any point to refresh the comparison with
# whatever partial results exist so far.

function read_csv(path)
    rows = Dict{String,Vector{SubString{String}}}()
    if !isfile(path)
        return rows
    end
    open(path) do io
        readline(io) # header
        for line in eachline(io)
            if isempty(line)
                continue
            end
            fields = split(line, ',')
            rows[fields[1]] = fields
        end
    end
    return rows
end

function case_bus_count(case_name)
    m = match(r"case(\d+)_", case_name)
    if m === nothing
        return typemax(Int)
    else
        return parse(Int, m.captures[1])
    end
end

function merge_side(pf_path, pm_path, out_path)
    pf_rows = read_csv(pf_path)
    pm_rows = read_csv(pm_path)
    cases = sort(collect(union(keys(pf_rows), keys(pm_rows))); by = case_bus_count)

    open(out_path, "w") do io
        println(
            io,
            "case,buses,powerflow_cost,powerflow_time_s,powerflow_status,powermodels_cost,powermodels_time_s,powermodels_status",
        )
        for case in cases
            pf = get(pf_rows, case, nothing)
            pm = get(pm_rows, case, nothing)

            if pf !== nothing
                buses = pf[2]
            elseif pm !== nothing
                buses = pm[2]
            else
                buses = ""
            end

            if pf !== nothing
                pf_cost, pf_time, pf_status = pf[3], pf[4], pf[5]
            else
                pf_cost, pf_time, pf_status = "", "", ""
            end

            if pm !== nothing
                pm_cost, pm_time, pm_status = pm[3], pm[4], pm[5]
            else
                pm_cost, pm_time, pm_status = "", "", ""
            end

            println(io, "$case,$buses,$pf_cost,$pf_time,$pf_status,$pm_cost,$pm_time,$pm_status")
        end
    end
    println("Wrote $out_path ($(length(cases)) cases)")
    return nothing
end

dir_path = @__DIR__

merge_side(joinpath(dir_path, "pf_dcopf_partial.csv"), joinpath(dir_path, "pm_dcopf_partial.csv"), joinpath(dir_path, "dcopf.csv"))
merge_side(joinpath(dir_path, "pf_acopf_partial.csv"), joinpath(dir_path, "pm_acopf_partial.csv"), joinpath(dir_path, "acopf.csv"))
