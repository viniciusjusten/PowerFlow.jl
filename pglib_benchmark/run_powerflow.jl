import Pkg
Pkg.activate(@__DIR__)

using PowerSystemData
using PowerFlow
using JSON

const CASE_DIR = joinpath(@__DIR__, "pglib-opf")
const CACHE_DIR = joinpath(@__DIR__, "cache")
const MAX_BUSES = if length(ARGS) >= 1
    parse(Int, ARGS[1])
else
    1803
end
const DC_CSV = joinpath(@__DIR__, "pf_dcopf_partial.csv")
const AC_CSV = joinpath(@__DIR__, "pf_acopf_partial.csv")

# This environment does not depend on PowerModels (adding it here would
# conflict with PowerFlow's NonlinearSolve dependency: PowerModels pins
# NLsolve to 4.x while NonlinearSolve's NLsolve extension needs 5.x). So the
# per-unit PowerModels dict for each case is read from the JSON cache that
# run_powermodels.jl writes, instead of being parsed here.
function list_cases()
    cases = NamedTuple{(:name, :buses),Tuple{String,Int}}[]
    for f in readdir(CASE_DIR)
        if !endswith(f, ".m")
            continue
        end
        m = match(r"^pglib_opf_case(\d+)_", f)
        if m === nothing
            continue
        end
        buses = parse(Int, m.captures[1])
        if buses > MAX_BUSES
            continue
        end
        push!(cases, (name = f[1:end-2], buses = buses))
    end
    sort!(cases; by = c -> c.buses)
    return cases
end

function read_done_cases(csv_path)
    done = Set{String}()
    if !isfile(csv_path)
        return done
    end
    for (i, line) in enumerate(eachline(csv_path))
        if i == 1
            continue
        end
        if isempty(line)
            continue
        end
        push!(done, split(line, ',')[1])
    end
    return done
end

function ensure_header!(csv_path, header)
    if isfile(csv_path)
        return nothing
    end
    open(csv_path, "w") do io
        println(io, header)
    end
    return nothing
end

function append_row!(csv_path, row)
    open(csv_path, "a") do io
        println(io, row)
    end
    return nothing
end

function main()
    cases = list_cases()
    dc_done = read_done_cases(DC_CSV)
    ac_done = read_done_cases(AC_CSV)
    ensure_header!(DC_CSV, "case,buses,cost,time_s,status")
    ensure_header!(AC_CSV, "case,buses,cost,time_s,status")

    total = length(cases)
    for (i, case) in enumerate(cases)
        need_dc = !(case.name in dc_done)
        need_ac = !(case.name in ac_done)
        if !need_dc && !need_ac
            println("[$i/$total] $(case.name): already done, skipping")
            continue
        end

        cache_path = joinpath(CACHE_DIR, case.name * ".json")
        if !isfile(cache_path)
            println("[$i/$total] $(case.name): no cache found, run run_powermodels.jl first — skipping")
            continue
        end

        println("[$i/$total] $(case.name) ($(case.buses) buses)")

        local input
        try
            pm_data = JSON.parsefile(cache_path)
            input = PowerSystemData.from_powermodels(pm_data)
        catch e
            println("  from_powermodels error: ", sprint(showerror, e))
            continue
        end

        if need_dc
            try
                t = @elapsed sol = run_dcopf(input)
                println("  dc: cost=$(sol.objective_value) time=$(round(t; digits = 3))s status=$(sol.termination_status)")
                append_row!(DC_CSV, "$(case.name),$(case.buses),$(sol.objective_value),$t,$(sol.termination_status)")
            catch e
                println("  dc error: ", sprint(showerror, e))
                append_row!(DC_CSV, "$(case.name),$(case.buses),NaN,NaN,error")
            end
        end

        if need_ac
            try
                t = @elapsed sol = run_acopf(input)
                println("  ac: cost=$(sol.objective_value) time=$(round(t; digits = 3))s status=$(sol.termination_status)")
                append_row!(AC_CSV, "$(case.name),$(case.buses),$(sol.objective_value),$t,$(sol.termination_status)")
            catch e
                println("  ac error: ", sprint(showerror, e))
                append_row!(AC_CSV, "$(case.name),$(case.buses),NaN,NaN,error")
            end
        end
    end
end

main()
