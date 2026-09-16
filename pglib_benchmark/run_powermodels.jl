import Pkg
Pkg.activate(joinpath(@__DIR__, "pm_env"))

using PowerModels
using JuMP
using Ipopt
using HiGHS
using JSON

PowerModels.silence()

const CASE_DIR = joinpath(@__DIR__, "pglib-opf")
const CACHE_DIR = joinpath(@__DIR__, "cache")
const MAX_BUSES = if length(ARGS) >= 1
    parse(Int, ARGS[1])
else
    1803
end
const DC_CSV = joinpath(@__DIR__, "pm_dcopf_partial.csv")
const AC_CSV = joinpath(@__DIR__, "pm_acopf_partial.csv")

mkpath(CACHE_DIR)

# Only the "typical operating condition" cases (the .m files at the root of
# pglib-opf); the api/ and sad/ variants are not included.
function list_cases()
    cases = NamedTuple{(:name, :path, :buses),Tuple{String,String,Int}}[]
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
        push!(cases, (name = f[1:end-2], path = joinpath(CASE_DIR, f), buses = buses))
    end
    sort!(cases; by = c -> c.buses)
    return cases
end

# PowerFlow.jl / PowerSystemData.from_powermodels only keeps the linear term
# of a polynomial generator cost curve (ignores the quadratic and constant
# terms). To compare PowerModels' own OPF against PowerFlow's on the same
# cost function, zero out everything but that linear term here too.
function linearize_cost!(data)
    for gen in values(data["gen"])
        coeffs = get(gen, "cost", Float64[])
        model = Int(get(gen, "model", 2))
        if model != 2 || length(coeffs) < 2
            gen["cost"] = zeros(Float64, max(length(coeffs), 1))
            continue
        end
        linear = coeffs[end-1]
        new_coeffs = zeros(Float64, length(coeffs))
        new_coeffs[end-1] = linear
        gen["cost"] = new_coeffs
    end
    return data
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

    dc_optimizer = JuMP.optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false)
    ac_optimizer = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)

    total = length(cases)
    for (i, case) in enumerate(cases)
        need_dc = !(case.name in dc_done)
        need_ac = !(case.name in ac_done)
        cache_path = joinpath(CACHE_DIR, case.name * ".json")
        need_cache = !isfile(cache_path)

        if !need_dc && !need_ac && !need_cache
            println("[$i/$total] $(case.name): already done, skipping")
            continue
        end

        println("[$i/$total] $(case.name) ($(case.buses) buses)")

        local pm_data
        try
            pm_data = PowerModels.parse_matpower(case.path)
        catch e
            println("  parse error: ", sprint(showerror, e))
            continue
        end

        if need_cache
            open(cache_path, "w") do io
                JSON.print(io, pm_data)
            end
        end

        if need_dc
            try
                data = deepcopy(pm_data)
                linearize_cost!(data)
                t = @elapsed result = PowerModels.solve_dc_opf(data, dc_optimizer)
                println("  dc: cost=$(result["objective"]) time=$(round(t; digits = 3))s status=$(result["termination_status"])")
                append_row!(DC_CSV, "$(case.name),$(case.buses),$(result["objective"]),$t,$(result["termination_status"])")
            catch e
                println("  dc error: ", sprint(showerror, e))
                append_row!(DC_CSV, "$(case.name),$(case.buses),NaN,NaN,error")
            end
        end

        if need_ac
            try
                data = deepcopy(pm_data)
                linearize_cost!(data)
                t = @elapsed result = PowerModels.solve_ac_opf(data, ac_optimizer)
                println("  ac: cost=$(result["objective"]) time=$(round(t; digits = 3))s status=$(result["termination_status"])")
                append_row!(AC_CSV, "$(case.name),$(case.buses),$(result["objective"]),$t,$(result["termination_status"])")
            catch e
                println("  ac error: ", sprint(showerror, e))
                append_row!(AC_CSV, "$(case.name),$(case.buses),NaN,NaN,error")
            end
        end
    end
end

main()
