abstract type ACFormulation end

struct ACOPFFormulation <: ACFormulation end
struct ACOPFMinRedispatchFormulation <: ACFormulation end
struct ACOPFTapFormulation <: ACFormulation end
struct ACOPFMinVoltageViolationFormulation <: ACFormulation end
struct ACPFFeasibilityFormulation <: ACFormulation end
struct ACPFPmisVmisFormulation <: ACFormulation end

# standalone marker for the NonlinearSolve.jl path — not an ACFormulation, it never builds a JuMP model
struct ACPFFormulation end

allows_generation_redispatch(::ACFormulation) = false
allows_generation_redispatch(::ACOPFFormulation) = true
allows_generation_redispatch(::ACOPFMinRedispatchFormulation) = true
allows_generation_redispatch(::ACOPFTapFormulation) = true
allows_generation_redispatch(::ACOPFMinVoltageViolationFormulation) = true
allows_generation_redispatch(::ACPFPmisVmisFormulation) = true

has_voltage_setpoint(::ACFormulation) = false
has_voltage_setpoint(::ACPFFeasibilityFormulation) = true

# whether branch tap ratios are JuMP decision variables (bounded by each branch's
# min_tap_ratio/max_tap_ratio) instead of fixed data
optimizes_tap(::ACFormulation) = false
optimizes_tap(::ACOPFTapFormulation) = true

# whether the voltage limits of some buses are relaxed into slacks that the objective drives to
# zero, instead of being hard bounds on V. See `voltage_violation_buses` for which buses those are.
minimizes_voltage_violation(::ACFormulation) = false
minimizes_voltage_violation(::ACOPFMinVoltageViolationFormulation) = true
