abstract type ACFormulation end

struct ACOPFFormulation <: ACFormulation end
struct ACOPFMinRedispatchFormulation <: ACFormulation end
struct ACPFFeasibilityFormulation <: ACFormulation end
struct ACPFPmisVmisFormulation <: ACFormulation end

# standalone marker for the NonlinearSolve.jl path — not an ACFormulation, it never builds a JuMP model
struct ACPFFormulation end

allows_generation_redispatch(::ACFormulation) = false
allows_generation_redispatch(::ACOPFFormulation) = true
allows_generation_redispatch(::ACOPFMinRedispatchFormulation) = true
allows_generation_redispatch(::ACPFPmisVmisFormulation) = true

has_voltage_setpoint(::ACFormulation) = false
has_voltage_setpoint(::ACPFFeasibilityFormulation) = true
