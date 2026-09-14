# Verify that MOSEK is installed and licensed:  julia --project=. scripts/check_mosek.jl
using MomentMINLP, JuMP, LinearAlgebra
m = Model(mosek_optimizer()); set_silent(m)
@variable(m, X[1:3, 1:3], PSD)
C = [2.0 1 0; 1 2 1; 0 1 2]
@objective(m, Min, tr(C * X)); @constraint(m, tr(X) == 1)
optimize!(m)
println("MOSEK status: ", termination_status(m), "  obj = ", objective_value(m),
        "  (expected λmin = ", minimum(eigvals(C)), ")")
println("License file: ", get(ENV, "MOSEKLM_LICENSE_FILE", "~/mosek/mosek.lic (default)"))
