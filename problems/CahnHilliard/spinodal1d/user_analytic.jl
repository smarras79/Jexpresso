#---------------------------------------------------------------------------------
# Linear Cahn-Hilliard solution, overlaid on the 1D output by the plotter.
#
# Linearised about c = 0 (f'(c) ≈ -c), the mode grows as
#
#     c = δ e^{σt} cos(kx),     μ = (ε²k² - 1) c,     σ = k²(1 - ε²k²)
#
# This is the τ → 0 limit, so the gap includes the O(τσ) relaxation error.
# It stops being a solution once the cubic term matters: for amplitude A the
# c³ in μ is of relative size A² against the linear part. Past A = 0.05
# (0.25%) this returns NaN, which the plotter skips, and the late frames show
# the nonlinear phase separation alone.
#---------------------------------------------------------------------------------
function user_analytic_solution(x, t, outvar, inputs)

    p = ch_parameters()
    σ = ch_growth_rate(p.k)
    A = p.δ*exp(σ*t)

    nvar = length(outvar)
    q    = fill(NaN, length(x), nvar)
    A > 0.05 && return q

    for i in eachindex(x)
        c = A*cos(p.k*x[i])
        μ = (p.ε*p.ε*p.k*p.k - 1.0)*c

        for ivar = 1:nvar
            name = string(outvar[ivar])
            q[i, ivar] = name == "c" ? c :
                         name == "μ" ? μ : NaN
        end
    end

    return q
end
