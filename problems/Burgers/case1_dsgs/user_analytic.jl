#---------------------------------------------------------------------------------
# Exact entropy solution of inviscid Burgers from q(x,0) = sin(2πx) on the
# periodic unit interval, for overlay on the 1D output plots (see the header
# of problems/CompEuler/sod1d/user_analytic.jl for how the plotter uses it).
#
# Along characteristics q is constant, so q(x,t) = sin(2πξ) with ξ the foot
#
#     x = ξ + t sin(2πξ).
#
# The data are odd about x = 0.5, so the shock that forms at t = 1/(2π) sits
# at x = 0.5 for all time, and q(1 - x) = -q(x). For x in (0, 0.5) the foot is
# on the branch of x(ξ) that still increases from ξ = 0: x'(ξ) = 1 + 2πt
# cos(2πξ) > 0 for ξ < ξ*, cos(2πξ*) = -1/(2πt) (ξ* = 0.5 before the shock
# forms). Feet beyond ξ* land past x = 0.5 and have entered the shock. x(ξ) is
# monotone on [0, ξ*], so bisection finds the foot.
#---------------------------------------------------------------------------------
function burgers_sine_exact(x, t)

    xr = mod(x, 1.0)
    xr == 0.5 && return 0.0
    xr > 0.5 && return -burgers_sine_exact(1.0 - xr, t)

    ξmax = 2π*t > 1.0 ? acos(-1.0/(2π*t))/(2π) : 0.5
    lo, hi = 0.0, ξmax
    for _ = 1:60
        ξ = 0.5*(lo + hi)
        ξ + t*sin(2π*ξ) < xr ? (lo = ξ) : (hi = ξ)
    end
    return sin(2π*0.5*(lo + hi))
end

function user_analytic_solution(x, t, outvar, inputs)

    q = fill(NaN, length(x), length(outvar))
    for i in eachindex(x)
        q[i, 1] = burgers_sine_exact(x[i], t)
    end
    return q
end
